#!/usr/bin/env python3
"""Check that files referenced by a samples YAML file exist."""

from __future__ import annotations

import argparse
import copy
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from workflow.scripts.sample_inputs import missing_sample_files, validate_samples


CRAM_SUFFIX = ".sorted.dups.cram"
STEM_TOKEN_RE = re.compile(r"[_-]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check fastq_1, fastq_2, cram, and crai paths in a samples YAML file."
    )
    parser.add_argument(
        "samples_yaml",
        nargs="?",
        default="config/samples.yaml",
        help="Samples YAML file to check (default: config/samples.yaml).",
    )
    parser.add_argument(
        "--fix-output",
        type=Path,
        default=None,
        help="Write accepted CRAM filename fixes to this new YAML file.",
    )
    return parser.parse_args()


def sample_name_from_cram(path: str) -> str | None:
    name = Path(path).name
    if not name.endswith(CRAM_SUFFIX):
        return None
    return name[: -len(CRAM_SUFFIX)]


def sample_tokens(sample_name: str) -> list[str]:
    return [token for token in STEM_TOKEN_RE.split(sample_name) if token]


def normalized_token(token: str) -> str:
    if token.isdigit():
        stripped = token.lstrip("0")
        return stripped or "0"
    return token


def sample_signature(sample_name: str) -> tuple[str, str] | None:
    tokens = sample_tokens(sample_name)
    if len(tokens) < 2:
        return None
    return normalized_token(tokens[0]), tokens[-1]


def candidate_search_roots(missing_path: Path) -> list[Path]:
    roots = [missing_path.parent]
    if len(missing_path.parents) >= 3:
        batch_root = missing_path.parents[2]
        if batch_root.is_dir() and batch_root not in roots:
            roots.append(batch_root)
    return roots


def rename_sample_references(samples: dict, old_sample: str, new_sample: str) -> None:
    tumours = samples.get("tumours", {})
    rebuilt_tumours = {}
    for tumour, normal in tumours.items():
        if tumour == old_sample:
            tumour = new_sample
        if normal == old_sample:
            normal = new_sample
        rebuilt_tumours[tumour] = normal
    if "tumours" in samples:
        samples["tumours"] = rebuilt_tumours
    comparisons = samples.get("comparisons", {})
    rebuilt_comparisons = {}
    for organoid, values in comparisons.items():
        organoid = new_sample if organoid == old_sample else organoid
        values = copy.deepcopy(values)
        if values.get("baseline") == old_sample:
            values["baseline"] = new_sample
        rebuilt_comparisons[organoid] = values
    if comparisons:
        samples["comparisons"] = rebuilt_comparisons


def find_normal_cram_alternatives(path: str) -> list[str]:
    return find_cram_alternatives(path)


def find_cram_alternatives(path: str) -> list[str]:
    missing_path = Path(path)
    missing_sample = sample_name_from_cram(path)
    missing_signature = sample_signature(missing_sample or "")
    if not missing_sample or not missing_signature:
        return []

    candidates = []
    missing_path_str = str(missing_path)
    for root in candidate_search_roots(missing_path):
        if root == missing_path.parent:
            candidates_iter = root.glob(f"*{CRAM_SUFFIX}")
        else:
            candidates_iter = root.glob(f"*/out/*{CRAM_SUFFIX}")
        for candidate in sorted(candidates_iter):
            if str(candidate) == missing_path_str:
                continue
            candidate_sample = sample_name_from_cram(str(candidate))
            if candidate_sample and sample_signature(candidate_sample) == missing_signature:
                candidates.append(str(candidate))

    return sorted(dict.fromkeys(candidates))


def prompt_for_candidate(sample: str, missing_path: str, candidates: list[str]) -> str | None:
    print("")
    print(f"Missing CRAM for sample {sample}:")
    print(f"  {missing_path}")

    if len(candidates) == 1:
        candidate = candidates[0]
        print("Suggested replacement:")
        print(f"  {candidate}")
        answer = input("Accept this replacement? [y/N] ").strip().lower()
        if answer in {"y", "yes"}:
            return candidate
        return None

    print("Suggested replacements:")
    for index, candidate in enumerate(candidates, start=1):
        print(f"  {index}. {candidate}")

    answer = input("Choose a replacement number, or press Enter to skip: ").strip()
    if not answer:
        return None
    try:
        index = int(answer)
    except ValueError:
        print("Skipping: response was not a number.")
        return None
    if index < 1 or index > len(candidates):
        print("Skipping: response was outside the candidate range.")
        return None
    return candidates[index - 1]


def apply_cram_replacement(samples: dict, old_sample: str, new_cram: str) -> bool:
    new_sample = sample_name_from_cram(new_cram)
    if not new_sample:
        print(f"WARNING: cannot derive sample name from {new_cram}; skipping.")
        return False

    sample_table = samples["samples"]
    if new_sample != old_sample and new_sample in sample_table:
        del sample_table[old_sample]
        rename_sample_references(samples, old_sample, new_sample)

        existing_cram = sample_table[new_sample].get("cram")
        default_crai = f"{existing_cram or new_cram}.crai"
        if not Path(default_crai).expanduser().exists():
            print(f"WARNING: default CRAI is still missing: {default_crai}")

        return True

    sample_data = copy.deepcopy(sample_table[old_sample])
    sample_data["cram"] = new_cram
    sample_data.pop("crai", None)

    if new_sample == old_sample:
        sample_table[old_sample] = sample_data
    else:
        rebuilt_samples = {}
        for sample, data in sample_table.items():
            if sample == old_sample:
                rebuilt_samples[new_sample] = sample_data
            else:
                rebuilt_samples[sample] = data
        samples["samples"] = rebuilt_samples

        rename_sample_references(samples, old_sample, new_sample)

    default_crai = f"{new_cram}.crai"
    if not Path(default_crai).expanduser().exists():
        print(f"WARNING: default CRAI is still missing: {default_crai}")

    return True


def repair_missing_crams(samples: dict, missing: list[tuple[str, str, str]]) -> bool:
    changed = False
    for sample, key, path in missing:
        if key != "cram":
            continue

        candidates = find_cram_alternatives(path)
        if not candidates:
            continue

        selected = prompt_for_candidate(sample, path, candidates)
        if not selected:
            continue

        if apply_cram_replacement(samples, sample, selected):
            changed = True

    return changed


def main() -> int:
    args = parse_args()
    samples_path = Path(args.samples_yaml)

    if not samples_path.exists():
        print(f"ERROR: sample config not found: {samples_path}", file=sys.stderr)
        return 2

    with samples_path.open() as handle:
        samples = yaml.safe_load(handle) or {}

    if "comparisons" in samples:
        try:
            validate_samples(samples, str(samples_path))
        except ValueError as error:
            print(f"ERROR: invalid sample manifest: {error}", file=sys.stderr)
            return 2

    missing = missing_sample_files(samples, str(samples_path))
    if missing and args.fix_output:
        if repair_missing_crams(samples, missing):
            args.fix_output.parent.mkdir(parents=True, exist_ok=True)
            with args.fix_output.open("w") as handle:
                yaml.safe_dump(samples, handle, sort_keys=False)
            print("")
            print(f"Wrote corrected sample config to {args.fix_output}.")
            missing = missing_sample_files(samples, str(args.fix_output))
        else:
            print("No fixes accepted; no output file written.")

    if missing:
        print("Missing sample files:")
        for sample, key, path in missing:
            print(f"  {sample}.{key}: {path}")
        return 1

    print(f"All sample files exist in {samples_path}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
