#!/usr/bin/env python3
"""Create expected Somalier groups and stable cohort mismatch/swap flags."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import yaml


def validate_somalier_config(config: dict, reference_build: str) -> None:
    somalier = config.get("somalier", {})
    if not somalier.get("enabled", False):
        return
    if reference_build != "grch37":
        raise ValueError("the configured Somalier workflow currently requires GRCh37")
    for key in ("sites_vcf", "sites_vcf_index"):
        path = somalier.get(key)
        if not path or not Path(path).is_file():
            raise ValueError(f"somalier.{key} must name an existing file")
    minimum = int(somalier.get("minimum_depth", 20))
    expected = float(somalier.get("expected_pair_concordance_min", 0.6))
    unexpected = float(somalier.get("unexpected_pair_concordance_max", 0.4))
    if minimum < 1 or not 0 <= unexpected < expected <= 1:
        raise ValueError("invalid Somalier depth or concordance thresholds")
    ancestry = somalier.get("ancestry", {})
    if ancestry.get("enabled", False):
        if not ancestry.get("labels") or not Path(ancestry["labels"]).is_file():
            raise ValueError("Somalier ancestry labels must exist")
        directory = Path(ancestry.get("reference_somalier_dir", ""))
        if not directory.is_dir() or not any(directory.glob("*.somalier")):
            raise ValueError("Somalier ancestry reference directory must contain .somalier files")


def expected_pairs(samples: dict) -> set[frozenset[str]]:
    return {frozenset((str(tumour), str(normal))) for tumour, normal in samples.get("tumours", {}).items()}


def write_groups(samples_path: Path, output: Path) -> None:
    samples = yaml.safe_load(samples_path.read_text()) or {}
    groups = {}
    for tumour, normal in samples.get("tumours", {}).items():
        groups.setdefault(str(normal), set()).add(str(tumour))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for normal in sorted(groups):
            handle.write(",".join([normal, *sorted(groups[normal])]) + "\n")


def write_flags(samples_path: Path, pairs_path: Path, output: Path, expected_min: float, unexpected_max: float) -> None:
    if not 0 <= unexpected_max < expected_min <= 1:
        raise ValueError("Somalier thresholds must satisfy 0 ≤ unexpected maximum < expected minimum ≤ 1")
    samples = yaml.safe_load(samples_path.read_text()) or {}
    expected = expected_pairs(samples)
    rows = []
    with pairs_path.open(newline="") as handle:
        first = handle.readline().lstrip("#")
        reader = csv.DictReader([first, *handle], delimiter="\t")
        required = {"sample_a", "sample_b", "relatedness"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Somalier pairs table lacks sample_a, sample_b, or relatedness")
        for row in reader:
            pair = frozenset((row["sample_a"], row["sample_b"]))
            relatedness = float(row["relatedness"])
            is_expected = pair in expected
            flag = "pass"
            if is_expected and relatedness < expected_min:
                flag = "expected_pair_mismatch"
            elif not is_expected and relatedness >= unexpected_max:
                flag = "unexpected_pair_possible_swap"
            rows.append({"sample_a": row["sample_a"], "sample_b": row["sample_b"],
                         "expected_pair": str(is_expected).lower(), "relatedness": relatedness, "flag": flag})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("sample_a", "sample_b", "expected_pair", "relatedness", "flag"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["sample_a"], row["sample_b"])))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    groups = sub.add_parser("groups")
    groups.add_argument("--samples", required=True, type=Path)
    groups.add_argument("--output", required=True, type=Path)
    flags = sub.add_parser("flags")
    flags.add_argument("--samples", required=True, type=Path)
    flags.add_argument("--pairs", required=True, type=Path)
    flags.add_argument("--expected-min", required=True, type=float)
    flags.add_argument("--unexpected-max", required=True, type=float)
    flags.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "groups":
        write_groups(args.samples, args.output)
    else:
        write_flags(args.samples, args.pairs, args.output, args.expected_min, args.unexpected_max)


if __name__ == "__main__":
    main()
