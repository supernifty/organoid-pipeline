#!/usr/bin/env python3
"""Validate and atomically import legacy duplicate-marked BAMs as CRAM 3.0."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import uuid
from pathlib import Path
from typing import Any

import yaml


SCHEMA_VERSION = 1
GRCH37_PRIMARY_LENGTHS = {
    **{str(index): length for index, length in enumerate((
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
        159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
        115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566,
    ), 1)},
    "X": 155270560,
    "Y": 59373566,
}


class ImportError(RuntimeError):
    pass


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def identity(path: Path, checksum: bool = True) -> dict[str, Any]:
    stat = path.stat()
    value = {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        value["sha256"] = sha256(path)
    return value


def run(command: list[str], commands: list[list[str]], capture: bool = False) -> str:
    commands.append(command)
    completed = subprocess.run(command, check=False, text=True, capture_output=capture)
    if completed.returncode:
        detail = completed.stderr.strip() if capture else ""
        raise ImportError(f"command failed ({completed.returncode}): {' '.join(command)}{': ' + detail if detail else ''}")
    return completed.stdout if capture else ""


def load_legacy_samples(path: Path) -> tuple[list[str], dict[str, str]]:
    value = yaml.safe_load(path.read_text()) or {}
    if not isinstance(value, dict):
        raise ImportError("legacy sample YAML root must be a mapping")
    raw_samples = value.get("samples")
    if isinstance(raw_samples, dict):
        sample_ids = list(raw_samples)
    elif isinstance(raw_samples, list):
        sample_ids = []
        for entry in raw_samples:
            if isinstance(entry, str):
                sample_ids.append(entry)
            elif isinstance(entry, dict) and isinstance(entry.get("sample"), str):
                sample_ids.append(entry["sample"])
            elif isinstance(entry, dict) and len(entry) == 1:
                sample_ids.append(str(next(iter(entry))))
            else:
                raise ImportError("list-style samples must contain IDs or mappings with a sample key")
    else:
        raise ImportError("legacy sample YAML must define samples as a mapping or list")
    if not sample_ids or len(sample_ids) != len(set(sample_ids)):
        raise ImportError("sample IDs must be non-empty and unique")

    raw_pairs = value.get("tumours", {})
    if isinstance(raw_pairs, dict):
        pairs = {str(tumour): str(normal) for tumour, normal in raw_pairs.items()}
    elif isinstance(raw_pairs, list):
        pairs = {}
        for entry in raw_pairs:
            if not isinstance(entry, dict):
                raise ImportError("list-style tumour pairs must be mappings")
            tumour = entry.get("tumour")
            normal = entry.get("normal")
            if not isinstance(tumour, str) or not isinstance(normal, str) or tumour in pairs:
                raise ImportError("each list-style pair requires unique tumour and normal IDs")
            pairs[tumour] = normal
    else:
        raise ImportError("tumours must be a mapping or list")
    unknown = sorted((set(pairs) | set(pairs.values())) - set(sample_ids))
    if unknown:
        raise ImportError("pair mapping references unknown samples: " + ", ".join(unknown))
    if any(tumour == normal for tumour, normal in pairs.items()):
        raise ImportError("a tumour cannot be paired with itself")
    return sample_ids, pairs


def resolve_bam_index(bam: Path) -> Path:
    candidates = (Path(f"{bam}.bai"), bam.with_suffix(".bai"))
    present = [path for path in candidates if path.is_file()]
    if not present:
        raise ImportError(f"missing BAM index for {bam}; expected {candidates[0]} or {candidates[1]}")
    if len(present) == 2 and (present[0].stat().st_size != present[1].stat().st_size or sha256(present[0]) != sha256(present[1])):
        raise ImportError(f"conflicting BAM indexes for {bam}")
    return present[0]


def parse_header(text: str) -> dict[str, Any]:
    sort_orders, sequences, read_groups = set(), {}, []
    for line in text.splitlines():
        fields = line.split("\t")
        tags = dict(field.split(":", 1) for field in fields[1:] if ":" in field)
        if fields[0] == "@HD" and "SO" in tags:
            sort_orders.add(tags["SO"])
        elif fields[0] == "@SQ":
            if "SN" not in tags or "LN" not in tags or tags["SN"] in sequences:
                raise ImportError("invalid or duplicate @SQ header record")
            sequences[tags["SN"]] = int(tags["LN"])
        elif fields[0] == "@RG":
            read_groups.append(tags)
    if sort_orders != {"coordinate"}:
        raise ImportError("alignment header must declare SO:coordinate")
    if not read_groups or any(not group.get("ID") or not group.get("SM") for group in read_groups):
        raise ImportError("every input must contain read groups with ID and SM tags")
    sample_names = {group["SM"] for group in read_groups}
    if len(sample_names) != 1:
        raise ImportError("all read groups in one input must have one unique SM value")
    return {"sequences": sequences, "read_groups": read_groups, "sample_name": next(iter(sample_names))}


def validate_primary_contigs(sequences: dict[str, int]) -> None:
    missing = [name for name in GRCH37_PRIMARY_LENGTHS if name not in sequences]
    prefixed = [name for name in sequences if name.startswith("chr") and name[3:] in GRCH37_PRIMARY_LENGTHS]
    wrong = [name for name, length in GRCH37_PRIMARY_LENGTHS.items() if sequences.get(name) not in (None, length)]
    if missing or prefixed or wrong:
        details = []
        if missing:
            details.append("missing " + ",".join(missing))
        if prefixed:
            details.append("UCSC-style names " + ",".join(prefixed))
        if wrong:
            details.append("length mismatch " + ",".join(wrong))
        raise ImportError("input is not GRCh37 primary-contig compatible: " + "; ".join(details))


def reference_lengths(reference: Path) -> dict[str, int]:
    result = {}
    with Path(f"{reference}.fai").open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            if len(fields) < 2 or fields[0] in result:
                raise ImportError(f"malformed reference index line {number}")
            result[fields[0]] = int(fields[1])
    return result


def idxstats_totals(text: str) -> dict[str, int]:
    mapped = unmapped = 0
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            raise ImportError("malformed samtools idxstats output")
        mapped += int(fields[2])
        unmapped += int(fields[3])
    return {"mapped": mapped, "unmapped": unmapped, "total": mapped + unmapped}


def alignment_checksum(text: str) -> list[str]:
    rows = [line for line in text.splitlines() if line and not line.startswith("#")]
    if not rows:
        raise ImportError("samtools checksum produced no alignment checksum rows")
    return rows


def validate_alignment(path: Path, reference: Path, commands: list[list[str]], require_grch37: bool) -> dict[str, Any]:
    run(["samtools", "quickcheck", "-v", str(path)], commands)
    header = parse_header(run(["samtools", "view", "-H", "-T", str(reference), str(path)], commands, True))
    if require_grch37:
        validate_primary_contigs(header["sequences"])
    lengths = reference_lengths(reference)
    mismatched = [name for name, length in header["sequences"].items() if lengths.get(name) != length]
    if mismatched:
        raise ImportError("alignment header does not match the reference index: " + ", ".join(mismatched))
    totals = idxstats_totals(run(["samtools", "idxstats", "-@", "1", "--input-fmt-option", f"reference={reference}", str(path)], commands, True))
    checksum = alignment_checksum(run([
        "samtools", "checksum", "-a", "-T", "--input-fmt-option", f"reference={reference}", str(path)
    ], commands, True))
    return {"header": header, "totals": totals, "alignment_checksum": checksum, "quickcheck": "pass"}


def comparable_header(header: dict[str, Any]) -> dict[str, Any]:
    return {"sequences": header["sequences"], "read_groups": header["read_groups"], "sample_name": header["sample_name"]}


def validate_output(source: dict[str, Any], output: dict[str, Any]) -> None:
    if comparable_header(source["header"]) != comparable_header(output["header"]):
        raise ImportError("source and CRAM headers differ")
    if source["totals"] != output["totals"]:
        raise ImportError("source and CRAM indexed mapped/unmapped totals differ")
    if source["alignment_checksum"] != output["alignment_checksum"]:
        raise ImportError("source and CRAM format-independent alignment checksums differ")


def destination_paths(destination: Path, sample: str) -> tuple[Path, Path]:
    cram = destination / f"{sample}.sorted.dups.cram"
    return cram, Path(f"{cram}.crai")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--samples", required=True, type=Path, help="Legacy sample YAML")
    result.add_argument("--legacy-output", required=True, type=Path, help="Directory containing duplicate-marked BAMs")
    result.add_argument("--destination", required=True, type=Path, help="CRAM destination directory")
    result.add_argument("--reference", required=True, type=Path, help="GRCh37 FASTA")
    result.add_argument("--output-samples", required=True, type=Path, help="Modern sample YAML to publish")
    result.add_argument("--manifest", type=Path, help="Manifest path (default: DESTINATION/import_manifest.json)")
    result.add_argument("--threads", type=int, default=8)
    result.add_argument("--preflight-only", action="store_true")
    result.add_argument("--allow-noncanonical-reference", action="store_true", help=argparse.SUPPRESS)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.threads < 1:
        raise ImportError("--threads must be at least 1")
    sample_ids, pairs = load_legacy_samples(args.samples)
    reference = args.reference.resolve()
    if not reference.is_file() or not Path(f"{reference}.fai").is_file():
        raise ImportError("reference FASTA and .fai must exist")
    if not args.allow_noncanonical_reference:
        validate_primary_contigs(reference_lengths(reference))
    args.destination.mkdir(parents=True, exist_ok=True)
    manifest_path = args.manifest or args.destination / "import_manifest.json"

    sources = []
    required = 0
    for sample in sample_ids:
        bam = args.legacy_output / f"{sample}.sorted.dups.bam"
        if not bam.is_file():
            raise ImportError(f"missing legacy BAM: {bam}")
        index = resolve_bam_index(bam)
        required += bam.stat().st_size + index.stat().st_size
        cram, crai = destination_paths(args.destination, sample)
        if cram.exists() != crai.exists():
            raise ImportError(f"partial destination for {sample}; CRAM and CRAI must both exist or both be absent")
        sources.append((sample, bam, index, cram, crai))
    if shutil.disk_usage(args.destination).free < required:
        raise ImportError(f"insufficient destination capacity: need at least {required} bytes")
    if args.preflight_only:
        print(json.dumps({"samples": len(sources), "source_bytes": required, "free_bytes": shutil.disk_usage(args.destination).free}, sort_keys=True))
        return 0

    started = utcnow()
    reference_identity = {"fasta": identity(reference), "fai": identity(Path(f"{reference}.fai"))}
    records, sample_names = [], {}
    modern = {"samples": {}, "tumours": pairs}
    for sample, bam, index, cram, crai in sources:
        commands: list[list[str]] = []
        source_validation = validate_alignment(bam, reference, commands, not args.allow_noncanonical_reference)
        sm = source_validation["header"]["sample_name"]
        if sm in sample_names:
            raise ImportError(f"duplicate SM value {sm!r} in samples {sample_names[sm]!r} and {sample!r}")
        sample_names[sm] = sample
        reused = cram.exists()
        if reused:
            output_validation = validate_alignment(cram, reference, commands, not args.allow_noncanonical_reference)
            validate_output(source_validation, output_validation)
        else:
            staging = args.destination / f".{sample}.{uuid.uuid4().hex}.staging"
            staging.mkdir()
            staged_cram = staging / cram.name
            staged_crai = Path(f"{staged_cram}.crai")
            try:
                run(["samtools", "view", "-@", str(args.threads), "-C", "-T", str(reference),
                     "--output-fmt-option", "version=3.0", "-o", str(staged_cram), str(bam)], commands)
                run(["samtools", "index", "-@", str(args.threads), str(staged_cram), str(staged_crai)], commands)
                output_validation = validate_alignment(staged_cram, reference, commands, not args.allow_noncanonical_reference)
                validate_output(source_validation, output_validation)
                os.replace(staged_cram, cram)
                os.replace(staged_crai, crai)
            finally:
                if staging.exists():
                    shutil.rmtree(staging)
        values = {"cram": str(cram.resolve()), "crai": str(crai.resolve())}
        if sm != sample:
            values["bam_sample"] = sm
        modern["samples"][sample] = values
        records.append({
            "sample": sample, "bam_sample": sm, "pair_role": "tumour" if sample in pairs else "normal",
            "source": identity(bam), "source_index": identity(index), "output": identity(cram),
            "output_index": identity(crai), "reference": reference_identity, "commands": commands,
            "validation": {"source": source_validation, "output": output_validation, "comparison": "pass"},
            "reused": reused,
        })

    args.output_samples.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_yaml = args.output_samples.with_name(f".{args.output_samples.name}.{uuid.uuid4().hex}.tmp")
    temp_manifest = manifest_path.with_name(f".{manifest_path.name}.{uuid.uuid4().hex}.tmp")
    temp_yaml.write_text(yaml.safe_dump(modern, sort_keys=False))
    temp_manifest.write_text(json.dumps({
        "schema_version": SCHEMA_VERSION, "started_at": started, "completed_at": utcnow(),
        "legacy_samples": identity(args.samples), "legacy_output": str(args.legacy_output.resolve()),
        "destination": str(args.destination.resolve()), "reference": reference_identity,
        "tumour_normal_pairs": pairs, "samples": records,
    }, indent=2, sort_keys=True) + "\n")
    os.replace(temp_yaml, args.output_samples)
    os.replace(temp_manifest, manifest_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ImportError as exc:
        raise SystemExit(f"ERROR: {exc}")
