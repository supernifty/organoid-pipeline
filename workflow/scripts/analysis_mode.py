#!/usr/bin/env python3
"""Analysis-mode validation, territory preparation, and run safety helpers."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Iterable


VALID_ANALYSIS_TYPES = {"wes", "wgs"}
VALID_REFERENCE_BUILDS = {"grch37", "grch38"}

PRIMARY_CONTIG_LENGTHS = {
    "grch37": {
        "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
        "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
        "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
        "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
        "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    },
    "grch38": {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
        "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    },
}


def analysis_settings(config: dict) -> dict:
    analysis = config.get("analysis") or {}
    mode = str(analysis.get("type", "wes")).lower()
    if mode not in VALID_ANALYSIS_TYPES:
        raise ValueError("analysis.type must be one of: wes, wgs")
    wgs = analysis.get("wgs") or {}
    contigs = [str(value) for value in wgs.get(
        "contigs", [str(i) for i in range(1, 23)] + ["X", "Y"]
    )]
    if not contigs or len(contigs) != len(set(contigs)):
        raise ValueError("analysis.wgs.contigs must be a non-empty list without duplicates")
    target = int(wgs.get("target_bases_per_shard", 20_000_000))
    scatter = wgs.get("scatter_count")
    concurrency = int(wgs.get("max_concurrent_mutect2_shards", 32))
    if target <= 0 or (scatter is not None and int(scatter) <= 0) or concurrency <= 0:
        raise ValueError("WGS shard size, scatter count, and concurrency must be positive")
    return {
        "type": mode,
        "contigs": contigs,
        "target_bases_per_shard": target,
        "scatter_count": None if scatter is None else int(scatter),
        "max_concurrent_mutect2_shards": concurrency,
    }


def reference_settings(config: dict) -> dict:
    reference = config.get("reference") or {}
    build = str(reference.get("build", "grch37")).lower()
    if build not in VALID_REFERENCE_BUILDS:
        raise ValueError("reference.build must be one of: grch37, grch38")
    capture_source = str(reference.get("regions_source", "custom")).lower()
    if capture_source not in {"vendor", "remapped", "custom"}:
        raise ValueError("reference.regions_source must be one of: vendor, remapped, custom")
    source_build = str(reference.get("regions_source_build", build)).lower()
    if source_build not in VALID_REFERENCE_BUILDS:
        raise ValueError("reference.regions_source_build must be one of: grch37, grch38")
    if source_build != build:
        raise ValueError(
            "reference.regions must use the configured reference.build; remap the capture BED first"
        )
    if capture_source == "remapped" and not reference.get("regions_metadata"):
        raise ValueError("Remapped capture territory requires reference.regions_metadata")
    return {
        "build": build,
        "capture_source": capture_source,
        "capture_source_build": source_build,
        "capture_metadata": reference.get("regions_metadata"),
    }


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_fai(path: Path) -> tuple[list[str], dict[str, int]]:
    order, lengths = [], {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed FASTA index {path} line {line_number}")
            name, length = fields[0], int(fields[1])
            order.append(name)
            lengths[name] = length
    if not order:
        raise ValueError(f"Empty FASTA index: {path}")
    return order, lengths


def read_sequence_dictionary(path: Path) -> tuple[list[str], dict[str, int]]:
    order, lengths = [], {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith("@SQ"):
                continue
            values = dict(field.split(":", 1) for field in line.rstrip().split("\t")[1:])
            try:
                name, length = values["SN"], int(values["LN"])
            except (KeyError, ValueError) as error:
                raise ValueError(f"Malformed sequence dictionary {path} line {line_number}") from error
            order.append(name)
            lengths[name] = length
    if not order:
        raise ValueError(f"Empty sequence dictionary: {path}")
    return order, lengths


def read_vcf_dictionary(path: Path) -> dict[str, int | None]:
    dictionary: dict[str, int | None] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##contig=<"):
                body = line.rstrip()[10:-1]
                values = dict(
                    field.split("=", 1) for field in body.split(",") if "=" in field
                )
                if "ID" in values:
                    dictionary[values["ID"]] = (
                        int(values["length"]) if values.get("length", "").isdigit() else None
                    )
            elif line.startswith("#CHROM"):
                break
    if not dictionary:
        raise ValueError(f"VCF has no sequence dictionary: {path}")
    return dictionary


def read_interval_list(path: Path) -> tuple[dict[str, int], list[tuple[str, int, int]]]:
    dictionary, intervals = {}, []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith("@SQ"):
                values = dict(field.split(":", 1) for field in line.rstrip().split("\t")[1:])
                dictionary[values["SN"]] = int(values["LN"])
            elif line.startswith("@"):
                continue
            else:
                fields = line.rstrip().split("\t")
                if len(fields) < 3:
                    raise ValueError(f"Malformed interval list {path} line {line_number}")
                intervals.append((fields[0], int(fields[1]) - 1, int(fields[2])))
    return dictionary, intervals


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    intervals = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED {path} line {line_number}")
            intervals.append((fields[0], int(fields[1]), int(fields[2])))
    return intervals


def validate_intervals(
    intervals: Iterable[tuple[str, int, int]], order: list[str], lengths: dict[str, int]
) -> list[tuple[str, int, int]]:
    rank = {name: index for index, name in enumerate(order)}
    result, previous = [], None
    for interval in intervals:
        contig, start, end = interval
        if contig not in rank:
            raise ValueError(f"Territory contig {contig!r} is absent from the reference index")
        if start < 0 or end <= start or end > lengths[contig]:
            raise ValueError(f"Invalid territory interval {contig}:{start + 1}-{end}")
        key = (rank[contig], start, end)
        if previous is not None and key < previous:
            raise ValueError("Territory intervals are not in reference dictionary order")
        if result and result[-1][0] == contig and start < result[-1][2]:
            raise ValueError(f"Overlapping territory intervals on {contig}")
        result.append(interval)
        previous = key
    if not result:
        raise ValueError("Analysis territory is empty")
    return result


def merge_adjacent(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for contig, start, end in intervals:
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            old = merged[-1]
            merged[-1] = (contig, old[1], max(old[2], end))
        else:
            merged.append((contig, start, end))
    return merged


def subtract_intervals(
    territory: list[tuple[str, int, int]], exclusions: list[tuple[str, int, int]]
) -> list[tuple[str, int, int]]:
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, start, end in exclusions:
        by_contig.setdefault(contig, []).append((start, end))
    output = []
    for contig, start, end in territory:
        pieces = [(start, end)]
        for exc_start, exc_end in by_contig.get(contig, []):
            updated = []
            for piece_start, piece_end in pieces:
                if exc_end <= piece_start or exc_start >= piece_end:
                    updated.append((piece_start, piece_end))
                else:
                    if piece_start < exc_start:
                        updated.append((piece_start, exc_start))
                    if exc_end < piece_end:
                        updated.append((exc_end, piece_end))
            pieces = updated
        output.extend((contig, piece_start, piece_end) for piece_start, piece_end in pieces)
    return output


def prepare_wgs_territory(
    source: Path, fai: Path, contigs: list[str], exclude: Path | None,
    bed_output: Path, interval_output: Path, metadata_output: Path,
) -> None:
    reference_order, reference_lengths = read_fai(fai)
    missing = [contig for contig in contigs if contig not in reference_lengths]
    if missing:
        raise ValueError(f"WGS contigs absent from reference index: {', '.join(missing)}")
    source_dictionary, intervals = read_interval_list(source)
    for contig in contigs:
        if contig not in source_dictionary:
            raise ValueError(f"WGS calling interval dictionary lacks contig {contig}")
        if source_dictionary[contig] != reference_lengths[contig]:
            raise ValueError(
                f"Reference length mismatch for {contig}: interval list "
                f"{source_dictionary[contig]}, FASTA {reference_lengths[contig]}"
            )
    selected = [item for item in intervals if item[0] in set(contigs)]
    selected = validate_intervals(selected, reference_order, reference_lengths)
    if exclude:
        exclusions = merge_adjacent(sorted(
            read_bed(exclude), key=lambda item: (reference_order.index(item[0]), item[1], item[2])
        ))
        validate_intervals(exclusions, reference_order, reference_lengths)
        selected = subtract_intervals(selected, exclusions)
    selected = validate_intervals(selected, reference_order, reference_lengths)
    bed_output.parent.mkdir(parents=True, exist_ok=True)
    with bed_output.open("w") as bed:
        for contig, start, end in selected:
            bed.write(f"{contig}\t{start}\t{end}\n")
    with interval_output.open("w") as out:
        out.write("@HD\tVN:1.6\tSO:coordinate\n")
        for contig in reference_order:
            out.write(f"@SQ\tSN:{contig}\tLN:{reference_lengths[contig]}\n")
        for contig, start, end in selected:
            out.write(f"{contig}\t{start + 1}\t{end}\t+\tCALLABLE\n")
    bases = sum(end - start for _, start, end in selected)
    metadata_output.write_text(json.dumps({
        "callable_bases": bases,
        "contigs": contigs,
        "interval_count": len(selected),
    }, indent=2, sort_keys=True) + "\n")


def validate_reference_profile(
    build: str,
    reference_fai: Path,
    reference_dict: Path,
    territory: Path,
    chromosomes: list[str],
    vcfs: list[Path],
    output: Path,
    capture_source: str = "custom",
    capture_metadata: Path | None = None,
) -> None:
    if build not in VALID_REFERENCE_BUILDS:
        raise ValueError(f"Unsupported reference build: {build}")
    fai_order, fai_lengths = read_fai(reference_fai)
    dict_order, dict_lengths = read_sequence_dictionary(reference_dict)
    if fai_order != dict_order or fai_lengths != dict_lengths:
        raise ValueError("Reference FASTA index and sequence dictionary do not match")

    expected = PRIMARY_CONTIG_LENGTHS[build]
    expected_contigs = list(expected)
    if chromosomes != expected_contigs:
        raise ValueError(
            f"Configured chromosomes do not match {build}: expected {', '.join(expected_contigs)}"
        )
    for contig, expected_length in expected.items():
        if fai_lengths.get(contig) != expected_length:
            raise ValueError(
                f"Reference does not match {build}: {contig} length is "
                f"{fai_lengths.get(contig)!r}, expected {expected_length}"
            )

    rank = {name: index for index, name in enumerate(fai_order)}
    intervals = read_bed(territory)
    intervals = validate_intervals(intervals, fai_order, fai_lengths)
    territory_contigs = {contig for contig, _, _ in intervals}
    unexpected = sorted(territory_contigs - set(expected), key=lambda name: rank[name])
    if unexpected:
        raise ValueError(
            "Analysis territory contains non-primary contigs: " + ", ".join(unexpected)
        )

    for vcf in vcfs:
        vcf_dictionary = read_vcf_dictionary(vcf)
        missing = sorted(territory_contigs - set(vcf_dictionary), key=lambda name: rank[name])
        if missing:
            raise ValueError(f"VCF {vcf} lacks territory contigs: {', '.join(missing)}")
        for contig in territory_contigs:
            vcf_length = vcf_dictionary[contig]
            if vcf_length is not None and vcf_length != fai_lengths[contig]:
                raise ValueError(
                    f"VCF/reference length mismatch for {contig} in {vcf}: "
                    f"{vcf_length} != {fai_lengths[contig]}"
                )

    if capture_source == "remapped":
        if capture_metadata is None:
            raise ValueError("Remapped capture territory requires capture metadata")
        metadata = json.loads(capture_metadata.read_text())
        if metadata.get("target_build") != build:
            raise ValueError(
                f"Capture metadata target build is {metadata.get('target_build')!r}, expected {build!r}"
            )
        if metadata.get("output_sha256") != sha256(territory):
            raise ValueError("Capture metadata checksum does not match reference.regions")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({
        "reference_build": build,
        "reference_contig_count": len(fai_order),
        "territory_bases": sum(end - start for _, start, end in intervals),
        "territory_contigs": sorted(territory_contigs, key=lambda name: rank[name]),
        "territory_interval_count": len(intervals),
        "validated_vcfs": [str(path) for path in vcfs],
        "capture_source": capture_source,
        "capture_metadata": str(capture_metadata) if capture_metadata else None,
    }, indent=2, sort_keys=True) + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_payload(
    mode: str, reference_fai: Path, territory: Path, reference_dict: Path | None = None,
    reference_build: str = "grch37", capture_source: str = "custom",
    capture_metadata: Path | None = None,
) -> dict:
    intervals = read_bed(territory)
    return {
        "schema_version": 2,
        "analysis_type": mode,
        "reference_build": reference_build,
        "capture_source": capture_source,
        "capture_metadata_sha256": sha256(capture_metadata) if capture_metadata else None,
        "reference_fai_sha256": sha256(reference_fai),
        "reference_dict_sha256": sha256(reference_dict) if reference_dict else None,
        "territory_sha256": sha256(territory),
        "territory_bases": sum(end - start for _, start, end in intervals),
    }


def pon_fingerprint_payload(
    mode: str, reference_fai: Path, territory: Path, normals: list[str],
    reference_dict: Path | None = None, samples: Path | None = None,
    reference_build: str = "grch37",
) -> dict:
    payload = manifest_payload(
        mode, reference_fai, territory, reference_dict, reference_build=reference_build
    )
    payload["normal_samples"] = sorted(normals)
    payload["samples_config_sha256"] = sha256(samples) if samples else None
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return {"fingerprint": hashlib.sha256(encoded).hexdigest(), **payload}


def validate_existing_manifest(
    results_dir: Path, requested_mode: str, requested_build: str | None = None
) -> None:
    manifest = results_dir / "analysis_manifest.json"
    if manifest.exists():
        current = json.loads(manifest.read_text()).get("analysis_type")
        if current != requested_mode:
            raise ValueError(
                f"results/ contains {current!r} outputs but analysis.type is {requested_mode!r}; "
                "use a clean output directory"
            )
        current_build = json.loads(manifest.read_text()).get("reference_build")
        if requested_build and current_build and current_build != requested_build:
            raise ValueError(
                f"results/ contains {current_build!r} outputs but reference.build is "
                f"{requested_build!r}; use a clean output directory"
            )


def shard_count(callable_bases: int, target: int, configured: int | None) -> int:
    if callable_bases <= 0:
        raise ValueError("Callable territory must contain at least one base")
    return configured if configured is not None else math.ceil(callable_bases / target)


def validate_contamination(pileup: Path, contamination: Path) -> None:
    rows = [line for line in pileup.read_text().splitlines() if line and not line.startswith("#")]
    if len(rows) <= 1:
        raise ValueError(f"No usable contamination sites were produced in {pileup}")
    values = []
    for line in contamination.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > 1 and fields[0].lower() == "sample" and fields[1].lower() == "contamination":
            continue
        try:
            values.append(float(fields[1]))
        except (IndexError, ValueError) as error:
            raise ValueError(f"Malformed contamination table row: {line}") from error
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("CalculateContamination did not produce a finite estimate")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    territory = sub.add_parser("prepare-territory")
    territory.add_argument("--source", required=True, type=Path)
    territory.add_argument("--fai", required=True, type=Path)
    territory.add_argument("--contig", action="append", required=True)
    territory.add_argument("--exclude", type=Path)
    territory.add_argument("--bed-output", required=True, type=Path)
    territory.add_argument("--interval-output", required=True, type=Path)
    territory.add_argument("--metadata-output", required=True, type=Path)
    validate = sub.add_parser("validate-contamination")
    validate.add_argument("--pileup", required=True, type=Path)
    validate.add_argument("--contamination", required=True, type=Path)
    manifest = sub.add_parser("write-manifest")
    manifest.add_argument("--mode", choices=sorted(VALID_ANALYSIS_TYPES), required=True)
    manifest.add_argument("--build", choices=sorted(VALID_REFERENCE_BUILDS), required=True)
    manifest.add_argument("--capture-source", choices=["vendor", "remapped", "custom"], required=True)
    manifest.add_argument("--capture-metadata", type=Path)
    manifest.add_argument("--reference-fai", required=True, type=Path)
    manifest.add_argument("--reference-dict", required=True, type=Path)
    manifest.add_argument("--territory", required=True, type=Path)
    manifest.add_argument("--output", required=True, type=Path)
    fingerprint = sub.add_parser("write-pon-fingerprint")
    fingerprint.add_argument("--mode", choices=sorted(VALID_ANALYSIS_TYPES), required=True)
    fingerprint.add_argument("--build", choices=sorted(VALID_REFERENCE_BUILDS), required=True)
    fingerprint.add_argument("--reference-fai", required=True, type=Path)
    fingerprint.add_argument("--reference-dict", required=True, type=Path)
    fingerprint.add_argument("--territory", required=True, type=Path)
    fingerprint.add_argument("--samples", required=True, type=Path)
    fingerprint.add_argument("--normal", action="append", default=[])
    fingerprint.add_argument("--output", required=True, type=Path)
    profile = sub.add_parser("validate-reference-profile")
    profile.add_argument("--build", choices=sorted(VALID_REFERENCE_BUILDS), required=True)
    profile.add_argument("--reference-fai", required=True, type=Path)
    profile.add_argument("--reference-dict", required=True, type=Path)
    profile.add_argument("--territory", required=True, type=Path)
    profile.add_argument("--chromosome", action="append", required=True)
    profile.add_argument("--vcf", action="append", default=[], type=Path)
    profile.add_argument("--capture-source", choices=["vendor", "remapped", "custom"], required=True)
    profile.add_argument("--capture-metadata", type=Path)
    profile.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare-territory":
        prepare_wgs_territory(args.source, args.fai, args.contig, args.exclude,
                              args.bed_output, args.interval_output, args.metadata_output)
    elif args.command == "validate-contamination":
        validate_contamination(args.pileup, args.contamination)
    elif args.command == "write-manifest":
        payload = manifest_payload(
            args.mode, args.reference_fai, args.territory, args.reference_dict,
            args.build, args.capture_source, args.capture_metadata,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, args.output)
    elif args.command == "write-pon-fingerprint":
        payload = pon_fingerprint_payload(
            args.mode, args.reference_fai, args.territory, args.normal,
            args.reference_dict, args.samples, args.build,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, args.output)
    else:
        validate_reference_profile(
            args.build, args.reference_fai, args.reference_dict, args.territory,
            args.chromosome, args.vcf, args.output, args.capture_source, args.capture_metadata,
        )


if __name__ == "__main__":
    main()
