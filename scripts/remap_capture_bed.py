#!/usr/bin/env python3
"""Remap an hg19 capture BED to GRCh38 and publish a validated canonical BED."""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fai(path: Path) -> tuple[list[str], dict[str, int]]:
    order, lengths = [], {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            order.append(fields[0])
            lengths[fields[0]] = int(fields[1])
    return order, lengths


def source_contig(contig: str) -> str:
    if contig.startswith("chr"):
        return contig
    if contig in {"M", "MT"}:
        return "chrM"
    return f"chr{contig}"


def read_source(path: Path) -> list[tuple[str, int, int, str]]:
    intervals = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED line {line_number}: {line.rstrip()}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval on line {line_number}")
            intervals.append((source_contig(fields[0]), start, end, f"capture_{len(intervals):08d}"))
    if not intervals:
        raise ValueError("Capture BED is empty")
    return intervals


def merge_intervals(
    intervals: list[tuple[str, int, int]], order: list[str]
) -> list[tuple[str, int, int]]:
    rank = {contig: index for index, contig in enumerate(order)}
    intervals.sort(key=lambda value: (rank[value[0]], value[1], value[2]))
    merged: list[tuple[str, int, int]] = []
    for contig, start, end in intervals:
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1] = (contig, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((contig, start, end))
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-bed", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--target-fai", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="Canonical .bed.gz output")
    parser.add_argument("--mapping", required=True, type=Path, help="Raw mapped BED with source IDs")
    parser.add_argument("--unmapped", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--minimum-mapped-base-fraction", type=float, default=0.99)
    parser.add_argument("--liftover", default="liftOver")
    parser.add_argument("--bgzip", default="bgzip")
    parser.add_argument("--tabix", default="tabix")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output.suffixes[-2:] != [".bed", ".gz"]:
        raise ValueError("--output must end in .bed.gz")
    if not 0 < args.minimum_mapped_base_fraction <= 1:
        raise ValueError("--minimum-mapped-base-fraction must be in (0, 1]")
    for executable in (args.liftover, args.bgzip, args.tabix):
        if not shutil.which(executable):
            raise ValueError(f"Required executable is unavailable: {executable}")

    source = read_source(args.source_bed)
    source_bases = sum(end - start for _, start, end, _ in source)
    source_ids = {identifier for *_, identifier in source}
    source_lengths = {identifier: end - start for _, start, end, identifier in source}
    target_order, target_lengths = read_fai(args.target_fai)
    allowed = {f"chr{value}" for value in range(1, 23)} | {"chrX", "chrY"}

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.mapping.parent.mkdir(parents=True, exist_ok=True)
    args.unmapped.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="capture-remap.", dir=args.output.parent) as directory:
        work = Path(directory)
        prepared = work / "source.bed"
        mapped = work / "mapped.bed"
        raw_unmapped = work / "unmapped.bed"
        with prepared.open("w") as handle:
            for interval in source:
                handle.write("\t".join(map(str, interval)) + "\n")

        subprocess.run(
            [args.liftover, str(prepared), str(args.chain), str(mapped), str(raw_unmapped)],
            check=True,
        )

        mapped_intervals: list[tuple[str, int, int]] = []
        mapped_ids: list[str] = []
        with mapped.open() as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip().split("\t")
                if len(fields) < 4:
                    raise ValueError(f"Malformed liftOver output line {line_number}")
                contig, start, end, identifier = fields[0], int(fields[1]), int(fields[2]), fields[3]
                if identifier not in source_ids:
                    raise ValueError(f"Unexpected interval identifier from liftOver: {identifier}")
                if contig not in allowed:
                    raise ValueError(f"Mapped capture interval is not on a primary contig: {contig}")
                if contig not in target_lengths or start < 0 or end <= start or end > target_lengths[contig]:
                    raise ValueError(f"Mapped interval is incompatible with target reference: {line.rstrip()}")
                mapped_intervals.append((contig, start, end))
                mapped_ids.append(identifier)

        mapped_source_bases = sum(source_lengths[identifier] for identifier in set(mapped_ids))
        mapped_base_fraction = mapped_source_bases / source_bases
        if mapped_base_fraction < args.minimum_mapped_base_fraction:
            raise ValueError(
                f"Only {mapped_base_fraction:.4%} of capture bases mapped; required "
                f"{args.minimum_mapped_base_fraction:.4%}"
            )
        canonical = merge_intervals(mapped_intervals, target_order)
        canonical_bed = work / "canonical.bed"
        with canonical_bed.open("w") as handle:
            for interval in canonical:
                handle.write("\t".join(map(str, interval)) + "\n")

        compressed = work / "canonical.bed.gz"
        with compressed.open("wb") as handle:
            subprocess.run([args.bgzip, "-c", str(canonical_bed)], stdout=handle, check=True)
        subprocess.run([args.tabix, "-f", "-p", "bed", str(compressed)], check=True)

        shutil.copyfile(mapped, args.mapping)
        shutil.copyfile(raw_unmapped, args.unmapped)
        os.replace(compressed, args.output)
        os.replace(Path(f"{compressed}.tbi"), Path(f"{args.output}.tbi"))

    missing_ids = sorted(source_ids - set(mapped_ids))
    mapped_id_counts = Counter(mapped_ids)
    split_ids = sorted(identifier for identifier, count in mapped_id_counts.items() if count > 1)
    report = {
        "schema_version": 1,
        "source_build": "hg19",
        "target_build": "grch38",
        "source_bed": str(args.source_bed),
        "source_bed_sha256": checksum(args.source_bed),
        "chain_sha256": checksum(args.chain),
        "output_sha256": checksum(args.output),
        "mapping": str(args.mapping),
        "source_interval_count": len(source),
        "mapped_interval_count": len(mapped_intervals),
        "canonical_interval_count": len(canonical),
        "source_bases": source_bases,
        "mapped_source_bases": mapped_source_bases,
        "canonical_bases": sum(end - start for _, start, end in canonical),
        "mapped_base_fraction": mapped_base_fraction,
        "unmapped_interval_ids": missing_ids,
        "split_interval_ids": split_ids,
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
