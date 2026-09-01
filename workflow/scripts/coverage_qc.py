#!/usr/bin/env python3
"""Validate exon resources and summarize mosdepth region/threshold output."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path


EXON_COLUMNS = ("CHROM", "START", "END", "EXON_ID", "GENE_ID", "GENE_SYMBOL")
THRESHOLDS = (10, 20, 50, 100)


def validate_exon_config(config: dict, reference_build: str) -> None:
    coverage = config.get("coverage", {})
    if not coverage.get("exon_enabled", False):
        return
    resource = coverage.get("exon_bed", {})
    required = ("build", "version", "path", "expected_sha256", "source", "licence", "access_date")
    missing = [key for key in required if not resource.get(key)]
    if missing:
        raise ValueError("coverage.exon_bed missing metadata: " + ", ".join(missing))
    if str(resource["build"]).lower() != reference_build:
        raise ValueError("coverage.exon_bed build does not match the reference")
    path = Path(resource["path"])
    if not path.is_file():
        raise ValueError(f"coverage.exon_bed does not exist: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != str(resource["expected_sha256"]).lower():
        raise ValueError("coverage.exon_bed checksum mismatch")
    read_exons(path)


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_exons(path: Path) -> list[dict[str, str]]:
    rows = []
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n").lstrip("#").split("\t")
        if tuple(header) != EXON_COLUMNS:
            raise ValueError("exon BED header must be: " + ", ".join(EXON_COLUMNS))
        seen = set()
        for number, line in enumerate(handle, 2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(EXON_COLUMNS):
                raise ValueError(f"malformed exon BED line {number}")
            row = dict(zip(EXON_COLUMNS, fields))
            start, end = int(row["START"]), int(row["END"])
            if start < 0 or end <= start or row["EXON_ID"] in seen:
                raise ValueError(f"invalid or duplicate exon on line {number}")
            seen.add(row["EXON_ID"])
            row["START"], row["END"] = start, end
            rows.append(row)
    if not rows:
        raise ValueError("exon BED is empty")
    return rows


def prepare_bed(source: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for row in read_exons(source):
            name = "|".join(str(row[key]) for key in ("EXON_ID", "GENE_ID", "GENE_SYMBOL"))
            handle.write(f"{row['CHROM']}\t{row['START']}\t{row['END']}\t{name}\n")


def load_mosdepth(regions: Path, thresholds: Path) -> list[dict[str, object]]:
    rows = {}
    with open_text(regions) as handle:
        for line in handle:
            chrom, start, end, name, mean = line.rstrip().split("\t")[:5]
            key = (chrom, int(start), int(end), name)
            rows[key] = {"chrom": chrom, "start": int(start), "end": int(end), "name": name, "mean_depth": float(mean)}
    with open_text(thresholds) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), int(fields[2]), fields[3])
            if key not in rows or len(fields) != 4 + len(THRESHOLDS):
                raise ValueError("mosdepth region/threshold outputs do not match")
            for depth, value in zip(THRESHOLDS, fields[4:]):
                rows[key][f"covered_{depth}x"] = int(value)
    if any(f"covered_{depth}x" not in row for row in rows.values() for depth in THRESHOLDS):
        raise ValueError("incomplete mosdepth threshold output")
    return list(rows.values())


def summarize(sample: str, regions: Path, thresholds: Path, exon_output: Path, gene_output: Path, warning_depth: int) -> None:
    if warning_depth not in THRESHOLDS:
        raise ValueError("warning depth must be one of 10, 20, 50, or 100")
    rows = load_mosdepth(regions, thresholds)
    exon_output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "exon_id", "gene_id", "gene_symbol", "chrom", "start", "end", "bases", "mean_depth"]
    for depth in THRESHOLDS:
        fields.extend((f"covered_bases_{depth}x", f"covered_fraction_{depth}x", f"complete_{depth}x"))
    fields.append("coverage_warning")
    genes = defaultdict(lambda: {"bases": 0, "depth_bases": 0.0, **{f"covered_{depth}x": 0 for depth in THRESHOLDS}})
    with exon_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["chrom"], item["start"], item["end"], item["name"])):
            exon_id, gene_id, gene_symbol = str(row["name"]).split("|", 2)
            bases = int(row["end"]) - int(row["start"])
            output = {"sample": sample, "exon_id": exon_id, "gene_id": gene_id, "gene_symbol": gene_symbol,
                      "chrom": row["chrom"], "start": row["start"], "end": row["end"], "bases": bases,
                      "mean_depth": f"{float(row['mean_depth']):.6f}"}
            gene = genes[(gene_id, gene_symbol)]
            gene["bases"] += bases
            gene["depth_bases"] += float(row["mean_depth"]) * bases
            for depth in THRESHOLDS:
                covered = int(row[f"covered_{depth}x"])
                output[f"covered_bases_{depth}x"] = covered
                output[f"covered_fraction_{depth}x"] = f"{covered / bases:.6f}"
                output[f"complete_{depth}x"] = str(covered == bases).lower()
                gene[f"covered_{depth}x"] += covered
            output["coverage_warning"] = str(int(row[f"covered_{warning_depth}x"]) != bases).lower()
            writer.writerow(output)
    with gene_output.open("w", newline="") as handle:
        gene_fields = ["sample", "gene_id", "gene_symbol", "bases", "mean_depth"]
        for depth in THRESHOLDS:
            gene_fields.extend((f"covered_bases_{depth}x", f"covered_fraction_{depth}x", f"complete_{depth}x"))
        gene_fields.append("coverage_warning")
        writer = csv.DictWriter(handle, fieldnames=gene_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (gene_id, symbol), value in sorted(genes.items()):
            bases = value["bases"]
            output = {"sample": sample, "gene_id": gene_id, "gene_symbol": symbol, "bases": bases,
                      "mean_depth": f"{value['depth_bases'] / bases:.6f}"}
            for depth in THRESHOLDS:
                covered = value[f"covered_{depth}x"]
                output[f"covered_bases_{depth}x"] = covered
                output[f"covered_fraction_{depth}x"] = f"{covered / bases:.6f}"
                output[f"complete_{depth}x"] = str(covered == bases).lower()
            output["coverage_warning"] = str(value[f"covered_{warning_depth}x"] != bases).lower()
            writer.writerow(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare-bed")
    prepare.add_argument("--source", required=True, type=Path)
    prepare.add_argument("--output", required=True, type=Path)
    summary = sub.add_parser("summarize")
    summary.add_argument("--sample", required=True)
    summary.add_argument("--regions", required=True, type=Path)
    summary.add_argument("--thresholds", required=True, type=Path)
    summary.add_argument("--exon-output", required=True, type=Path)
    summary.add_argument("--gene-output", required=True, type=Path)
    summary.add_argument("--warning-depth", required=True, type=int)
    args = parser.parse_args()
    if args.command == "prepare-bed":
        prepare_bed(args.source, args.output)
    else:
        summarize(args.sample, args.regions, args.thresholds, args.exon_output, args.gene_output, args.warning_depth)


if __name__ == "__main__":
    main()
