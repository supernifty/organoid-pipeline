#!/usr/bin/env python3
"""Summarise final intersect VCF record counts for QC reporting."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


def count_variants(vcf: Path) -> int:
    count = 0
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.strip():
                count += 1
    return count


def parse_input(value: str) -> tuple[str, Path]:
    fields = value.split("=", maxsplit=1)
    if len(fields) != 2 or not fields[0] or not fields[1]:
        raise argparse.ArgumentTypeError(
            "--input values must be formatted as tumour=final_intersect_vcf"
        )
    return fields[0], Path(fields[1])


def write_multiqc_table(rows: list[tuple[str, int]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# id: final_variant_counts\n")
        handle.write("# section_name: Final variant counts\n")
        handle.write("# description: Final intersect VCF record counts. Zero final variants are reported as QC warnings.\n")
        handle.write("# plot_type: table\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["tumour", "final_variant_count", "zero_final_variants"])
        for tumour, count in rows:
            status = "WARN_ZERO_FINAL_VARIANTS" if count == 0 else "OK"
            writer.writerow([tumour, count, status])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=parse_input,
        help="Tumour and final VCF path formatted as tumour=path. May be repeated.",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    rows = [(tumour, count_variants(vcf)) for tumour, vcf in args.input]
    write_multiqc_table(rows, args.output)


if __name__ == "__main__":
    main()
