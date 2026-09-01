#!/usr/bin/env python3
"""Validate final intersect VCF count QC output."""

from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_vcf(path: Path, records: list[str]) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        *records,
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_gzipped_vcf(path: Path, records: list[str]) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        *records,
        "",
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        data_lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(data_lines, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        nonzero = tmp / "nonzero.intersect.vcf"
        zero = tmp / "zero.intersect.vcf"
        gzipped = tmp / "gzipped.intersect.vcf.gz"
        output = tmp / "final_variant_counts_mqc.tsv"

        write_vcf(
            nonzero,
            [
                "1\t10\t.\tA\tC\t.\tPASS\t.",
                "1\t20\t.\tG\tT\t.\tPASS\t.",
            ],
        )
        write_vcf(zero, [])
        write_gzipped_vcf(gzipped, ["2\t30\t.\tC\tA\t.\tPASS\t."])

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "final_variant_counts.py"),
                "--input",
                f"tumour_nonzero={nonzero}",
                "--input",
                f"tumour_zero={zero}",
                "--input",
                f"tumour_gzipped={gzipped}",
                "--output",
                str(output),
            ],
            cwd=ROOT,
            check=True,
        )

        rows = {row["tumour"]: row for row in read_rows(output)}
        expected = {
            "tumour_nonzero": {
                "tumour": "tumour_nonzero",
                "final_variant_count": "2",
                "zero_final_variants": "OK",
            },
            "tumour_zero": {
                "tumour": "tumour_zero",
                "final_variant_count": "0",
                "zero_final_variants": "WARN_ZERO_FINAL_VARIANTS",
            },
            "tumour_gzipped": {
                "tumour": "tumour_gzipped",
                "final_variant_count": "1",
                "zero_final_variants": "OK",
            },
        }
        if rows != expected:
            raise AssertionError(f"unexpected rows: {rows}")

        header = output.read_text(encoding="utf-8").splitlines()[:4]
        if "# id: final_variant_counts" not in header:
            raise AssertionError(f"missing MultiQC id header: {header}")
        if "# plot_type: table" not in header:
            raise AssertionError(f"missing MultiQC plot type header: {header}")


if __name__ == "__main__":
    main()
