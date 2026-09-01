#!/usr/bin/env python3
"""Validate extended-context aggregate counting."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_fasta_and_index(tmp: Path) -> tuple[Path, Path]:
    bases = list("C" * 60)
    bases[0] = "A"
    bases[1] = "A"
    bases[4] = "T"
    bases[14] = "A"
    bases[17] = "T"
    bases[18] = "T"
    bases[24] = "A"
    bases[34] = "C"
    bases[39] = "G"
    sequence = "".join(bases)

    fasta = tmp / "genome.fa"
    fasta.write_text(f">chr1\n{sequence}\n")

    fai = tmp / "genome.fa.fai"
    offset = len(">chr1\n")
    fai.write_text(f"chr1\t{len(sequence)}\t{offset}\t{len(sequence)}\t{len(sequence) + 1}\n")
    return fasta, fai


def write_vcf(tmp: Path) -> Path:
    vcf = tmp / "sample.intersect.vcf"
    vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr1\t5\t.\tT\tC\t.\tPASS\t.",
                "chr1\t15\t.\tA\tG\t.\tPASS\t.",
                "chr1\t25\t.\tA\tG\t.\tPASS\t.",
                "chr1\t35\t.\tC\tT\t.\tPASS\t.",
                "chr1\t40\t.\tG\tGT\t.\tPASS\t.",
                "",
            ]
        )
    )
    return vcf


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        fasta, fai = write_fasta_and_index(tmp)
        vcf = write_vcf(tmp)
        output = tmp / "extended_contexts.tsv"

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "extended_contexts.py"),
                "--reference",
                str(fasta),
                "--reference-index",
                str(fai),
                "--input",
                f"tumour_a={vcf}",
                "--output",
                str(output),
            ],
            cwd=ROOT,
            check=True,
        )

        with output.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        if len(rows) != 1:
            raise AssertionError(f"expected one row, got {len(rows)}")
        row = rows[0]
        expected = {
            "tumour": "tumour_a",
            "variant_count": "5",
            "snv_count": "4",
            "indel_count": "1",
            "colibactin_aat_count": "2",
            "colibactin_at_snv_count": "3",
            "colibactin_snv_proportion": "0.666667",
        }
        if row != expected:
            raise AssertionError(f"unexpected row: {row}")


if __name__ == "__main__":
    main()
