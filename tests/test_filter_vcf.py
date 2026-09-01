#!/usr/bin/env python3
"""Validate final VCF depth/AF filtering."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


VCF = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR
1\t100\t.\tA\tC\t.\tPASS\t.\tGT:AD:AF\t0/0:20,0:0\t0/1:27,3:0.1
1\t101\t.\tA\tG\t.\tPASS\t.\tGT:AD:AF\t0/0:19,0:0\t0/1:27,3:0.1
1\t102\t.\tA\tT\t.\tPASS\t.\tGT:AD:AF\t0/0:20,0:0\t0/1:26,3:0.1
1\t103\t.\tA\tG\t.\tPASS\t.\tGT:AD:AF\t0/0:20,0:0\t0/1:28,2:0.09
"""


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        input_vcf = tmp / "input.vcf"
        output_vcf = tmp / "output.vcf"
        input_vcf.write_text(VCF)

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "filter_vcf.py"),
                "--input",
                str(input_vcf),
                "--output",
                str(output_vcf),
                "--tumour-sample",
                "TUMOUR",
                "--normal-sample",
                "NORMAL",
                "--af",
                "0.1",
                "--tumour-dp",
                "30",
                "--normal-dp",
                "20",
            ],
            cwd=ROOT,
            check=True,
        )

        records = [
            line.split("\t")[1]
            for line in output_vcf.read_text().splitlines()
            if not line.startswith("#")
        ]
        if records != ["100"]:
            raise AssertionError(f"unexpected passing records: {records}")


if __name__ == "__main__":
    main()
