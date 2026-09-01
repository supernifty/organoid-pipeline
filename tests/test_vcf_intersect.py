#!/usr/bin/env python3
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "vcf_intersect.py"
HEADER = """##fileformat=VCFv4.2
##contig=<ID=1,length=1000>
##contig=<ID=2,length=1000>
##contig=<ID=10,length=1000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def write_vcf(path, records):
    text = HEADER + "".join("\t".join(record) + "\n" for record in records)
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def variants(path):
    return [line.rstrip("\n").split("\t") for line in path.read_text().splitlines() if not line.startswith("#")]


def main():
    with tempfile.TemporaryDirectory() as directory:
        directory = Path(directory)
        mutect2 = directory / "mutect2.vcf.gz"
        snvs = directory / "snvs.vcf.gz"
        indels = directory / "indels.vcf.gz"
        output = directory / "output.vcf"

        write_vcf(mutect2, [
            ("1", "10", ".", "A", "C", ".", "PASS", "."),
            ("2", "20", ".", "G", "A", ".", "artifact", "."),
            ("2", "20", ".", "G", "T", ".", "LowDepth", "."),
            ("10", "30", ".", "AT", "A", ".", "PASS", "."),
        ])
        write_vcf(snvs, [
            ("1", "10", ".", "A", "C", ".", "PASS", "."),
            ("2", "20", ".", "G", "A", ".", "PASS", "."),
            ("2", "20", ".", "G", "T", ".", "LowDepth", "."),
        ])
        write_vcf(indels, [
            ("10", "30", ".", "AT", "A", ".", "PASS", "."),
        ])

        completed = subprocess.run(
            [
                sys.executable, str(SCRIPT),
                "--mutect2-vcf", str(mutect2),
                "--strelka-vcf", str(snvs), str(indels),
                "--allowed-filters", "str_contraction", "LowDepth",
                "--output-vcf", str(output),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "Intersected 3 variants from 4 Mutect2 records" in completed.stdout
        observed = [(record[0], record[1], record[3], record[4]) for record in variants(output)]
        assert observed == [
            ("1", "10", "A", "C"),
            ("2", "20", "G", "T"),
            ("10", "30", "AT", "A"),
        ]

        empty = directory / "empty.vcf.gz"
        empty_output = directory / "empty-output.vcf"
        write_vcf(empty, [])
        subprocess.run(
            [
                sys.executable, str(SCRIPT),
                "--mutect2-vcf", str(mutect2),
                "--strelka-vcf", str(empty),
                "--output-vcf", str(empty_output),
            ],
            check=True,
        )
        assert variants(empty_output) == []


if __name__ == "__main__":
    main()
