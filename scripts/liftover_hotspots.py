#!/usr/bin/env python3
"""Lift a hotspot TSV between assemblies with Picard LiftoverVcf validation."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def read_fai(path: Path) -> list[tuple[str, int]]:
    with path.open() as handle:
        return [(fields[0], int(fields[1])) for fields in map(lambda line: line.rstrip().split("\t"), handle)]


def chain_contig(contig: str) -> str:
    if contig.startswith("chr"):
        return contig
    if contig in {"M", "MT"}:
        return "chrM"
    return f"chr{contig}"


def reference_base(samtools: str, reference: Path, chrom: str, pos: str) -> str:
    result = subprocess.run(
        [samtools, "faidx", str(reference), f"{chrom}:{pos}-{pos}"],
        check=True,
        text=True,
        capture_output=True,
    )
    sequence = "".join(line.strip() for line in result.stdout.splitlines() if not line.startswith(">"))
    if len(sequence) != 1:
        raise ValueError(f"Could not retrieve reference base for {chrom}:{pos}")
    return sequence.upper()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--source-reference", required=True, type=Path)
    parser.add_argument("--source-fai", required=True, type=Path)
    parser.add_argument("--target-reference", required=True, type=Path)
    parser.add_argument("--chain", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rejected-vcf", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--picard", default="picard")
    parser.add_argument("--samtools", default="samtools")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for executable in (args.picard, args.samtools):
        if not shutil.which(executable):
            raise ValueError(f"Required executable is unavailable: {executable}")

    with args.source.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "pos", "ref", "alt", "name"}
        if reader.fieldnames is None or set(reader.fieldnames) < required:
            raise ValueError(f"Hotspot TSV must contain: {', '.join(sorted(required))}")
        rows = list(reader)
    if not rows:
        raise ValueError("Hotspot TSV is empty")

    source_dictionary = dict(read_fai(args.source_fai))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.rejected_vcf.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hotspot-liftover.", dir=args.output.parent) as directory:
        work = Path(directory)
        source_vcf = work / "source.vcf"
        lifted_vcf = work / "lifted.vcf"
        rejected_vcf = work / "rejected.vcf"
        names: dict[str, str] = {}
        with source_vcf.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            for contig, length in source_dictionary.items():
                handle.write(f"##contig=<ID={chain_contig(contig)},length={length}>\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for index, row in enumerate(rows):
                identifier = f"HOTSPOT{index:08d}"
                chrom, pos = row["chrom"].strip(), row["pos"].strip()
                ref, alt = row["ref"].strip().upper(), row["alt"].strip().upper()
                if chrom not in source_dictionary:
                    raise ValueError(f"Hotspot contig is absent from source reference: {chrom}")
                if reference_base(args.samtools, args.source_reference, chrom, pos) != ref:
                    raise ValueError(f"Hotspot REF does not match source reference: {chrom}:{pos} {ref}")
                names[identifier] = row["name"].strip()
                handle.write(
                    f"{chain_contig(chrom)}\t{pos}\t{identifier}\t{ref}\t{alt}\t.\tPASS\t.\n"
                )

        subprocess.run(
            [
                args.picard, "LiftoverVcf",
                f"I={source_vcf}", f"O={lifted_vcf}", f"CHAIN={args.chain}",
                f"REJECT={rejected_vcf}", f"R={args.target_reference}",
                "CREATE_INDEX=false",
            ],
            check=True,
        )

        lifted = []
        with lifted_vcf.open() as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip().split("\t")
                chrom, pos, identifier, ref, alt = fields[:5]
                if identifier not in names:
                    raise ValueError(f"Unexpected lifted hotspot identifier: {identifier}")
                if reference_base(args.samtools, args.target_reference, chrom, pos) != ref.upper():
                    raise ValueError(f"Lifted REF does not match target reference: {chrom}:{pos} {ref}")
                lifted.append((chrom, pos, ref.upper(), alt.upper(), names[identifier]))

        with args.output.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "pos", "ref", "alt", "name"])
            writer.writerows(lifted)
        shutil.copyfile(rejected_vcf, args.rejected_vcf)

    report = {
        "schema_version": 1,
        "source_hotspot_count": len(rows),
        "lifted_hotspot_count": len(lifted),
        "rejected_hotspot_count": len(rows) - len(lifted),
        "source": str(args.source),
        "output": str(args.output),
        "rejected_vcf": str(args.rejected_vcf),
    }
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
