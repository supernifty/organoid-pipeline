#!/usr/bin/env python3
"""Deterministically downsample a name-paired alignment and measure achieved depth."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
from pathlib import Path


def weighted_depth(path):
    total, bases = 0.0, 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            length = int(fields[2]) - int(fields[1])
            total += length * float(fields[-1])
            bases += length
    if not bases:
        raise ValueError("mosdepth reported no callable bases")
    return total / bases


def fai_lengths(path):
    return {
        fields[0]: int(fields[1])
        for fields in (line.split("\t") for line in Path(path).read_text().splitlines())
    }


def parse_alignment_header(text):
    sort_order = None
    sequences = {}
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
            sort_order = values.get("SO")
        elif fields[0] == "@SQ":
            values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
            if "SN" in values and "LN" in values:
                sequences[values["SN"]] = int(values["LN"])
    return sort_order, sequences


def territory_contigs(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("@", "#", "track", "browser")):
                continue
            contig = line.split("\t", 1)[0]
            if contig not in contigs:
                contigs.append(contig)
    if not contigs:
        raise ValueError(f"analysis territory has no contigs: {path}")
    return contigs


def validate_dictionaries(alignment_lengths, reference_lengths, required_contigs, mapped):
    for contig in required_contigs:
        if alignment_lengths.get(contig) != reference_lengths.get(contig):
            raise ValueError(
                f"alignment and configured reference disagree for required contig {contig}: "
                f"{alignment_lengths.get(contig)!r} != {reference_lengths.get(contig)!r}"
            )
    incompatible = [
        contig
        for contig, count in mapped.items()
        if count > 0 and alignment_lengths.get(contig) != reference_lengths.get(contig)
    ]
    if incompatible:
        preview = ", ".join(incompatible[:10])
        suffix = " ..." if len(incompatible) > 10 else ""
        raise ValueError(
            "mapped reads use contigs absent from or incompatible with the configured "
            f"reference: {preview}{suffix}"
        )


def idxstats_command(alignment, reference):
    command = ["samtools", "idxstats"]
    if Path(alignment).suffix.lower() == ".cram":
        command.extend(("--input-fmt-option", f"reference={reference}"))
    command.append(str(alignment))
    return command


def checked_output(command, operation):
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or "no stderr was produced"
        raise RuntimeError(f"{operation} failed (exit {completed.returncode}): {detail}")
    return completed.stdout


def validate_input_reference(alignment, reference, territory):
    fai = Path(f"{reference}.fai")
    if not fai.is_file():
        raise ValueError(f"reference FASTA index is missing: {fai}")
    subprocess.run(["samtools", "quickcheck", "-v", alignment], check=True)
    header = checked_output(
        ["samtools", "view", "-H", alignment], "samtools alignment-header validation"
    )
    sort_order, alignment_lengths = parse_alignment_header(header)
    if sort_order != "coordinate":
        raise ValueError(f"input alignment is not coordinate sorted (SO={sort_order!r})")
    idxstats = checked_output(
        idxstats_command(alignment, reference), "samtools alignment-index validation"
    )
    mapped = {}
    for line in idxstats.splitlines():
        fields = line.split("\t")
        if len(fields) >= 4 and fields[0] != "*":
            mapped[fields[0]] = int(fields[2])
    reference_lengths = fai_lengths(fai)
    required = territory_contigs(territory)
    validate_dictionaries(alignment_lengths, reference_lengths, required, mapped)
    return {
        "sort_order": sort_order,
        "required_contigs": required,
        "mapped_contig_count": sum(count > 0 for count in mapped.values()),
        "reference_fai_sha256": hashlib.sha256(fai.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--territory", required=True)
    parser.add_argument("--input-depth", required=True, type=float)
    parser.add_argument("--target-depth", required=True, type=float)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    compatibility = validate_input_reference(args.input, args.reference, args.territory)
    fraction = args.target_depth / args.input_depth
    if not 0 < fraction <= 1:
        raise ValueError("target depth must be positive and no greater than measured input depth")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp.cram")
    # samtools -s hashes the template name, so mates receive the same deterministic decision.
    command = ["samtools", "view", "-@", str(args.threads), "-T", args.reference, "-C"]
    if fraction < 1:
        subsample = f"{args.seed}.{str(fraction).split('.', 1)[1][:8]}"
        command.extend(("-s", subsample))
    command.extend(("-o", str(temporary), args.input))
    subprocess.run(command, check=True)
    subprocess.run(["samtools", "index", "-@", str(args.threads), str(temporary)], check=True)
    prefix = output.with_suffix(output.suffix + ".mosdepth")
    subprocess.run(
        [
            "mosdepth",
            "--threads",
            str(args.threads),
            "--no-per-base",
            "--by",
            args.territory,
            "--fasta",
            args.reference,
            str(prefix),
            str(temporary),
        ],
        check=True,
    )
    achieved = weighted_depth(Path(f"{prefix}.regions.bed.gz"))
    os.replace(temporary, output)
    os.replace(Path(f"{temporary}.crai"), Path(f"{output}.crai"))
    Path(args.report).write_text(
        json.dumps(
            {
                "input": args.input,
                "seed": args.seed,
                "input_depth": args.input_depth,
                "target_depth": args.target_depth,
                "fraction": fraction,
                "achieved_depth": achieved,
                "reference_compatibility": compatibility,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    for suffix in (
        ".mosdepth.global.dist.txt",
        ".mosdepth.region.dist.txt",
        ".mosdepth.summary.txt",
        ".regions.bed.gz",
        ".regions.bed.gz.csi",
    ):
        Path(f"{prefix}{suffix}").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
