#!/usr/bin/env python3
"""Deterministically downsample a name-paired alignment and measure achieved depth."""

from __future__ import annotations

import argparse
import gzip
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
