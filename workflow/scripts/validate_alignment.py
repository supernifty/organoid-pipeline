#!/usr/bin/env python3
"""Fail-closed validation of CRAM/index/header/reference compatibility."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


def fai_lengths(path):
    return {
        fields[0]: int(fields[1])
        for fields in (line.split("\t") for line in Path(path).read_text().splitlines())
    }


def parse_header(text):
    sort_order = None
    sequence, samples = {}, set()
    for line in text.splitlines():
        fields = line.split("\t")
        if fields[0] == "@HD":
            values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
            sort_order = values.get("SO")
        elif fields[0] == "@SQ":
            values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
            if "SN" in values and "LN" in values:
                sequence[values["SN"]] = int(values["LN"])
        elif fields[0] == "@RG":
            values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
            if values.get("SM"):
                samples.add(values["SM"])
    return sort_order, sequence, samples


def checked_output(command, operation):
    completed = subprocess.run(command, text=True, capture_output=True)
    if completed.returncode:
        detail = completed.stderr.strip() or "no stderr was produced"
        raise ValueError(f"{operation} failed (exit {completed.returncode}): {detail}")
    return completed.stdout


def idxstats_command(cram, reference):
    return [
        "samtools",
        "idxstats",
        "--input-fmt-option",
        f"reference={reference}",
        cram,
    ]


def cram_version(path):
    with Path(path).open("rb") as handle:
        header = handle.read(6)
    if len(header) != 6 or header[:4] != b"CRAM":
        raise ValueError(f"Alignment does not have a valid CRAM file header: {path}")
    return f"{header[4]}.{header[5]}"


def validate(cram, crai, reference, fai, expected_sample, contigs, expected_cram_version):
    for path in (cram, crai, reference, fai):
        if not Path(path).is_file() or not os.access(path, os.R_OK):
            raise ValueError(f"Required alignment resource is not readable: {path}")
    observed_cram_version = cram_version(cram)
    if observed_cram_version != expected_cram_version:
        raise ValueError(
            f"{cram} uses CRAM {observed_cram_version}, but the configured tools require "
            f"CRAM {expected_cram_version}; regenerate the alignment in the configured format"
        )
    checked_output(["samtools", "quickcheck", "-v", cram], "samtools CRAM quickcheck")
    header = checked_output(
        ["samtools", "view", "-H", "-T", reference, cram],
        "samtools CRAM header validation",
    )
    sort_order, sequence, samples = parse_header(header)
    if sort_order != "coordinate":
        raise ValueError(f"{cram} is not coordinate sorted (SO={sort_order!r})")
    if samples != {expected_sample}:
        raise ValueError(
            f"{cram} read-group SM values {sorted(samples)} do not equal {expected_sample!r}"
        )
    reference_lengths = fai_lengths(fai)
    for contig in contigs:
        if sequence.get(contig) != reference_lengths.get(contig):
            raise ValueError(f"{cram} and reference disagree for contig {contig}")
    checked_output(idxstats_command(cram, reference), "samtools CRAM index validation")
    return {
        "cram": cram,
        "crai": crai,
        "sample": expected_sample,
        "sort_order": sort_order,
        "cram_version": observed_cram_version,
        "validated_contigs": contigs,
        "reference_fai_sha256": hashlib.sha256(Path(fai).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cram", required=True)
    parser.add_argument("--crai", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--contig", action="append", required=True)
    parser.add_argument("--cram-version", default="3.0")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = validate(
        args.cram,
        args.crai,
        args.reference,
        args.fai,
        args.sample,
        args.contig,
        args.cram_version,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, output)


if __name__ == "__main__":
    main()
