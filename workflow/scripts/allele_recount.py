#!/usr/bin/env python3
"""Recount normalized exact alleles with one samtools mpileup pass per sample.

Overlapping mates are collapsed by samtools's default overlap handling. Reads and
bases below the configured mapping/base quality are excluded before counting.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
from collections import defaultdict
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def candidates(path):
    result = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            result.append((fields[0], int(fields[1]), fields[3].upper(), fields[4].upper()))
    return result


def decode_bases(text, reference):
    """Return (observed base, strand, attached indel) observations."""
    observations = []
    index = 0
    while index < len(text):
        char = text[index]
        if char == "^":
            index += 2
            continue
        if char == "$":
            index += 1
            continue
        if char in "<>*#":
            index += 1
            continue
        if char in ".,":
            base, strand = reference.upper(), "+" if char == "." else "-"
        elif char.upper() in "ACGTN":
            base, strand = char.upper(), "+" if char.isupper() else "-"
        else:
            index += 1
            continue
        index += 1
        indel = ""
        if index < len(text) and text[index] in "+-":
            sign = text[index]
            index += 1
            start = index
            while index < len(text) and text[index].isdigit():
                index += 1
            if start == index:
                raise ValueError(f"Malformed mpileup indel length in {text!r}")
            length = int(text[start:index])
            indel = sign + text[index : index + length].upper()
            index += length
        observations.append((base, strand, indel))
    return observations


def event_for(ref, alt):
    if len(ref) == len(alt) == 1:
        return "snv", alt
    if alt.startswith(ref):
        return "insertion", "+" + alt[len(ref) :]
    if ref.startswith(alt):
        return "deletion", "-" + ref[len(alt) :]
    return "complex", alt


def count_observations(ref, alt, bases):
    kind, event = event_for(ref, alt)
    observations = decode_bases(bases, ref[0])
    ref_fwd = ref_rev = alt_fwd = alt_rev = 0
    for base, strand, indel in observations:
        is_alt = (kind == "snv" and base == event and not indel) or (
            kind in {"insertion", "deletion"} and base == ref[0] and indel == event
        )
        is_ref = base == ref[0] and not indel
        if is_alt:
            alt_fwd += strand == "+"
            alt_rev += strand == "-"
        elif is_ref:
            ref_fwd += strand == "+"
            ref_rev += strand == "-"
    return ref_fwd, ref_rev, alt_fwd, alt_rev, len(observations)


def parse_pileup(lines, wanted):
    result = {}
    for line in lines:
        fields = line.rstrip().split("\t")
        if len(fields) < 6:
            continue
        locus = (fields[0], int(fields[1]))
        for ref, alt in wanted.get(locus, []):
            counts = count_observations(ref, alt, fields[4])
            qualities = [ord(char) - 33 for char in fields[5]]
            mappings = [ord(char) - 33 for char in fields[6]] if len(fields) > 6 else []
            result[(fields[0], int(fields[1]), ref, alt)] = (
                *counts,
                sum(qualities) / len(qualities) if qualities else 0.0,
                sum(mappings) / len(mappings) if mappings else 0.0,
            )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", action="append", required=True, help="SAMPLE=CRAM")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    args = parser.parse_args()
    alleles = candidates(args.candidates)
    wanted = defaultdict(list)
    for chrom, pos, ref, alt in alleles:
        wanted[(chrom, pos)].append((ref, alt))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    bed = output.with_suffix(output.suffix + ".positions.bed")
    with bed.open("w", encoding="utf-8") as handle:
        for chrom, pos in wanted:
            handle.write(f"{chrom}\t{pos - 1}\t{pos}\n")
    rows = []
    try:
        for item in args.sample:
            sample, alignment = item.split("=", 1)
            command = [
                "samtools",
                "mpileup",
                "-aa",
                "-B",
                "-s",
                "-q",
                str(args.min_mapq),
                "-Q",
                str(args.min_baseq),
                "-l",
                str(bed),
                "-f",
                args.reference,
                alignment,
            ]
            process = subprocess.run(command, text=True, capture_output=True, check=True)
            counts = parse_pileup(process.stdout.splitlines(), wanted)
            for allele in alleles:
                ref_fwd, ref_rev, alt_fwd, alt_rev, observed, baseq, mapq = counts.get(
                    allele, (0, 0, 0, 0, 0, 0.0, 0.0)
                )
                ref_count, alt_count = ref_fwd + ref_rev, alt_fwd + alt_rev
                depth = ref_count + alt_count
                rows.append(
                    (
                        *allele,
                        sample,
                        depth,
                        ref_count,
                        alt_count,
                        alt_count / depth if depth else 0.0,
                        ref_fwd,
                        ref_rev,
                        alt_fwd,
                        alt_rev,
                        baseq,
                        mapq,
                        observed,
                    )
                )
    finally:
        bed.unlink(missing_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "chrom",
                "pos",
                "ref",
                "alt",
                "sample",
                "depth",
                "ref_count",
                "alt_count",
                "vaf",
                "ref_forward",
                "ref_reverse",
                "alt_forward",
                "alt_reverse",
                "mean_baseq",
                "mean_mapq",
                "observed_bases",
            )
        )
        writer.writerows(rows)
    os.replace(temporary, output)


if __name__ == "__main__":
    main()
