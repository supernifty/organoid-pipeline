#!/usr/bin/env python3
"""Count SBS96 channels from named normalized SNV VCF catalogs."""

from __future__ import annotations

import argparse
import os
import subprocess
from collections import Counter
from pathlib import Path

from caller_tiers import read_vcf

COMPLEMENT = str.maketrans("ACGT", "TGCA")
SUBSTITUTIONS = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


def canonical_channel(context, ref, alt):
    context, ref, alt = context.upper(), ref.upper(), alt.upper()
    if len(context) != 3 or any(base not in "ACGT" for base in context + ref + alt):
        return None
    if context[1] != ref:
        raise ValueError(f"Reference allele {ref} disagrees with FASTA context {context}")
    if ref in "AG":
        context = context.translate(COMPLEMENT)[::-1]
        ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
    return f"{context[0]}[{ref}>{alt}]{context[2]}"


def fetch_contexts(reference, keys, region_file):
    ordered = list(keys)
    region_file.write_text(
        "".join(f"{chrom}:{pos - 1}-{pos + 1}\n" for chrom, pos, _, _ in ordered)
    )
    result = subprocess.run(
        ["samtools", "faidx", "-r", str(region_file), reference],
        text=True,
        capture_output=True,
        check=True,
    )
    sequences, current = [], []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    if len(sequences) != len(ordered):
        raise ValueError("samtools faidx did not return every requested SBS context")
    return dict(zip(ordered, sequences, strict=True))


def channels():
    return [f"{left}[{sub}]{right}" for sub in SUBSTITUTIONS for left in "ACGT" for right in "ACGT"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--catalog", action="append", required=True, help="NAME=VCF")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    catalogs, all_keys = {}, []
    for item in args.catalog:
        name, path = item.split("=", 1)
        _, records = read_vcf(path)
        keys = [key for key in records if len(key[2]) == len(key[3]) == 1]
        catalogs[name] = keys
        all_keys.extend(keys)
    unique = list(dict.fromkeys(all_keys))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    regions = output.with_suffix(output.suffix + ".regions")
    try:
        contexts = fetch_contexts(args.reference, unique, regions)
    finally:
        regions.unlink(missing_ok=True)
    counts = {}
    for name, keys in catalogs.items():
        counts[name] = Counter(
            filter(None, (canonical_channel(contexts[key], key[2], key[3]) for key in keys))
        )
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write("channel\t" + "\t".join(catalogs) + "\n")
        for channel in channels():
            handle.write(
                channel + "\t" + "\t".join(str(counts[name][channel]) for name in catalogs) + "\n"
            )
    os.replace(temporary, output)


if __name__ == "__main__":
    main()
