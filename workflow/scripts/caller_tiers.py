#!/usr/bin/env python3
"""Build deterministic exact-allele caller tiers without dropping provenance."""

from __future__ import annotations

import argparse
import gzip
import os
import re
from pathlib import Path


def open_text(path, mode="rt"):
    return (
        gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")
    )


def read_vcf(path):
    headers, records = [], {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record in {path}: {line.rstrip()}")
            alts = fields[4].split(",")
            if len(alts) != 1:
                raise ValueError(f"VCF must be split before caller comparison: {path}")
            key = (fields[0], int(fields[1]), fields[3], alts[0])
            if key in records:
                raise ValueError(f"Duplicate normalized allele in {path}: {key}")
            records[key] = fields
    if not any(line.startswith("#CHROM") for line in headers):
        raise ValueError(f"Missing VCF column header: {path}")
    return headers, records


def passed(fields):
    return fields[6] in {"PASS", "."}


def add_support(fields, support):
    result = list(fields)
    tag = f"CALLER_SUPPORT={support}"
    result[7] = tag if result[7] in {"", "."} else f"{result[7]};{tag}"
    return result


def write_vcf(path, headers, records, support, sites_only=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        inserted = False
        for line in headers:
            if line.startswith("#CHROM") and not inserted:
                handle.write(
                    '##INFO=<ID=CALLER_SUPPORT,Number=.,Type=String,Description="PASS-supporting callers">\n'
                )
                inserted = True
            if sites_only and line.startswith("#CHROM"):
                handle.write("\t".join(line.rstrip("\n").split("\t")[:8]) + "\n")
            else:
                handle.write(line)
        contigs = [
            re.search(r"ID=([^,>]+)", line).group(1)
            for line in headers
            if line.startswith("##contig=<")
        ]
        rank = {contig: index for index, contig in enumerate(contigs)}
        for key in sorted(
            records,
            key=lambda item: (rank.get(item[0], len(rank)), item[0], item[1], item[2], item[3]),
        ):
            fields = add_support(records[key], support[key])
            handle.write("\t".join(fields[:8] if sites_only else fields) + "\n")
    os.replace(temporary, path)


def build_tiers(mutect2, strelka_paths):
    headers, all_mutect2 = read_vcf(mutect2)
    mutect2_pass = {key: fields for key, fields in all_mutect2.items() if passed(fields)}
    strelka_pass = {}
    for path in strelka_paths:
        _, records = read_vcf(path)
        for key, fields in records.items():
            if passed(fields):
                strelka_pass.setdefault(key, fields)
    both = set(mutect2_pass) & set(strelka_pass)
    union = set(mutect2_pass) | set(strelka_pass)
    combined = {
        key: mutect2_pass[key] if key in mutect2_pass else strelka_pass[key] for key in union
    }
    support = {
        key: "Mutect2,Strelka2" if key in both else "Mutect2" if key in mutect2_pass else "Strelka2"
        for key in union
    }
    return headers, mutect2_pass, strelka_pass, both, combined, support


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutect2", required=True)
    parser.add_argument("--strelka", nargs="+", required=True)
    for name in (
        "mutect2-pass",
        "strelka-pass",
        "intersection",
        "union",
        "mutect2-only",
        "strelka-only",
    ):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    headers, m2, strelka, both, union, support = build_tiers(args.mutect2, args.strelka)
    strelka_headers, _ = read_vcf(args.strelka[0])
    outputs = {
        args.mutect2_pass: (headers, m2, {key: "Mutect2" for key in m2}, False),
        args.strelka_pass: (strelka_headers, strelka, {key: "Strelka2" for key in strelka}, False),
        args.intersection: (
            headers,
            {key: union[key] for key in both},
            {key: "Mutect2,Strelka2" for key in both},
            False,
        ),
        args.union: (headers, union, support, True),
        args.mutect2_only: (
            headers,
            {key: m2[key] for key in set(m2) - set(strelka)},
            {key: "Mutect2" for key in set(m2) - set(strelka)},
            False,
        ),
        args.strelka_only: (
            strelka_headers,
            {key: strelka[key] for key in set(strelka) - set(m2)},
            {key: "Strelka2" for key in set(strelka) - set(m2)},
            False,
        ),
    }
    for path, (output_headers, records, labels, sites_only) in outputs.items():
        write_vcf(path, output_headers, records, labels, sites_only)


if __name__ == "__main__":
    main()
