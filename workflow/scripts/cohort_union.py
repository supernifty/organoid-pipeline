#!/usr/bin/env python3
"""Merge normalized per-organoid unions into one exact-allele cohort candidate VCF."""

import argparse

from caller_tiers import read_vcf, write_vcf


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", required=True, help="SAMPLE=VCF")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    headers = None
    records, carriers = {}, {}
    for item in args.input:
        sample, path = item.split("=", 1)
        current_headers, current = read_vcf(path)
        headers = headers or current_headers
        for key, fields in current.items():
            records.setdefault(key, fields)
            carriers.setdefault(key, []).append(sample)
    write_vcf(
        args.output,
        headers,
        records,
        {key: "COHORT:" + ",".join(sorted(value)) for key, value in carriers.items()},
        sites_only=True,
    )


if __name__ == "__main__":
    main()
