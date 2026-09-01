#!/usr/bin/env python3
"""Exact-allele benchmark metrics for synthetic and SEQC2 truth comparisons."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from caller_tiers import read_vcf, write_vcf


def safe_ratio(numerator, denominator):
    return numerator / denominator if denominator else 0.0


def metrics(truth, calls, callable_bases):
    truth, calls = set(truth), set(calls)
    tp, fp, fn = len(truth & calls), len(calls - truth), len(truth - calls)
    precision = safe_ratio(tp, tp + fp)
    recall = safe_ratio(tp, tp + fn)
    return {
        "truth": len(truth),
        "calls": len(calls),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "f1": safe_ratio(2 * precision * recall, precision + recall),
        "false_positives_per_callable_gb": safe_ratio(fp, callable_bases / 1_000_000_000),
    }


def cosine_similarity(left, right):
    keys = set(left) | set(right)
    dot = sum(left.get(key, 0) * right.get(key, 0) for key in keys)
    left_norm = math.sqrt(sum(left.get(key, 0) ** 2 for key in keys))
    right_norm = math.sqrt(sum(right.get(key, 0) ** 2 for key in keys))
    return safe_ratio(dot, left_norm * right_norm)


def info_number(fields, key):
    for item in fields[7].split(";"):
        name, separator, value = item.partition("=")
        if separator and name == key:
            try:
                return float(value.split(",")[0])
            except ValueError:
                return None
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", required=True)
    parser.add_argument("--calls", required=True)
    parser.add_argument("--callable-bases", required=True, type=int)
    parser.add_argument("--output", required=True)
    parser.add_argument("--truth-vaf-field", default="VAF")
    parser.add_argument("--false-positive-vcf")
    args = parser.parse_args()
    _, truth = read_vcf(args.truth)
    calls_headers, calls = read_vcf(args.calls)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result = metrics(truth, calls, args.callable_bases)
    high_vaf = {
        key
        for key, fields in truth.items()
        if (info_number(fields, args.truth_vaf_field) or 0) >= 0.25
    }
    # Keep genuine lower-VAF truth matches out of the high-VAF stratum without
    # forgiving calls that are absent from the complete truth catalog.
    high_vaf_calls = (set(calls) & high_vaf) | (set(calls) - set(truth))
    result["truth_vaf_ge_0.25"] = metrics(high_vaf, high_vaf_calls, args.callable_bases)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.false_positive_vcf:
        false_positives = {key: fields for key, fields in calls.items() if key not in truth}
        write_vcf(
            args.false_positive_vcf,
            calls_headers,
            false_positives,
            {key: "FALSE_POSITIVE" for key in false_positives},
            sites_only=True,
        )


if __name__ == "__main__":
    main()
