#!/usr/bin/env python3
"""Annotate cohort candidates and emit reversible, reason-coded catalogs."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import yaml
from caller_tiers import open_text, read_vcf


def parse_info(value):
    result = {}
    if value not in {"", "."}:
        for item in value.split(";"):
            key, separator, data = item.partition("=")
            result[key] = data if separator else True
    return result


def population_values(path, wanted, field):
    values = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip().split("\t")
            chrom, pos, ref = columns[0], int(columns[1]), columns[3]
            alts = columns[4].split(",")
            info = parse_info(columns[7])
            frequencies = str(info.get(field, ".")).split(",")
            for index, alt in enumerate(alts):
                key = chrom, pos, ref, alt
                if key not in wanted:
                    continue
                raw = frequencies[index] if index < len(frequencies) else "."
                try:
                    values[key] = float(raw)
                except ValueError:
                    values[key] = 0.0
    return values


def read_counts(path):
    result = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chrom"], int(row["pos"]), row["ref"], row["alt"], row["sample"])
            result[key] = {
                name: float(row[name])
                for name in (
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
                )
            }
    return result


def read_masks(items):
    masks = {}
    for item in items:
        label, path = item.split("=", 1)
        intervals = {}
        with open_text(path) as handle:
            for line in handle:
                if not line.startswith("#") and line.strip():
                    chrom, start, end = line.split()[:3]
                    intervals.setdefault(chrom, []).append((int(start), int(end)))
        masks[label] = intervals
    return masks


def masked(key, intervals):
    chrom, pos, _, _ = key
    point = pos - 1
    return any(start <= point < end for start, end in intervals.get(chrom, []))


def add_annotations(fields, annotations):
    result = list(fields)
    additions = [f"{key}={value}" for key, value in annotations.items()]
    result[7] = ";".join(([result[7]] if result[7] not in {"", "."} else []) + additions)
    return result


def write_vcf(path, headers, rows):
    metadata = (
        '##INFO=<ID=FILTER_REASONS,Number=.,Type=String,Description="Reason-coded organoid filters">\n',
        '##INFO=<ID=POP_AF,Number=1,Type=Float,Description="Exact-allele population frequency">\n',
        '##INFO=<ID=RECURRENCE_TOTAL,Number=1,Type=Integer,Description="Later organoids carrying exact allele">\n',
        '##INFO=<ID=RECURRENCE_LINEAGE,Number=1,Type=Integer,Description="Same-lineage later organoids carrying exact allele">\n',
        '##INFO=<ID=RECURRENCE_UNRELATED,Number=1,Type=Integer,Description="Other-donor organoids carrying exact allele">\n',
        '##INFO=<ID=CARRIER_SAMPLES,Number=.,Type=String,Description="Later samples carrying exact allele">\n',
        '##INFO=<ID=LATER_EVIDENCE,Number=1,Type=String,Description="DP,REF,ALT,VAF for later sample">\n',
        '##INFO=<ID=BASELINE_EVIDENCE,Number=1,Type=String,Description="DP,REF,ALT,VAF for baseline">\n',
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        for line in headers:
            if line.startswith("#CHROM"):
                handle.writelines(metadata)
            handle.write(line)
        for _, fields, annotations, _ in rows:
            handle.write("\t".join(add_annotations(fields, annotations)) + "\n")
    os.replace(temporary, path)


def write_tsv(path, rows):
    columns = (
        "chrom",
        "pos",
        "ref",
        "alt",
        "caller_support",
        "later_depth",
        "later_ref",
        "later_alt",
        "later_vaf",
        "baseline_depth",
        "baseline_ref",
        "baseline_alt",
        "baseline_vaf",
        "population_af",
        "recurrence_total",
        "recurrence_lineage",
        "recurrence_unrelated",
        "carrier_samples",
        "filter_reasons",
    )
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row[3] for row in rows)
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--union", action="append", required=True, help="SAMPLE=VCF")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--population-vcf", required=True)
    parser.add_argument("--population-af-field", default="AF")
    parser.add_argument("--population-af-threshold", type=float, default=0.001)
    parser.add_argument("--min-later-alt", type=int, default=2)
    parser.add_argument("--min-later-vaf", type=float, default=0.20)
    parser.add_argument("--min-baseline-depth", type=int, default=6)
    parser.add_argument("--max-baseline-alt", type=int, default=1)
    parser.add_argument("--mask", action="append", default=[])
    for name in (
        "audit-vcf",
        "audit-tsv",
        "stringent-vcf",
        "stringent-tsv",
        "sensitivity-vcf",
        "sensitivity-tsv",
        "rejected-vcf",
        "rejected-tsv",
        "shared-lineage-vcf",
        "shared-lineage-tsv",
        "stage-counts",
    ):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(Path(args.manifest).read_text())
    baseline = manifest["comparisons"][args.sample]["baseline"]
    metadata = manifest["samples"]
    cohort = {}
    headers = None
    for item in args.union:
        sample, path = item.split("=", 1)
        current_headers, records = read_vcf(path)
        headers = headers or current_headers
        cohort[sample] = records
    records = cohort[args.sample]
    wanted = set().union(*(set(value) for value in cohort.values()))
    population = population_values(args.population_vcf, wanted, args.population_af_field)
    counts = read_counts(args.counts)
    masks = read_masks(args.mask)
    audit, stringent, sensitivity, rejected, shared = [], [], [], [], []
    stage = {"caller_union": len(records)}
    reason_sets = []
    for key, fields in records.items():
        carriers = sorted(sample for sample, variants in cohort.items() if key in variants)
        same_lineage = [
            sample
            for sample in carriers
            if metadata[sample]["lineage"] == metadata[args.sample]["lineage"]
        ]
        unrelated = [
            sample
            for sample in carriers
            if metadata[sample]["donor"] != metadata[args.sample]["donor"]
        ]
        later = counts.get((*key, args.sample), {})
        base = counts.get((*key, baseline), {})
        reasons = []
        if later.get("alt_count", 0) < args.min_later_alt:
            reasons.append("LOW_LATER_ALT")
        if later.get("vaf", 0) < args.min_later_vaf:
            reasons.append("LOW_LATER_VAF")
        if base.get("depth", 0) < args.min_baseline_depth:
            reasons.append("LOW_BASELINE_DEPTH")
        if base.get("alt_count", 0) > args.max_baseline_alt:
            reasons.append("BASELINE_ALT_EVIDENCE")
        if unrelated:
            reasons.append("RECURRENT_UNRELATED")
        if len(same_lineage) > 1:
            reasons.append("SHARED_LINEAGE")
        for label, intervals in masks.items():
            if masked(key, intervals):
                reasons.append(f"MASK_{label.upper()}")
        pop_af = population.get(key)
        common_reasons = list(reasons)
        if pop_af is not None:
            reasons.append("POPULATION_ANY")
        if pop_af is not None and pop_af > args.population_af_threshold:
            common_reasons.append("POPULATION_AF")
        info = parse_info(fields[7])
        annotations = {
            "FILTER_REASONS": ",".join(reasons) if reasons else "PASS",
            "POP_AF": "." if pop_af is None else f"{pop_af:.8g}",
            "RECURRENCE_TOTAL": len(carriers),
            "RECURRENCE_LINEAGE": len(same_lineage),
            "RECURRENCE_UNRELATED": len(unrelated),
            "CARRIER_SAMPLES": ",".join(carriers),
            "LATER_EVIDENCE": f"{int(later.get('depth', 0))},{int(later.get('ref_count', 0))},{int(later.get('alt_count', 0))},{later.get('vaf', 0):.6g}",
            "BASELINE_EVIDENCE": f"{int(base.get('depth', 0))},{int(base.get('ref_count', 0))},{int(base.get('alt_count', 0))},{base.get('vaf', 0):.6g}",
        }
        row = {
            "chrom": key[0],
            "pos": key[1],
            "ref": key[2],
            "alt": key[3],
            "caller_support": info.get("CALLER_SUPPORT", "."),
            "later_depth": int(later.get("depth", 0)),
            "later_ref": int(later.get("ref_count", 0)),
            "later_alt": int(later.get("alt_count", 0)),
            "later_vaf": later.get("vaf", 0),
            "baseline_depth": int(base.get("depth", 0)),
            "baseline_ref": int(base.get("ref_count", 0)),
            "baseline_alt": int(base.get("alt_count", 0)),
            "baseline_vaf": base.get("vaf", 0),
            "population_af": "." if pop_af is None else pop_af,
            "recurrence_total": len(carriers),
            "recurrence_lineage": len(same_lineage),
            "recurrence_unrelated": len(unrelated),
            "carrier_samples": ",".join(carriers),
            "filter_reasons": ",".join(reasons) if reasons else "PASS",
        }
        item = (key, fields, annotations, row)
        reason_sets.append(set(reasons))
        audit.append(item)
        (stringent if not reasons else rejected).append(item)
        if not common_reasons:
            sensitivity.append(item)
        if len(same_lineage) > 1:
            shared.append(item)
    for prefix, rows in (
        ("audit", audit),
        ("stringent", stringent),
        ("sensitivity", sensitivity),
        ("rejected", rejected),
        ("shared_lineage", shared),
    ):
        write_vcf(getattr(args, f"{prefix}_vcf"), headers, rows)
        write_tsv(getattr(args, f"{prefix}_tsv"), rows)
    active = list(reason_sets)
    stages = (
        ("later_alt_support", {"LOW_LATER_ALT"}),
        ("later_vaf", {"LOW_LATER_VAF"}),
        ("baseline_depth", {"LOW_BASELINE_DEPTH"}),
        ("baseline_alt_evidence", {"BASELINE_ALT_EVIDENCE"}),
        ("population_absent", {"POPULATION_ANY"}),
        ("recurrence", {"RECURRENT_UNRELATED", "SHARED_LINEAGE"}),
    )
    for name, blocked in stages:
        active = [reasons for reasons in active if not reasons & blocked]
        stage[name] = len(active)
    active = [
        reasons for reasons in active if not any(reason.startswith("MASK_") for reason in reasons)
    ]
    stage["region_masks"] = len(active)
    stage.update(
        {
            "stringent": len(stringent),
            "sensitivity": len(sensitivity),
            "rejected": len(rejected),
            "shared_lineage": len(shared),
        }
    )
    path = Path(args.stage_counts)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        "stage\tcount\n" + "".join(f"{name}\t{value}\n" for name, value in stage.items())
    )
    os.replace(temporary, path)


if __name__ == "__main__":
    main()
