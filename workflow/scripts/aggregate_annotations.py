#!/usr/bin/env python3
"""Build cohort tables directly from detailed annotated VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterator

import yaml

from annotation import parse_csq, parse_csq_fields, select_pick, validate_annotation_config


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def records(path: Path) -> Iterator[tuple[list[str], dict[str, str], list[str]]]:
    csq_fields: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_fields(line)
            elif line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"malformed VCF record in {path}")
                info = {}
                for item in fields[7].split(";"):
                    key, _, value = item.partition("=")
                    info[key] = value if value else "1"
                yield fields, info, csq_fields


def callable_bases(manifest: Path) -> int:
    value = json.loads(manifest.read_text()).get("territory_bases")
    if not isinstance(value, int) or value <= 0:
        raise ValueError("analysis manifest must contain positive territory_bases")
    return value


def write_burden(inputs: list[tuple[str, Path]], manifest: Path, output: Path) -> None:
    denominator = callable_bases(manifest)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sample", "callable_bases", "pass_snv_count", "pass_indel_count", "pass_total_count", "snv_per_mb", "indel_per_mb", "total_per_mb"))
        for sample, path in sorted(inputs):
            snvs = indels = 0
            for fields, _, _ in records(path):
                if fields[6] != "PASS":
                    continue
                for alt in fields[4].split(","):
                    if len(fields[3]) == len(alt) == 1:
                        snvs += 1
                    else:
                        indels += 1
            scale = 1_000_000 / denominator
            writer.writerow((sample, denominator, snvs, indels, snvs + indels,
                             f"{snvs * scale:.12g}", f"{indels * scale:.12g}", f"{(snvs + indels) * scale:.12g}"))


def write_recurrence(inputs: list[tuple[str, Path]], output: Path, minimum_carriers: int) -> None:
    if minimum_carriers < 1:
        raise ValueError("minimum carriers must be positive")
    variants: dict[tuple[str, int, str, str], dict[str, Any]] = defaultdict(lambda: {"samples": set(), "pick": {}})
    for sample, path in sorted(inputs):
        for fields, info, csq_fields in records(path):
            if fields[6] != "PASS":
                continue
            alts = fields[4].split(",")
            consequences = parse_csq(info.get("CSQ"), csq_fields) if csq_fields else []
            for alt in alts:
                key = (fields[0], int(fields[1]), fields[3], alt)
                variants[key]["samples"].add(sample)
                allele_csq = [item for item in consequences if item.get("Allele", alt) == alt]
                picked = select_pick(allele_csq or consequences)
                if picked and not variants[key]["pick"]:
                    variants[key]["pick"] = picked
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("variant_key", "chrom", "pos", "ref", "alt", "carrier_count", "cohort_frequency", "samples", "pick_gene", "pick_transcript", "pick_consequence"))
        cohort_size = len(inputs)
        for (chrom, pos, ref, alt), value in sorted(variants.items(), key=lambda item: (item[0][0], item[0][1], item[0][2], item[0][3])):
            carriers = sorted(value["samples"])
            if len(carriers) < minimum_carriers:
                continue
            pick = value["pick"]
            writer.writerow((f"{chrom}:{pos}:{ref}:{alt}", chrom, pos, ref, alt, len(carriers),
                             f"{len(carriers) / cohort_size:.12g}", ",".join(carriers),
                             pick.get("SYMBOL", pick.get("Gene", "")), pick.get("Feature", ""), pick.get("Consequence", "")))


def write_resources(config_path: Path, output: Path) -> None:
    config = yaml.safe_load(config_path.read_text()) or {}
    rows = validate_annotation_config(config, str(config.get("reference", {}).get("build", "")).lower())
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ("name", "build", "version", "path", "size", "checksum", "source", "licence", "access_date", "validation_status")
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_inputs(values: list[str]) -> list[tuple[str, Path]]:
    result = []
    for value in values:
        sample, separator, path = value.partition("=")
        if not separator or not sample or not path:
            raise ValueError("inputs must use SAMPLE=VCF")
        result.append((sample, Path(path)))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    burden = subparsers.add_parser("burden")
    burden.add_argument("--input", action="append", required=True)
    burden.add_argument("--manifest", required=True, type=Path)
    burden.add_argument("--output", required=True, type=Path)
    recurrence = subparsers.add_parser("recurrence")
    recurrence.add_argument("--input", action="append", required=True)
    recurrence.add_argument("--minimum-carriers", type=int, default=2)
    recurrence.add_argument("--output", required=True, type=Path)
    resources = subparsers.add_parser("resources")
    resources.add_argument("--config", required=True, type=Path)
    resources.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "burden":
        write_burden(parse_inputs(args.input), args.manifest, args.output)
    elif args.command == "recurrence":
        write_recurrence(parse_inputs(args.input), args.output, args.minimum_carriers)
    else:
        write_resources(args.config, args.output)


if __name__ == "__main__":
    main()
