#!/usr/bin/env python3
"""Aggregate exact hotspot calls across final and caller-specific VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


PASSABLE_FILTERS = {"", ".", "PASS"}


@dataclass(frozen=True)
class Hotspot:
    chrom: str
    pos: str
    ref: str
    alt: str
    name: str

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (normalize_chrom(self.chrom), self.pos, self.ref.upper(), self.alt.upper())


@dataclass
class Hit:
    caller: str
    filter: str
    af: str = ""
    depth: str = ""

    @property
    def passable(self) -> bool:
        return self.filter in PASSABLE_FILTERS


@dataclass(frozen=True)
class SampleInput:
    tumour: str
    intersect: Path
    mutect2: Path | None = None
    strelka: Path | None = None


def normalize_chrom(chrom: str) -> str:
    normalized = chrom.split(".", maxsplit=1)[0].split("_", maxsplit=1)[-1]
    if normalized.lower().startswith("chr"):
        normalized = normalized[3:]
    normalized = normalized.upper()
    if normalized == "MT":
        return "M"
    return normalized


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open(encoding="utf-8")


def read_hotspots(path: Path) -> list[Hotspot]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "pos", "ref", "alt", "name"}
        if reader.fieldnames is None or set(reader.fieldnames) < required:
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        return [
            Hotspot(
                chrom=row["chrom"].strip(),
                pos=row["pos"].strip(),
                ref=row["ref"].strip(),
                alt=row["alt"].strip(),
                name=row["name"].strip(),
            )
            for row in reader
            if row.get("chrom") and row.get("pos") and row.get("ref") and row.get("alt")
        ]


def parse_info(info_field: str) -> dict[str, str]:
    info = {}
    if info_field in ("", "."):
        return info
    for entry in info_field.split(";"):
        if "=" in entry:
            key, value = entry.split("=", 1)
            info[key] = value
    return info


def first_number(value: str | None) -> str:
    if value is None or value in ("", "."):
        return ""
    first = value.split(",", maxsplit=1)[0]
    return "" if first in ("", ".") else first


def parse_float(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.6g}"


def sum_number_list(value: str | None) -> float | None:
    if value is None or value in ("", "."):
        return None
    total = 0.0
    found = False
    for part in value.split(","):
        if part in ("", "."):
            continue
        parsed = parse_float(part)
        if parsed is None:
            continue
        total += parsed
        found = True
    return total if found else None


def sample_metrics(
    fields: list[str],
    alt_index: int,
    sample_names: list[str],
    preferred_sample: str | None = None,
) -> tuple[str, str]:
    info = parse_info(fields[7])
    info_af = first_number(info.get("AF"))
    info_depth = first_number(info.get("BAM_DEPTH")) or first_number(info.get("DP"))

    if len(fields) < 10:
        return info_af, info_depth

    format_keys = fields[8].split(":")
    samples = fields[9:]
    sample_index = 0
    if preferred_sample and preferred_sample in sample_names:
        sample_index = sample_names.index(preferred_sample)
    elif "TUMOR" in sample_names:
        sample_index = sample_names.index("TUMOR")
    elif len(samples) > 1:
        sample_index = len(samples) - 1

    format_map = dict(zip(format_keys, samples[sample_index].split(":")))
    af = first_number(format_map.get("AF")) or info_af
    depth = first_number(format_map.get("DP")) or info_depth

    if not af and "AD" in format_map:
        counts = [parse_float(part) for part in format_map["AD"].split(",")]
        if len(counts) > alt_index and counts[0] is not None and counts[alt_index] is not None:
            total = counts[0] + counts[alt_index]
            af = format_float(None if total == 0 else counts[alt_index] / total)
            depth = depth or format_float(total)

    if not af:
        ref_count = sum_number_list(format_map.get(f"{fields[3]}U"))
        alt = fields[4].split(",")[alt_index - 1]
        alt_count = sum_number_list(format_map.get(f"{alt}U"))
        if ref_count is not None and alt_count is not None:
            total = ref_count + alt_count
            af = format_float(None if total == 0 else alt_count / total)
            depth = depth or format_float(total)

    return af, depth


def collect_hits(path: Path, caller: str, wanted: set[tuple[str, str, str, str]]) -> dict[tuple[str, str, str, str], Hit]:
    hits: dict[tuple[str, str, str, str], Hit] = {}
    sample_names: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, ref = normalize_chrom(fields[0]), fields[1], fields[3].upper()
            alts = [alt.upper() for alt in fields[4].split(",")]
            for index, alt in enumerate(alts, start=1):
                key = (chrom, pos, ref, alt)
                if key not in wanted:
                    continue
                af, depth = sample_metrics(fields, index, sample_names)
                hit = Hit(caller=caller, filter=fields[6], af=af, depth=depth)
                previous = hits.get(key)
                if previous is None or (not previous.passable and hit.passable):
                    hits[key] = hit
    return hits


def choose_metric(hits: Iterable[Hit], attr: str) -> str:
    for hit in hits:
        value = getattr(hit, attr)
        if value:
            return value
    return ""


def classify(intersect: Hit | None, mutect2: Hit | None, strelka: Hit | None) -> str:
    if intersect is not None and intersect.passable:
        return "HIGH_CONFIDENCE"
    if intersect is not None:
        return "INTERSECT_NONPASS"
    if mutect2 is not None and strelka is not None:
        return "BOTH_CALLERS_LOW_CONFIDENCE"
    if strelka is not None:
        return "STRELKA_ONLY"
    if mutect2 is not None:
        return "MUTECT2_ONLY"
    return "ABSENT"


def classify_final(final: Hit | None) -> str:
    if final is not None and final.passable:
        return "HIGH_CONFIDENCE"
    if final is not None:
        return "INTERSECT_NONPASS"
    return "ABSENT"


def build_row(
    tumour: str,
    hotspot: Hotspot,
    intersect: Hit | None,
    mutect2: Hit | None,
    strelka: Hit | None,
) -> dict[str, str]:
    ordered_hits = [hit for hit in (intersect, mutect2, strelka) if hit is not None]
    return {
        "tumour": tumour,
        "hotspot": hotspot.name,
        "chrom": hotspot.chrom,
        "pos": hotspot.pos,
        "ref": hotspot.ref,
        "alt": hotspot.alt,
        "status": classify(intersect, mutect2, strelka),
        "callers": ";".join(hit.caller for hit in ordered_hits),
        "filters": ";".join(f"{hit.caller}:{hit.filter or '.'}" for hit in ordered_hits),
        "af": choose_metric(ordered_hits, "af"),
        "depth": choose_metric(ordered_hits, "depth"),
    }


def build_final_only_row(tumour: str, hotspot: Hotspot, final: Hit | None) -> dict[str, str]:
    hits = [final] if final is not None else []
    return {
        "tumour": tumour,
        "hotspot": hotspot.name,
        "chrom": hotspot.chrom,
        "pos": hotspot.pos,
        "ref": hotspot.ref,
        "alt": hotspot.alt,
        "status": classify_final(final),
        "callers": ";".join(hit.caller for hit in hits),
        "filters": ";".join(f"{hit.caller}:{hit.filter or '.'}" for hit in hits),
        "af": choose_metric(hits, "af"),
        "depth": choose_metric(hits, "depth"),
    }


def parse_sample_inputs(values: list[str]) -> list[SampleInput]:
    sample_inputs = []
    for value in values:
        fields = value.split("=", maxsplit=3)
        if len(fields) == 2 and fields[0] and fields[1]:
            sample_inputs.append(SampleInput(tumour=fields[0], intersect=Path(fields[1])))
            continue
        if len(fields) != 4:
            raise ValueError(
                "--input values must be tumour=final_vcf or tumour=intersect_vcf=mutect2_vcf=strelka_vcf"
            )
        sample_inputs.append(
            SampleInput(
                tumour=fields[0],
                intersect=Path(fields[1]),
                mutect2=Path(fields[2]),
                strelka=Path(fields[3]),
            )
        )
    return sample_inputs


def write_rows(rows: Iterable[dict[str, str]], output: Path) -> None:
    fieldnames = ["tumour", "hotspot", "chrom", "pos", "ref", "alt", "status", "callers", "filters", "af", "depth"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate exact hotspot calls.")
    parser.add_argument("--hotspots", required=True, type=Path)
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    hotspots = read_hotspots(args.hotspots)
    wanted = {hotspot.key for hotspot in hotspots}
    rows = []
    for sample in parse_sample_inputs(args.input):
        if sample.mutect2 is None and sample.strelka is None:
            final_hits = collect_hits(sample.intersect, "final", wanted)
            for hotspot in hotspots:
                rows.append(build_final_only_row(sample.tumour, hotspot, final_hits.get(hotspot.key)))
            continue

        intersect_hits = collect_hits(sample.intersect, "intersect", wanted)
        mutect2_hits = collect_hits(sample.mutect2, "mutect2", wanted)
        strelka_hits = collect_hits(sample.strelka, "strelka", wanted)
        for hotspot in hotspots:
            key = hotspot.key
            rows.append(
                build_row(
                    sample.tumour,
                    hotspot,
                    intersect_hits.get(key),
                    mutect2_hits.get(key),
                    strelka_hits.get(key),
                )
            )
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
