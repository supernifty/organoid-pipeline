#!/usr/bin/env python3
"""Count extended mutation contexts across tumour VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


VALID_BASES = set("ACGT")


@dataclass(frozen=True)
class FastaIndexRecord:
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFasta:
    def __init__(self, fasta: Path, fai: Path) -> None:
        self.fasta = fasta
        self.records = self._read_index(fai)
        self.handle = fasta.open("rb")

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "IndexedFasta":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    @staticmethod
    def _read_index(fai: Path) -> dict[str, FastaIndexRecord]:
        records: dict[str, FastaIndexRecord] = {}
        with fai.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                record = FastaIndexRecord(
                    name=fields[0],
                    length=int(fields[1]),
                    offset=int(fields[2]),
                    line_bases=int(fields[3]),
                    line_width=int(fields[4]),
                )
                for alias in chrom_aliases(fields[0]):
                    records.setdefault(alias, record)
        return records

    def fetch(self, chrom: str, start: int, end: int) -> str | None:
        """Return 1-based inclusive sequence, or None if the interval is invalid."""
        key = normalize_chrom(chrom)
        record = self.records.get(key)
        if record is None or start < 1 or end > record.length or start > end:
            return None

        chunks = []
        zero_based = start - 1
        remaining = end - start + 1
        while remaining > 0:
            line_offset = zero_based % record.line_bases
            bases_this_line = min(remaining, record.line_bases - line_offset)
            file_offset = (
                record.offset
                + (zero_based // record.line_bases) * record.line_width
                + line_offset
            )
            self.handle.seek(file_offset)
            chunks.append(self.handle.read(bases_this_line).decode("ascii"))
            zero_based += bases_this_line
            remaining -= bases_this_line
        return "".join(chunks).upper()


def normalize_chrom(chrom: str) -> str:
    normalized = chrom.split(".", maxsplit=1)[0].split("_", maxsplit=1)[-1]
    if normalized.lower().startswith("chr"):
        normalized = normalized[3:]
    normalized = normalized.upper()
    if normalized == "MT":
        return "M"
    return normalized


def chrom_aliases(chrom: str) -> set[str]:
    normalized = normalize_chrom(chrom)
    aliases = {normalized, chrom}
    if normalized == "M":
        aliases.add("MT")
    return {normalize_chrom(alias) for alias in aliases}


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def is_colibactin_aat_context(ref: str, context: str) -> bool:
    if ref == "T":
        return context[4] == "T" and context[0] == "A" and context[1] == "A"
    if ref == "A":
        return context[4] == "A" and context[7] == "T" and context[8] == "T"
    return False


def empty_counts(tumour: str) -> dict[str, int | str]:
    return {
        "tumour": tumour,
        "variant_count": 0,
        "snv_count": 0,
        "indel_count": 0,
        "colibactin_aat_count": 0,
        "colibactin_at_snv_count": 0,
    }


def count_vcf(tumour: str, vcf: Path, reference: IndexedFasta) -> dict[str, int | str]:
    counts = empty_counts(tumour)
    with open_text(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos_text, ref, alt_field = fields[0], fields[1], fields[3].upper(), fields[4].upper()
            pos = int(pos_text)
            for alt in alt_field.split(","):
                counts["variant_count"] += 1
                if len(ref) == 1 and len(alt) == 1:
                    if ref not in VALID_BASES or alt not in VALID_BASES:
                        continue
                    counts["snv_count"] += 1
                    if ref not in {"A", "T"}:
                        continue
                    context = reference.fetch(chrom, pos - 4, pos + 4)
                    if context is None or len(context) != 9 or set(context) - VALID_BASES:
                        continue
                    if context[4] != ref:
                        continue
                    counts["colibactin_at_snv_count"] += 1
                    if is_colibactin_aat_context(ref, context):
                        counts["colibactin_aat_count"] += 1
                else:
                    counts["indel_count"] += 1
    return counts


def add_proportion(row: dict[str, int | str]) -> dict[str, int | str]:
    denominator = int(row["colibactin_at_snv_count"])
    numerator = int(row["colibactin_aat_count"])
    row["colibactin_snv_proportion"] = f"{numerator / max(1, denominator):.6f}"
    return row


def write_counts(rows: Iterable[dict[str, int | str]], output: Path) -> None:
    fieldnames = [
        "tumour",
        "variant_count",
        "snv_count",
        "indel_count",
        "colibactin_aat_count",
        "colibactin_at_snv_count",
        "colibactin_snv_proportion",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(add_proportion(row))


def parse_inputs(values: list[str]) -> list[tuple[str, Path]]:
    inputs = []
    for value in values:
        if "=" in value:
            tumour, path = value.split("=", 1)
        else:
            path = value
            tumour = Path(path).name.split(".intersect.vcf", maxsplit=1)[0]
        inputs.append((tumour, Path(path)))
    return inputs


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate extended mutation context counts.")
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reference-index", required=True, type=Path)
    parser.add_argument("--input", required=True, nargs="+", help="VCF paths or tumour=VCF pairs")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    with IndexedFasta(args.reference, args.reference_index) as reference:
        rows = [
            count_vcf(tumour, vcf, reference)
            for tumour, vcf in parse_inputs(args.input)
        ]
    write_counts(rows, args.output)


if __name__ == "__main__":
    main()
