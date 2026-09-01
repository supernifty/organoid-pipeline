#!/usr/bin/env python3
"""Tests for exon resource validation and mosdepth arithmetic."""

from __future__ import annotations

import csv
import gzip
import hashlib
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from coverage_qc import prepare_bed, summarize, validate_exon_config  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        source = root / "exons.v1.bed"
        source.write_text(
            "CHROM\tSTART\tEND\tEXON_ID\tGENE_ID\tGENE_SYMBOL\n"
            "1\t0\t10\tE1\tG1\tGENE1\n1\t10\t20\tE2\tG1\tGENE1\n"
        )
        config = {"coverage": {"exon_enabled": True, "exon_bed": {
            "build": "grch37", "version": "v1", "path": str(source),
            "expected_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "source": "test", "licence": "test", "access_date": "2026-08-25",
        }}}
        validate_exon_config(config, "grch37")
        bed = root / "prepared.bed"
        prepare_bed(source, bed)
        regions, thresholds = root / "sample.regions.bed.gz", root / "sample.thresholds.bed.gz"
        with gzip.open(regions, "wt") as handle:
            handle.write("1\t0\t10\tE1|G1|GENE1\t20\n1\t10\t20\tE2|G1|GENE1\t0\n")
        with gzip.open(thresholds, "wt") as handle:
            handle.write("#chrom\tstart\tend\tregion\t10X\t20X\t50X\t100X\n")
            handle.write("1\t0\t10\tE1|G1|GENE1\t10\t10\t0\t0\n")
            handle.write("1\t10\t20\tE2|G1|GENE1\t0\t0\t0\t0\n")
        exon, gene = root / "exon.tsv", root / "gene.tsv"
        summarize("T", regions, thresholds, exon, gene, 20)
        exon_rows = list(csv.DictReader(exon.open(), delimiter="\t"))
        gene_row = next(csv.DictReader(gene.open(), delimiter="\t"))
        assert exon_rows[0]["complete_20x"] == "true" and exon_rows[1]["coverage_warning"] == "true"
        assert gene_row["mean_depth"] == "10.000000"
        assert gene_row["covered_fraction_20x"] == "0.500000" and gene_row["complete_20x"] == "false"


if __name__ == "__main__":
    main()
