#!/usr/bin/env python3
"""Validate provenance TSV generation."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["key", "value", "source", "note"]:
            raise AssertionError(f"unexpected header: {reader.fieldnames}")
        return {row["key"]: row for row in reader}


def main() -> None:
    samples = ROOT / "config" / "samples.yaml"
    if not samples.exists():
        samples = ROOT / "config" / "samples.example.yaml"

    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "provenance.tsv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "write_provenance.py"),
                "--config",
                str(ROOT / "config" / "config.yaml"),
                "--samples",
                str(samples),
                "--pixi",
                str(ROOT / "pixi.toml"),
                "--slurm-config",
                str(ROOT / "config" / "slurm" / "config.yaml"),
                "--output",
                str(output),
            ],
            cwd=ROOT,
            check=True,
        )

        text = output.read_text()
        fastq_markers = ("fastq_1", "fastq_2", "_R1.fastq", "_R2.fastq", "_R1.fq", "_R2.fq")
        if any(marker in text for marker in fastq_markers):
            raise AssertionError("provenance output includes FASTQ details")

        rows = load_rows(output)
        required_keys = {
            "run.path",
            "pipeline.name",
            "pipeline.version",
            "analysis.type",
            "samples.ids",
            "samples.tumour_normal_pairs",
            "samples.matched_normal_ids",
            "reference.genome",
            "container_runtime",
            "mutect2.interval_padding",
            "germline.hard_filters.snp.qd_min",
            "annotation.enabled",
            "coverage.exon_enabled",
            "somalier.minimum_depth",
            "filtering.af_threshold",
            "software.gatk.version",
            "software.bwa.declared_dependency",
        }
        missing = sorted(required_keys - set(rows))
        if missing:
            raise AssertionError(f"missing provenance keys: {missing}")

        for key in ("samples.ids", "samples.tumour_normal_pairs", "chromosomes"):
            json.loads(rows[key]["value"])

        overlay = Path(tmpdir) / "config.local.yaml"
        overlay.write_text("analysis:\n  type: wgs\n")
        overlay_output = Path(tmpdir) / "provenance.overlay.tsv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "write_provenance.py"),
                "--config",
                str(ROOT / "config" / "config.yaml"),
                "--config-overlay",
                str(overlay),
                "--samples",
                str(samples),
                "--pixi",
                str(ROOT / "pixi.toml"),
                "--output",
                str(overlay_output),
            ],
            cwd=ROOT,
            check=True,
        )
        overlay_rows = load_rows(overlay_output)
        assert overlay_rows["analysis.type"]["value"] == "wgs"
        assert overlay_rows["provenance.config_overlay"]["value"] == str(overlay)


if __name__ == "__main__":
    main()
