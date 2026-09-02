"""Regression tests for mosdepth WGS coverage summarization."""

from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/wgs_coverage_mqc.py"


def write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_mosdepth_threshold_header_is_accepted(tmp_path: Path) -> None:
    regions = tmp_path / "sample.regions.bed.gz"
    thresholds = tmp_path / "sample.thresholds.bed.gz"
    output = tmp_path / "coverage.tsv"
    write_gzip(regions, "1\t0\t100\tCALLABLE\t6\nX\t0\t20\tCALLABLE\t4\n")
    write_gzip(
        thresholds,
        "#chrom\tstart\tend\tregion\t1X\t3X\t5X\n"
        "1\t0\t100\tCALLABLE\t100\t90\t70\n"
        "X\t0\t20\tCALLABLE\t20\t15\t5\n",
    )

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--sample",
            "sample",
            "--role",
            "baseline",
            "--expected-depth",
            "6",
            "--regions",
            str(regions),
            "--thresholds",
            str(thresholds),
            "--output",
            str(output),
        ],
        check=True,
    )

    rows = output.read_text().splitlines()
    assert "pct_1x\tpct_3x\tpct_5x" in rows[-2]
    assert rows[-1].startswith("sample\tbaseline\t6.0\t6.0\tOK\t100.0\t87.5\t62.5\t")


def test_wgs_rule_does_not_move_mosdepth_outputs_onto_themselves() -> None:
    rules = (ROOT / "workflow/rules/qc.smk").read_text()
    rule = rules.split("rule wgs_per_contig_coverage:", 1)[1].split(
        "rule prepare_exon_coverage_bed:", 1
    )[0]
    assert "mv {params.prefix}.regions.bed.gz {output.regions}" not in rule
    assert "mv {params.prefix}.thresholds.bed.gz {output.thresholds}" not in rule
    assert "test -s {output.regions}" in rule
    assert "test -s {output.thresholds}" in rule
