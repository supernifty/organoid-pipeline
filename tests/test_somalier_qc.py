#!/usr/bin/env python3
"""Tests for shared-normal groups and stable Somalier flags."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from somalier_qc import write_flags, write_groups  # noqa: E402


def main() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        samples = root / "samples.yaml"
        samples.write_text(yaml.safe_dump({"samples": {name: {} for name in ("T1", "T2", "N", "U")}, "tumours": {"T1": "N", "T2": "N"}}))
        groups = root / "groups.txt"
        write_groups(samples, groups)
        assert groups.read_text() == "N,T1,T2\n"
        pairs = root / "pairs.tsv"
        pairs.write_text(
            "#sample_a\tsample_b\trelatedness\n"
            "T1\tN\t0.70\nT2\tN\t0.50\nT1\tU\t0.45\nN\tU\t0.10\n"
        )
        output = root / "flags.tsv"
        write_flags(samples, pairs, output, 0.6, 0.4)
        rows = list(csv.DictReader(output.open(), delimiter="\t"))
        flags = {(row["sample_a"], row["sample_b"]): row["flag"] for row in rows}
        assert flags[("T1", "N")] == "pass"
        assert flags[("T2", "N")] == "expected_pair_mismatch"
        assert flags[("T1", "U")] == "unexpected_pair_possible_swap"


if __name__ == "__main__":
    main()
