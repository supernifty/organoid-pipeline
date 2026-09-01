#!/usr/bin/env python3
"""Validate hotspot aggregate classification."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def write_hotspots(tmp: Path) -> Path:
    path = tmp / "hotspots.tsv"
    path.write_text(
        "\n".join(
            [
                "chrom\tpos\tref\talt\tname",
                "1\t100\tA\tC\thigh",
                "1\t101\tA\tG\tnonpass",
                "1\t102\tA\tT\tstrelka_only",
                "1\t103\tC\tT\tmutect2_only",
                "1\t104\tG\tA\tboth_low",
                "1\t105\tT\tG\tabsent",
                "",
            ]
        )
    )
    return path


def write_vcf(tmp: Path, name: str, records: list[tuple[int, str, str, str, str, str]]) -> Path:
    path = tmp / name
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR",
    ]
    for pos, ref, alt, filter_value, af, depth in records:
        alt_count = int(round(float(depth) * float(af)))
        ref_count = int(depth) - alt_count
        lines.append(
            f"chr1\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_value}\t."
            f"\tAD:AF:DP\t{ref_count},{alt_count}:{af}:{depth}"
        )
    lines.append("")
    path.write_text("\n".join(lines))
    return path


def main() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        hotspots = write_hotspots(tmp)
        intersect = write_vcf(
            tmp,
            "tumour.intersect.vcf",
            [
                (100, "A", "C", "PASS", "0.10", "100"),
                (101, "A", "G", "artifact", "0.20", "90"),
            ],
        )
        mutect2 = write_vcf(
            tmp,
            "tumour.mutect2.vcf",
            [
                (100, "A", "C", "PASS", "0.10", "100"),
                (101, "A", "G", "PASS", "0.20", "90"),
                (103, "C", "T", "PASS", "0.30", "80"),
                (104, "G", "A", "weak_evidence", "0.40", "70"),
            ],
        )
        strelka = write_vcf(
            tmp,
            "tumour.strelka.vcf",
            [
                (100, "A", "C", "PASS", "0.10", "100"),
                (102, "A", "T", "PASS", "0.25", "60"),
                (104, "G", "A", "LowDepth", "0.40", "70"),
            ],
        )
        output = tmp / "hotspots.out.tsv"

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "hotspots.py"),
                "--hotspots",
                str(hotspots),
                "--input",
                f"tumour={intersect}={mutect2}={strelka}",
                "--output",
                str(output),
            ],
            cwd=ROOT,
            check=True,
        )

        with output.open(newline="") as handle:
            rows = {
                row["hotspot"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }

        expected_statuses = {
            "high": "HIGH_CONFIDENCE",
            "nonpass": "INTERSECT_NONPASS",
            "strelka_only": "STRELKA_ONLY",
            "mutect2_only": "MUTECT2_ONLY",
            "both_low": "BOTH_CALLERS_LOW_CONFIDENCE",
            "absent": "ABSENT",
        }
        observed_statuses = {name: rows[name]["status"] for name in expected_statuses}
        if observed_statuses != expected_statuses:
            raise AssertionError(f"unexpected statuses: {observed_statuses}")

        if rows["high"]["callers"] != "intersect;mutect2;strelka":
            raise AssertionError(f"unexpected high-confidence callers: {rows['high']}")
        if rows["absent"]["callers"] or rows["absent"]["filters"]:
            raise AssertionError(f"absent hotspot should have empty evidence: {rows['absent']}")

        final_output = tmp / "hotspots.final.out.tsv"
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "workflow" / "scripts" / "hotspots.py"),
                "--hotspots",
                str(hotspots),
                "--input",
                f"tumour={intersect}",
                "--output",
                str(final_output),
            ],
            cwd=ROOT,
            check=True,
        )

        with final_output.open(newline="") as handle:
            final_rows = {
                row["hotspot"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }

        if final_rows["high"]["status"] != "HIGH_CONFIDENCE":
            raise AssertionError(f"unexpected final-only high row: {final_rows['high']}")
        if final_rows["high"]["callers"] != "final":
            raise AssertionError(f"unexpected final-only caller: {final_rows['high']}")
        if final_rows["nonpass"]["status"] != "INTERSECT_NONPASS":
            raise AssertionError(f"unexpected final-only nonpass row: {final_rows['nonpass']}")
        if final_rows["absent"]["status"] != "ABSENT":
            raise AssertionError(f"unexpected final-only absent row: {final_rows['absent']}")


if __name__ == "__main__":
    main()
