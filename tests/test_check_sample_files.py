#!/usr/bin/env python3
"""Test sample file existence checks."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.check_sample_files import find_cram_alternatives
from workflow.scripts.sample_inputs import missing_sample_files


def touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
    return str(path)


def test_missing_sample_files() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        samples = {
            "samples": {
                "FASTQ_OK": {
                    "fastq_1": touch(tmpdir / "FASTQ_OK_R1.fastq.gz"),
                    "fastq_2": touch(tmpdir / "FASTQ_OK_R2.fastq.gz"),
                },
                "CRAM_OK": {
                    "cram": touch(tmpdir / "CRAM_OK.sorted.dups.cram"),
                },
                "MISSING": {
                    "cram": str(tmpdir / "MISSING.sorted.dups.cram"),
                    "crai": str(tmpdir / "MISSING.sorted.dups.cram.crai"),
                },
            },
            "tumours": {"CRAM_OK": "FASTQ_OK"},
        }
        touch(tmpdir / "CRAM_OK.sorted.dups.cram.crai")

        assert missing_sample_files(samples, "test samples") == [
            ("MISSING", "cram", str(tmpdir / "MISSING.sorted.dups.cram")),
            ("MISSING", "crai", str(tmpdir / "MISSING.sorted.dups.cram.crai")),
        ]


def test_check_sample_files_cli() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        r1 = touch(tmpdir / "S_R1.fastq.gz")
        r2 = touch(tmpdir / "S_R2.fastq.gz")
        samples_yaml = tmpdir / "samples.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    "  S:",
                    f"    fastq_1: {r1}",
                    f"    fastq_2: {r2}",
                    "tumours: {}",
                    "",
                ]
            )
        )

        result = subprocess.run(
            [sys.executable, "scripts/check_sample_files.py", str(samples_yaml)],
            cwd=ROOT,
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 0
        assert "All sample files exist" in result.stdout


def test_find_cram_alternatives() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        missing = tmpdir / "0151020010_38190_BC.sorted.dups.cram"
        expected = touch(tmpdir / "0151020010_2983_BC.sorted.dups.cram")
        touch(tmpdir / "0151020010_38190_P1.sorted.dups.cram")
        touch(tmpdir / "0151020011_2983_BC.sorted.dups.cram")

        assert find_cram_alternatives(str(missing)) == [expected]


def test_find_cram_alternatives_across_batch_dirs() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        root = (
            tmpdir / "data" / "projects" / "project123" / "researcher" / "analyses" / "somatic-wes"
        )
        missing = root / "AGRF_CAGRF24030235_22KLWWLT3" / "out" / "0656012009_BC.sorted.dups.cram"
        expected = touch(
            root / "AGRF_CAGRF230514735_HC7HWDSX7" / "out" / "0656012009_1484_BC.sorted.dups.cram"
        )

        assert find_cram_alternatives(str(missing)) == [expected]


def test_find_cram_alternatives_hyphenated_prefix() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        root = (
            tmpdir / "data" / "projects" / "project123" / "researcher" / "analyses" / "somatic-wes"
        )
        missing = root / "PeterMac-WES-20231120" / "out" / "AUS509002-39399_BC.sorted.dups.cram"
        expected = touch(
            root / "PeterMac-WES-20231120" / "out" / "AUS509002-39448_BC.sorted.dups.cram"
        )

        assert find_cram_alternatives(str(missing)) == [expected]


def test_find_cram_alternatives_expands_middle_token() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        missing = tmpdir / "0656012009_P1.sorted.dups.cram"
        expected = touch(tmpdir / "0656012009_39233_P1.sorted.dups.cram")

        assert find_cram_alternatives(str(missing)) == [expected]


def test_find_cram_alternatives_leading_zero_variant() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        missing = tmpdir / "151007097_BC.sorted.dups.cram"
        expected = touch(tmpdir / "0151007097_BC.sorted.dups.cram")

        assert find_cram_alternatives(str(missing)) == [expected]


def test_find_cram_alternatives_reduces_middle_token() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        missing = tmpdir / "AUS469001_39277_BC.sorted.dups.cram"
        expected = touch(tmpdir / "AUS469001_BC.sorted.dups.cram")

        assert find_cram_alternatives(str(missing)) == [expected]


def test_check_sample_files_cli_writes_accepted_fix() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        tumour_cram = touch(outdir / "0151020010_38190_P1.sorted.dups.cram")
        touch(outdir / "0151020010_38190_P1.sorted.dups.cram.crai")
        normal_cram = touch(outdir / "0151020010_2983_BC.sorted.dups.cram")
        touch(outdir / "0151020010_2983_BC.sorted.dups.cram.crai")
        missing_normal = outdir / "0151020010_38190_BC.sorted.dups.cram"
        samples_yaml = tmpdir / "samples.yaml"
        fixed_yaml = tmpdir / "samples.fixed.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    '  "0151020010_38190_P1":',
                    f"    cram: {tumour_cram}",
                    '  "0151020010_38190_BC":',
                    f"    cram: {missing_normal}",
                    "tumours:",
                    '  "0151020010_38190_P1": "0151020010_38190_BC"',
                    "",
                ]
            )
        )

        result = subprocess.run(
            [
                sys.executable,
                "scripts/check_sample_files.py",
                str(samples_yaml),
                "--fix-output",
                str(fixed_yaml),
            ],
            cwd=ROOT,
            input="y\n",
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 0
        assert fixed_yaml.exists()
        assert str(normal_cram) in fixed_yaml.read_text()
        original = yaml.safe_load(samples_yaml.read_text())
        fixed = yaml.safe_load(fixed_yaml.read_text())
        assert "0151020010_38190_BC" in original["samples"]
        assert "0151020010_38190_BC" not in fixed["samples"]
        assert fixed["samples"]["0151020010_2983_BC"]["cram"] == str(normal_cram)
        assert fixed["tumours"]["0151020010_38190_P1"] == "0151020010_2983_BC"


def test_check_sample_files_cli_merges_existing_normal() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        tumour_1_cram = touch(outdir / "0151020010_38190_P1.sorted.dups.cram")
        touch(outdir / "0151020010_38190_P1.sorted.dups.cram.crai")
        tumour_2_cram = touch(outdir / "0151020010_38191_P1.sorted.dups.cram")
        touch(outdir / "0151020010_38191_P1.sorted.dups.cram.crai")
        normal_cram = touch(outdir / "0151020010_2983_BC.sorted.dups.cram")
        touch(outdir / "0151020010_2983_BC.sorted.dups.cram.crai")
        missing_normal = outdir / "0151020010_38190_BC.sorted.dups.cram"
        samples_yaml = tmpdir / "samples.yaml"
        fixed_yaml = tmpdir / "samples.fixed.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    '  "0151020010_38190_P1":',
                    f"    cram: {tumour_1_cram}",
                    '  "0151020010_38191_P1":',
                    f"    cram: {tumour_2_cram}",
                    '  "0151020010_2983_BC":',
                    f"    cram: {normal_cram}",
                    '  "0151020010_38190_BC":',
                    f"    cram: {missing_normal}",
                    "tumours:",
                    '  "0151020010_38190_P1": "0151020010_38190_BC"',
                    '  "0151020010_38191_P1": "0151020010_2983_BC"',
                    "",
                ]
            )
        )

        result = subprocess.run(
            [
                sys.executable,
                "scripts/check_sample_files.py",
                str(samples_yaml),
                "--fix-output",
                str(fixed_yaml),
            ],
            cwd=ROOT,
            input="y\n",
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 0
        fixed = yaml.safe_load(fixed_yaml.read_text())
        assert "0151020010_38190_BC" not in fixed["samples"]
        assert fixed["samples"]["0151020010_2983_BC"]["cram"] == str(normal_cram)
        assert fixed["tumours"]["0151020010_38190_P1"] == "0151020010_2983_BC"
        assert fixed["tumours"]["0151020010_38191_P1"] == "0151020010_2983_BC"


def test_check_sample_files_cli_keeps_existing_normal_match() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        root = (
            tmpdir / "data" / "projects" / "project123" / "researcher" / "analyses" / "somatic-wes"
        )
        existing_normal = touch(
            root / "AGRF_CAGRF230514735_HC7HWDSX7" / "out" / "AUS469001_BC.sorted.dups.cram"
        )
        touch(root / "AGRF_CAGRF230514735_HC7HWDSX7" / "out" / "AUS469001_BC.sorted.dups.cram.crai")
        touch(root / "AGRF_CAGRF230715533_HJGVMDSX7" / "out" / "AUS469001_BC.sorted.dups.cram")
        touch(root / "AGRF_CAGRF230715533_HJGVMDSX7" / "out" / "AUS469001_BC.sorted.dups.cram.crai")
        missing_normal = (
            root / "AGRF_CAGRF24030235_22KLWWLT3" / "out" / "AUS469001_39277_BC.sorted.dups.cram"
        )
        tumour_cram = touch(
            root / "AGRF_CAGRF24030235_22KLWWLT3" / "out" / "AUS469001_P5.sorted.dups.cram"
        )
        touch(root / "AGRF_CAGRF24030235_22KLWWLT3" / "out" / "AUS469001_P5.sorted.dups.cram.crai")
        samples_yaml = tmpdir / "samples.yaml"
        fixed_yaml = tmpdir / "samples.fixed.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    '  "AUS469001_P5":',
                    f"    cram: {tumour_cram}",
                    '  "AUS469001_39277_BC":',
                    f"    cram: {missing_normal}",
                    '  "AUS469001_BC":',
                    f"    cram: {existing_normal}",
                    "tumours:",
                    '  "AUS469001_P5": "AUS469001_39277_BC"',
                    "",
                ]
            )
        )

        result = subprocess.run(
            [
                sys.executable,
                "scripts/check_sample_files.py",
                str(samples_yaml),
                "--fix-output",
                str(fixed_yaml),
            ],
            cwd=ROOT,
            input="2\n",
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 0
        assert "WARNING: cannot merge" not in result.stdout
        fixed = yaml.safe_load(fixed_yaml.read_text())
        assert "AUS469001_39277_BC" not in fixed["samples"]
        assert fixed["samples"]["AUS469001_BC"]["cram"] == existing_normal
        assert fixed["tumours"]["AUS469001_P5"] == "AUS469001_BC"


def test_check_sample_files_cli_renames_tumour_key() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        tumour_cram = touch(outdir / "0656012009_39233_P1.sorted.dups.cram")
        touch(outdir / "0656012009_39233_P1.sorted.dups.cram.crai")
        normal_cram = touch(outdir / "0656012009_BC.sorted.dups.cram")
        touch(outdir / "0656012009_BC.sorted.dups.cram.crai")
        missing_tumour = outdir / "0656012009_P1.sorted.dups.cram"
        samples_yaml = tmpdir / "samples.yaml"
        fixed_yaml = tmpdir / "samples.fixed.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    '  "0656012009_P1":',
                    f"    cram: {missing_tumour}",
                    '  "0656012009_BC":',
                    f"    cram: {normal_cram}",
                    "tumours:",
                    '  "0656012009_P1": "0656012009_BC"',
                    "",
                ]
            )
        )

        result = subprocess.run(
            [
                sys.executable,
                "scripts/check_sample_files.py",
                str(samples_yaml),
                "--fix-output",
                str(fixed_yaml),
            ],
            cwd=ROOT,
            input="y\n",
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 0
        fixed = yaml.safe_load(fixed_yaml.read_text())
        assert "0656012009_P1" not in fixed["samples"]
        assert fixed["samples"]["0656012009_39233_P1"]["cram"] == str(tumour_cram)
        assert fixed["tumours"]["0656012009_39233_P1"] == "0656012009_BC"


def test_check_sample_files_cli_skip_writes_no_output() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        outdir = tmpdir / "out"
        touch(outdir / "0151020010_2983_BC.sorted.dups.cram")
        samples_yaml = tmpdir / "samples.yaml"
        fixed_yaml = tmpdir / "samples.fixed.yaml"
        samples_yaml.write_text(
            "\n".join(
                [
                    "samples:",
                    '  "0151020010_38190_BC":',
                    f"    cram: {outdir / '0151020010_38190_BC.sorted.dups.cram'}",
                    "tumours: {}",
                    "",
                ]
            )
        )

        result = subprocess.run(
            [
                sys.executable,
                "scripts/check_sample_files.py",
                str(samples_yaml),
                "--fix-output",
                str(fixed_yaml),
            ],
            cwd=ROOT,
            input="\n",
            check=False,
            text=True,
            capture_output=True,
        )

        assert result.returncode == 1
        assert not fixed_yaml.exists()
        assert "No fixes accepted" in result.stdout


if __name__ == "__main__":
    test_missing_sample_files()
    test_check_sample_files_cli()
    test_find_cram_alternatives()
    test_find_cram_alternatives_across_batch_dirs()
    test_find_cram_alternatives_hyphenated_prefix()
    test_find_cram_alternatives_expands_middle_token()
    test_find_cram_alternatives_leading_zero_variant()
    test_find_cram_alternatives_reduces_middle_token()
    test_check_sample_files_cli_writes_accepted_fix()
    test_check_sample_files_cli_merges_existing_normal()
    test_check_sample_files_cli_keeps_existing_normal_match()
    test_check_sample_files_cli_renames_tumour_key()
    test_check_sample_files_cli_skip_writes_no_output()
