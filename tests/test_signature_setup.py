#!/usr/bin/env python3
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "config.yaml"
FILTER_SCRIPT = ROOT / "workflow" / "scripts" / "filter_signature_definitions.py"


def load_config():
    with open(CONFIG_PATH, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def read_definition_signatures(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "Sig":
            raise AssertionError(f"Definition file has unexpected header: {path}")
        return {row[0] for row in reader if row}


def assert_exists(path):
    if not path.exists():
        raise AssertionError(f"Missing required file: {path}")


def validate_configured_signatures(config):
    signatures_cfg = config["mutational_signatures"]
    tool_root = ROOT / "tools" / "mutational_signature" / "mutational_signature"

    required_files = [
        tool_root / "__init__.py",
        tool_root / "count.py",
        tool_root / "decompose.py",
        ROOT / signatures_cfg["sbs_definition"],
        ROOT / signatures_cfg["id_definition"],
        ROOT / signatures_cfg["dbs_definition"],
    ]
    for path in required_files:
        assert_exists(path)

    for signature_type in ("sbs", "id", "dbs"):
        definition_path = ROOT / signatures_cfg[f"{signature_type}_definition"]
        available = read_definition_signatures(definition_path)
        missing = [
            name
            for name in signatures_cfg["tissue_signatures"][signature_type]
            if name not in available
        ]
        if missing:
            raise AssertionError(
                f"Configured {signature_type} signatures missing from {definition_path}: "
                + ", ".join(missing)
            )


def run_filter_helper_tests():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        definition = tmpdir_path / "definitions.tsv"
        success_output = tmpdir_path / "filtered.tsv"

        definition.write_text(
            "Sig\tCTX1\tCTX2\nSIG_A\t0.1\t0.9\nSIG_B\t0.3\t0.7\n",
            encoding="utf-8",
        )

        success = subprocess.run(
            [
                sys.executable,
                str(FILTER_SCRIPT),
                "--definition",
                str(definition),
                "--output",
                str(success_output),
                "--signatures",
                "SIG_B",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if success.returncode != 0:
            raise AssertionError(f"Helper success case failed: {success.stderr}")

        rows = success_output.read_text(encoding="utf-8").strip().splitlines()
        expected = ["Sig\tCTX1\tCTX2", "SIG_B\t0.3\t0.7"]
        if rows != expected:
            raise AssertionError(f"Unexpected helper output: {rows}")

        failure = subprocess.run(
            [
                sys.executable,
                str(FILTER_SCRIPT),
                "--definition",
                str(definition),
                "--output",
                str(tmpdir_path / "missing.tsv"),
                "--signatures",
                "SIG_MISSING",
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if failure.returncode == 0:
            raise AssertionError("Helper missing-signature case unexpectedly succeeded")
        if "Configured signatures missing" not in failure.stderr:
            raise AssertionError(
                "Helper missing-signature case did not report the expected error"
            )


def main():
    config = load_config()
    validate_configured_signatures(config)
    run_filter_helper_tests()
    print("Signature setup tests passed.")


if __name__ == "__main__":
    main()
