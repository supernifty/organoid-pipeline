#!/usr/bin/env python3
"""Exercise capture remapping publication and QC with a deterministic liftOver fixture."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    (ROOT / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ROOT / "tmp") as directory:
        root = Path(directory)
        source = root / "capture.hg19.bed"
        chain = root / "fixture.chain"
        fai = root / "target.fa.fai"
        output = root / "capture.grch38.bed.gz"
        mapping = root / "capture.mapping.bed"
        unmapped = root / "capture.unmapped.bed"
        report = root / "capture.json"
        fake_liftover = root / "liftOver"

        source.write_text("1\t10\t20\n1\t20\t30\nX\t5\t10\n")
        chain.write_text("fixture\n")
        fai.write_text("chr1\t1000\t0\t0\t0\nchrX\t500\t0\t0\t0\n")
        fake_liftover.write_text(
            "#!/bin/sh\n"
            "cp \"$1\" \"$3\"\n"
            ": > \"$4\"\n"
        )
        fake_liftover.chmod(0o755)

        subprocess.run(
            [
                sys.executable,
                str(ROOT / "scripts" / "remap_capture_bed.py"),
                "--source-bed", str(source),
                "--chain", str(chain),
                "--target-fai", str(fai),
                "--output", str(output),
                "--mapping", str(mapping),
                "--unmapped", str(unmapped),
                "--report", str(report),
                "--liftover", str(fake_liftover),
                "--bgzip", str(Path(sys.executable).with_name("bgzip")),
                "--tabix", str(Path(sys.executable).with_name("tabix")),
            ],
            cwd=ROOT,
            env={**os.environ, "TMPDIR": str(root)},
            check=True,
        )
        details = json.loads(report.read_text())
        assert details["mapped_base_fraction"] == 1.0
        assert details["source_interval_count"] == 3
        assert details["canonical_interval_count"] == 2
        assert output.is_file() and Path(f"{output}.tbi").is_file()
        assert len(mapping.read_text().splitlines()) == 3
        assert unmapped.read_text() == ""


if __name__ == "__main__":
    main()
