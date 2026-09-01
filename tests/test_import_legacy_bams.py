#!/usr/bin/env python3
"""Unit tests for legacy BAM importer validation and refusal behavior."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import import_legacy_bams as importer  # noqa: E402


def expect_error(function, text: str) -> None:
    try:
        function()
    except importer.ImportError as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError(f"expected ImportError containing {text!r}")


def main() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        mapping = root / "mapping.yaml"
        mapping.write_text(yaml.safe_dump({"samples": {"T": {}, "N": {}}, "tumours": {"T": "N"}}))
        assert importer.load_legacy_samples(mapping) == (["N", "T"], {"T": "N"}) or importer.load_legacy_samples(mapping) == (["T", "N"], {"T": "N"})
        listing = root / "listing.yaml"
        listing.write_text(yaml.safe_dump({"samples": [{"sample": "T"}, "N"], "tumours": [{"tumour": "T", "normal": "N"}]}))
        assert importer.load_legacy_samples(listing) == (["T", "N"], {"T": "N"})
        duplicate = root / "duplicate.yaml"
        duplicate.write_text(yaml.safe_dump({"samples": ["T", "T"]}))
        expect_error(lambda: importer.load_legacy_samples(duplicate), "unique")

        bam = root / "S.sorted.dups.bam"
        bam.write_bytes(b"bam")
        expect_error(lambda: importer.resolve_bam_index(bam), "missing BAM index")
        first, second = Path(f"{bam}.bai"), bam.with_suffix(".bai")
        first.write_bytes(b"one")
        second.write_bytes(b"two")
        expect_error(lambda: importer.resolve_bam_index(bam), "conflicting")
        second.write_bytes(b"one")
        assert importer.resolve_bam_index(bam) == first

        header = importer.parse_header(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:249250621\n"
            "@RG\tID:rg1\tSM:sample\tPL:ILLUMINA\n"
        )
        assert header["sample_name"] == "sample"
        expect_error(lambda: importer.parse_header("@HD\tSO:queryname\n@RG\tID:x\tSM:s\n"), "SO:coordinate")
        assert importer.idxstats_totals("1\t10\t3\t2\n*\t0\t0\t4\n") == {"mapped": 3, "unmapped": 6, "total": 9}
        mismatch = {"header": header, "totals": {"mapped": 1, "unmapped": 0, "total": 1}, "alignment_checksum": ["all\tall\t1\tx"]}
        expect_error(lambda: importer.validate_output(mismatch, {**mismatch, "totals": {"mapped": 2, "unmapped": 0, "total": 2}}), "totals differ")
        expect_error(lambda: importer.validate_output(mismatch, {**mismatch, "alignment_checksum": ["all\tall\t1\ty"]}), "checksums differ")

        legacy = root / "legacy"
        destination = root / "destination"
        legacy.mkdir()
        destination.mkdir()
        source_bam = legacy / "T.sorted.dups.bam"
        source_bam.write_bytes(b"bam")
        source_bam.with_suffix(".bai").write_bytes(b"index")
        reference = root / "reference.fa"
        reference.write_text(">1\nA\n")
        Path(f"{reference}.fai").write_text("1\t1\t0\t1\t2\n")
        (destination / "T.sorted.dups.cram").write_bytes(b"partial")
        partial_samples = root / "partial.yaml"
        partial_samples.write_text(yaml.safe_dump({"samples": ["T"]}))
        expect_error(lambda: importer.main([
            "--samples", str(partial_samples), "--legacy-output", str(legacy),
            "--destination", str(destination), "--reference", str(reference),
            "--output-samples", str(root / "modern.yaml"), "--allow-noncanonical-reference",
        ]), "partial destination")


if __name__ == "__main__":
    main()
