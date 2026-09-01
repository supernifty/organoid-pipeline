#!/usr/bin/env python3
"""Tests for annotation validation, CSQ parsing, burden, and recurrence."""

from __future__ import annotations

import csv
import gzip
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from aggregate_annotations import write_burden, write_recurrence  # noqa: E402
from annotation import checksum_path, parse_csq, parse_csq_fields, select_pick, validate_annotation_config  # noqa: E402


def write_vcf(path: Path, records: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write('##fileformat=VCFv4.2\n')
        handle.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence|SYMBOL|Gene|Feature|PICK|CANONICAL">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        handle.writelines(record + "\n" for record in records)


def main() -> None:
    fields = parse_csq_fields('##INFO=<ID=CSQ,Description="Format: Allele|Consequence|Feature|PICK">')
    consequences = parse_csq("A|missense_variant|ENST2|,A|synonymous_variant|ENST1|1", fields)
    assert select_pick(consequences)["Feature"] == "ENST1"
    assert select_pick(list(reversed(consequences)))["Feature"] == "ENST1"
    assert parse_csq(None, fields) == []

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        resources = {}
        for name in ("cache", "plugins", "clinvar", "gnomad", "revel", "cadd_snv", "cadd_indel", "alphamissense", "spliceai_snv", "spliceai_indel"):
            path = root / name
            if name in ("cache", "plugins"):
                path.mkdir()
                (path / "release.txt").write_text("116\n")
            else:
                path.write_text(name)
            resource = {"build": "grch37", "version": "116", "path": str(path),
                        "expected_sha256": checksum_path(path), "source": "test", "licence": "test", "access_date": "2026-08-25"}
            if path.is_file():
                index = Path(f"{path}.tbi")
                index.write_text("index")
                resource["index"] = {"path": str(index), "expected_sha256": checksum_path(index)}
                os.utime(index, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 1))
            resources[name] = resource
        config = {"reference": {"build": "grch37"}, "annotation": {
            "enabled": True, "version": 116, "cache_version": 116, "plugin_release": 116,
            "docker_image": "ensemblorg/ensembl-vep@sha256:" + "a" * 64,
            "pick_order": ["canonical", "rank"], "resources": resources,
        }}
        assert len(validate_annotation_config(config, "grch37")) == 10
        resources["clinvar"]["expected_sha256"] = "0" * 64
        try:
            validate_annotation_config(config, "grch37")
        except ValueError as exc:
            assert "checksum mismatch" in str(exc)
        else:
            raise AssertionError("expected checksum mismatch")
        resources["clinvar"]["expected_sha256"] = checksum_path(Path(resources["clinvar"]["path"]))
        clinvar = Path(resources["clinvar"]["path"])
        clinvar_index = Path(resources["clinvar"]["index"]["path"])
        os.utime(clinvar, ns=(clinvar.stat().st_atime_ns, clinvar_index.stat().st_mtime_ns + 1))
        try:
            validate_annotation_config(config, "grch37")
        except ValueError as exc:
            assert "index is stale" in str(exc)
        else:
            raise AssertionError("expected stale index refusal")

        manifest = root / "manifest.json"
        manifest.write_text(json.dumps({"territory_bases": 2_000_000}))
        first, second, empty = root / "A.vcf.gz", root / "B.vcf.gz", root / "empty.vcf.gz"
        write_vcf(first, [
            "1\t10\t.\tA\tC,G\t.\tPASS\tCSQ=C|missense_variant|GENE1|G1|ENST1|1|YES,G|synonymous_variant|GENE1|G1|ENST2||",
            "1\t20\t.\tAT\tA\t.\tPASS\tCSQ=A|frameshift_variant|GENE2|G2|ENST3|1|YES",
            "1\t30\t.\tA\tT\t.\tLowQual\tCSQ=T|missense_variant|GENE3|G3|ENST4|1|YES",
        ])
        write_vcf(second, ["1\t20\t.\tAT\tA\t.\tPASS\tCSQ=A|frameshift_variant|GENE2|G2|ENST3|1|YES"])
        write_vcf(empty, [])
        burden = root / "burden.tsv"
        write_burden([("A", first), ("empty", empty)], manifest, burden)
        rows = list(csv.DictReader(burden.open(), delimiter="\t"))
        assert rows[0]["pass_snv_count"] == "2" and rows[0]["pass_indel_count"] == "1"
        assert rows[0]["total_per_mb"] == "1.5"
        assert rows[1]["pass_total_count"] == "0"
        recurrence = root / "recurrence.tsv"
        write_recurrence([("A", first), ("B", second)], recurrence, 2)
        recurrent = list(csv.DictReader(recurrence.open(), delimiter="\t"))
        assert len(recurrent) == 1
        assert recurrent[0]["variant_key"] == "1:20:AT:A"
        assert recurrent[0]["samples"] == "A,B" and recurrent[0]["pick_transcript"] == "ENST3"


if __name__ == "__main__":
    main()
