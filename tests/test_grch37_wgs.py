"""Tests for native GRCh37 WGS provisioning and reference safety."""

from __future__ import annotations

import gzip
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import downsample_alignment  # noqa: E402
import provision_grch37  # noqa: E402
from analysis_mode import PRIMARY_CONTIG_LENGTHS  # noqa: E402


def test_grch37_overlay_selects_native_build(tmp_path):
    overlay = provision_grch37.config_text(tmp_path)
    assert "build: grch37" in overlay
    assert 'contigs: ["1", "2", "3"' in overlay
    assert "Homo_sapiens_assembly19.fasta" in overlay
    assert "af-only-gnomad.grch37.vcf.gz" in overlay
    assert "regions_source_build: grch37" in overlay
    assert "reference_build: grch37" in overlay


def test_grch37_interval_list_uses_primary_contigs(tmp_path):
    lengths = PRIMARY_CONTIG_LENGTHS["grch37"]
    dictionary = tmp_path / "genome.dict"
    dictionary.write_text(
        "@HD\tVN:1.6\n"
        + "".join(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in lengths.items())
        + "@SQ\tSN:GL000207.1\tLN:4262\n"
    )
    output = tmp_path / "wgs.interval_list"
    provision_grch37.create_interval_list(dictionary, output)
    records = [line for line in output.read_text().splitlines() if not line.startswith("@")]
    assert len(records) == 24
    assert records[0] == "1\t1\t249250621\t+\tCALLABLE"
    assert records[-1] == "Y\t1\t59373566\t+\tCALLABLE"
    assert not any(line.startswith("GL") for line in records)


def test_vcf_conversion_adds_reference_dictionary_and_index(tmp_path):
    source = tmp_path / "source.vcf"
    source.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t2\t.\tA\tC\t.\tPASS\tAF=0.1\n"
    )
    fai = tmp_path / "genome.fa.fai"
    fai.write_text("1\t100\t0\t0\t0\n")
    output = tmp_path / "derived.vcf.gz"
    bgzip = shutil.which("bgzip")
    tabix = shutil.which("tabix")
    assert bgzip and tabix
    provision_grch37.bgzip_with_reference_dictionary(
        source,
        output,
        fai,
        threads=1,
        bgzip=bgzip,
        tabix=tabix,
    )
    with gzip.open(output, "rt") as handle:
        text = handle.read()
    assert "##contig=<ID=1,length=100>" in text
    assert "1\t2\t.\tA\tC" in text
    assert Path(f"{output}.tbi").is_file()


def test_derived_resource_requires_manifest_checksums(tmp_path):
    output = tmp_path / "derived.vcf.gz"
    index = tmp_path / "derived.vcf.gz.tbi"
    output.write_bytes(b"vcf")
    index.write_bytes(b"index")
    manifest = {
        "generated": [
            {
                "filename": output.name,
                "size": output.stat().st_size,
                "sha256": provision_grch37.sha256(output),
            },
            {
                "filename": index.name,
                "size": index.stat().st_size,
                "sha256": provision_grch37.sha256(index),
            },
        ]
    }
    assert provision_grch37.derived_ready(tmp_path, output.name, manifest)
    output.write_bytes(b"changed")
    assert not provision_grch37.derived_ready(tmp_path, output.name, manifest)


def test_alignment_reference_mismatch_is_rejected():
    grch37 = {"1": 249250621, "GL000207.1": 4262}
    grch38 = {"chr1": 248956422}
    with pytest.raises(ValueError, match="required contig"):
        downsample_alignment.validate_dictionaries(grch37, grch38, ["chr1"], {"1": 10})
    with pytest.raises(ValueError, match="mapped reads use contigs"):
        downsample_alignment.validate_dictionaries(
            grch37, {"1": 249250621}, ["1"], {"1": 10, "GL000207.1": 5}
        )
    downsample_alignment.validate_dictionaries(
        grch37, {"1": 249250621}, ["1"], {"1": 10, "GL000207.1": 0}
    )


def test_idxstats_command_uses_reference_only_for_cram():
    assert downsample_alignment.idxstats_command("input.bam", "/refs/genome.fa") == [
        "samtools",
        "idxstats",
        "input.bam",
    ]
    assert downsample_alignment.idxstats_command("input.cram", "/refs/genome.fa") == [
        "samtools",
        "idxstats",
        "--input-fmt-option",
        "reference=/refs/genome.fa",
        "input.cram",
    ]


def test_checked_output_surfaces_external_stderr(monkeypatch):
    class Failed:
        returncode = 1
        stdout = ""
        stderr = "index is stale"

    monkeypatch.setattr(downsample_alignment.subprocess, "run", lambda *args, **kwargs: Failed())
    with pytest.raises(RuntimeError, match="alignment-index.*index is stale"):
        downsample_alignment.checked_output(
            ["samtools", "idxstats", "input.bam"], "alignment-index"
        )


def test_grch37_manifest_is_machine_readable(tmp_path):
    payload = {
        "bundle": "Broad GRCh37/hg19 v0 plus GATK somatic-b37",
        "resources": [],
        "generated": [],
    }
    path = tmp_path / "resource-manifest.json"
    path.write_text(json.dumps(payload))
    assert provision_grch37.load_manifest(tmp_path)["bundle"].startswith("Broad GRCh37")
