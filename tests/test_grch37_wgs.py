"""Tests for native GRCh37 WGS provisioning and reference safety."""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import downsample_alignment  # noqa: E402
import provision_grch37  # noqa: E402
import validate_alignment  # noqa: E402
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


def test_downsampling_forces_configured_cram_version():
    command = downsample_alignment.cram_command(
        "input.bam", "genome.fa", "output.cram", 8, 0.25, 1723, "3.0"
    )
    assert command == [
        "samtools",
        "view",
        "-@",
        "8",
        "-T",
        "genome.fa",
        "-C",
        "--output-fmt-option",
        "version=3.0",
        "-s",
        "1723.25",
        "-o",
        "output.cram",
        "input.bam",
    ]


def test_downsampling_command_writes_cram_3_0(tmp_path):
    reference = tmp_path / "genome.fa"
    reference.write_text(">1\n" + "A" * 100 + "\n>2\n" + "C" * 80 + "\n")
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    reference_dict = tmp_path / "genome.dict"
    subprocess.run(["samtools", "dict", "-o", str(reference_dict), str(reference)], check=True)
    sam = tmp_path / "input.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:1\tLN:100\n"
        "@RG\tID:rg\tSM:sample\n"
        "read1\t0\t1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tFFFFFFFFFF\tRG:Z:rg\n"
    )
    bam = tmp_path / "input.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    raw = tmp_path / "raw.cram"
    command = downsample_alignment.cram_command(bam, reference, raw, 1, 1.0, 1723, "3.0")
    subprocess.run(command, check=True)
    cram = tmp_path / "output.cram"
    metadata = downsample_alignment.expand_cram_dictionary(
        raw,
        cram,
        reference,
        reference_dict,
        tmp_path / "expanded.sam",
    )
    subprocess.run(["samtools", "index", str(cram)], check=True)
    assert validate_alignment.cram_version(cram) == "3.0"
    header = subprocess.run(
        ["samtools", "view", "-H", "-T", str(reference), str(cram)],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    assert list(validate_alignment.parse_header(header)[1]) == ["1", "2"]
    assert metadata["source_sequence_count"] == 1
    assert metadata["output_sequence_count"] == 2
    raw_records = subprocess.run(
        ["samtools", "view", "-T", str(reference), str(raw)],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    expanded_records = subprocess.run(
        ["samtools", "view", "-T", str(reference), str(cram)],
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    assert expanded_records == raw_records
    payload = validate_alignment.validate(
        str(cram),
        f"{cram}.crai",
        str(reference),
        f"{reference}.fai",
        "sample",
        ["1"],
        "3.0",
    )
    assert payload["cram_version"] == "3.0"


def test_header_expansion_rejects_reference_external_source_contig(tmp_path):
    reference_dict = tmp_path / "genome.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:1\tLN:100\n")
    source = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:other\tLN:100\n"
    with pytest.raises(ValueError, match="required contig other"):
        downsample_alignment.expanded_header(source, reference_dict)


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


def test_production_cram_idxstats_uses_supported_reference_option():
    assert validate_alignment.idxstats_command("input.cram", "/refs/genome.fa") == [
        "samtools",
        "idxstats",
        "--input-fmt-option",
        "reference=/refs/genome.fa",
        "input.cram",
    ]


def test_production_preflight_reads_cram_version(tmp_path):
    cram = tmp_path / "input.cram"
    cram.write_bytes(b"CRAM\x03\x00fixture")
    assert validate_alignment.cram_version(cram) == "3.0"
    cram.write_bytes(b"CRAM\x03\x01fixture")
    assert validate_alignment.cram_version(cram) == "3.1"
    cram.write_bytes(b"not-cram")
    with pytest.raises(ValueError, match="valid CRAM file header"):
        validate_alignment.cram_version(cram)


def test_production_preflight_surfaces_external_stderr(monkeypatch):
    class Failed:
        returncode = 1
        stdout = ""
        stderr = "reference could not be loaded"

    monkeypatch.setattr(validate_alignment.subprocess, "run", lambda *args, **kwargs: Failed())
    with pytest.raises(ValueError, match="CRAM index.*reference could not be loaded"):
        validate_alignment.checked_output(["samtools", "idxstats", "input.cram"], "CRAM index")


def test_grch37_manifest_is_machine_readable(tmp_path):
    payload = {
        "bundle": "Broad GRCh37/hg19 v0 plus GATK somatic-b37",
        "resources": [],
        "generated": [],
    }
    path = tmp_path / "resource-manifest.json"
    path.write_text(json.dumps(payload))
    assert provision_grch37.load_manifest(tmp_path)["bundle"].startswith("Broad GRCh37")
