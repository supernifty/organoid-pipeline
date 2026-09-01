#!/usr/bin/env python3
"""Tests for WES/WGS mode safety and canonical territory preparation."""

from __future__ import annotations

import json
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from analysis_mode import (  # noqa: E402
    PRIMARY_CONTIG_LENGTHS,
    analysis_settings,
    prepare_wgs_territory,
    reference_settings,
    shard_count,
    validate_contamination,
    validate_existing_manifest,
    validate_reference_profile,
)


def expect_error(function, text: str) -> None:
    try:
        function()
    except ValueError as error:
        assert text in str(error), str(error)
    else:
        raise AssertionError(f"expected ValueError containing {text!r}")


def main() -> None:
    assert analysis_settings({})["type"] == "wes"
    expect_error(lambda: analysis_settings({"analysis": {"type": "rna"}}), "wes, wgs")
    assert reference_settings({})["build"] == "grch37"
    expect_error(lambda: reference_settings({"reference": {"build": "t2t"}}), "grch37, grch38")
    expect_error(
        lambda: reference_settings({"reference": {
            "build": "grch38", "regions_source_build": "grch37"
        }}),
        "remap the capture BED",
    )
    assert shard_count(40_000_001, 20_000_000, None) == 3
    assert shard_count(40_000_001, 20_000_000, 7) == 7

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        fai = root / "genome.fa.fai"
        source = root / "calling.interval_list"
        exclude = root / "exclude.bed"
        bed = root / "callable.bed"
        intervals = root / "callable.interval_list"
        metadata = root / "callable.json"
        fai.write_text("1\t1000\t0\t0\t0\n2\t500\t0\t0\t0\nMT\t100\t0\t0\t0\n")
        source.write_text(
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:500\n@SQ\tSN:MT\tLN:100\n"
            "1\t1\t400\t+\tCALLABLE\n1\t501\t1000\t+\tCALLABLE\n"
            "2\t1\t500\t+\tCALLABLE\nMT\t1\t100\t+\tCALLABLE\n"
        )
        exclude.write_text("1\t100\t200\n2\t0\t500\n")
        prepare_wgs_territory(source, fai, ["1", "2"], exclude, bed, intervals, metadata)
        assert bed.read_text().splitlines() == ["1\t0\t100", "1\t200\t400", "1\t500\t1000"]
        details = json.loads(metadata.read_text())
        assert details == {"callable_bases": 800, "contigs": ["1", "2"], "interval_count": 3}
        assert not any(
            line.startswith("MT\t")
            for line in intervals.read_text().splitlines()
            if not line.startswith("@")
        )

        results = root / "results"
        results.mkdir()
        (results / "legacy.txt").write_text("legacy")
        validate_existing_manifest(results, "wes")
        validate_existing_manifest(results, "wgs")
        (results / "analysis_manifest.json").write_text(json.dumps({"analysis_type": "wgs"}))
        validate_existing_manifest(results, "wgs")
        expect_error(lambda: validate_existing_manifest(results, "wes"), "contains 'wgs' outputs")
        (results / "analysis_manifest.json").write_text(json.dumps({
            "analysis_type": "wgs", "reference_build": "grch38"
        }))
        validate_existing_manifest(results, "wgs", "grch38")
        expect_error(
            lambda: validate_existing_manifest(results, "wgs", "grch37"),
            "reference.build",
        )

        pileup = root / "pileup.table"
        contamination = root / "contamination.table"
        pileup.write_text("contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n")
        contamination.write_text("sample\tNaN\n")
        expect_error(lambda: validate_contamination(pileup, contamination), "No usable")
        pileup.write_text(pileup.read_text() + "1\t10\t10\t2\t0\t0.1\n")
        expect_error(lambda: validate_contamination(pileup, contamination), "finite")
        contamination.write_text("sample\t0.012\n")
        validate_contamination(pileup, contamination)

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        lengths = PRIMARY_CONTIG_LENGTHS["grch38"]
        fai = root / "genome.fa.fai"
        dictionary = root / "genome.dict"
        territory = root / "capture.bed"
        resource = root / "resource.vcf"
        report = root / "profile.json"
        fai.write_text("".join(f"{name}\t{length}\t0\t0\t0\n" for name, length in lengths.items()))
        dictionary.write_text(
            "@HD\tVN:1.6\n" + "".join(
                f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in lengths.items()
            )
        )
        territory.write_text("chr1\t100\t200\n")
        resource.write_text(
            "##fileformat=VCFv4.2\n"
            f"##contig=<ID=chr1,length={lengths['chr1']}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        validate_reference_profile(
            "grch38", fai, dictionary, territory, list(lengths), [resource], report
        )
        assert json.loads(report.read_text())["reference_build"] == "grch38"

    rules = (ROOT / "workflow" / "rules" / "variant_calling.smk").read_text()
    assert "--matched-normal" in rules
    assert "--interval-padding 0" in rules
    assert "strelka_mode_params()" in rules
    assert "--restrict-alleles-to BIALLELIC" in rules
    assert "MateOnSameContigOrNoMappedMateReadFilter" in rules
    pon_shard_rule = rules.split("rule mutect2_pon_shard:", 1)[1].split(
        "rule mutect2_panel_of_normals:", 1
    )[0]
    assert 'get_container_cmd(config["gatk"], bind_node_tmp=True)' in pon_shard_rule
    assert 'test -s "$workspace/callset.json"' in pon_shard_rule
    assert 'test -s "$workspace/vidmap.json"' in pon_shard_rule
    assert "--overwrite-existing-genomicsdb-workspace true" not in pon_shard_rule

    with tempfile.TemporaryDirectory() as directory:
        coverage = Path(directory) / "coverage_mqc.tsv"
        regions = Path(directory) / "sample.regions.bed.gz"
        with gzip.open(regions, "wt") as handle:
            handle.write("1\t0\t2\tregion1\t15\nX\t0\t1\tregionX\t6\n")
        subprocess.run(
            [sys.executable, str(ROOT / "workflow" / "scripts" / "wgs_coverage_mqc.py"),
             "--sample", "T1", "--regions", str(regions), "--output", str(coverage)],
            check=True,
        )
        assert coverage.read_text().splitlines()[-2:] == [
            "sample\t1\tX", "T1\t15.0\t6.0"
        ]

    # Parse and build the WGS preparation DAG against a tiny synthetic GRCh37-like fixture.
    if not (ROOT / "results").exists():
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "genome.fa"
            reference.write_text(">1\nA\n>2\nA\n")
            Path(f"{reference}.fai").write_text("1\t1000\t0\t0\t0\n2\t500\t0\t0\t0\n")
            reference_dict = root / "genome.dict"
            reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:500\n")
            source = root / "wgs.interval_list"
            source.write_text(
                "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:2\tLN:500\n"
                "1\t1\t1000\t+\tCALLABLE\n2\t1\t500\t+\tCALLABLE\n"
            )
            config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
            config["analysis"]["type"] = "wgs"
            config["analysis"]["wgs"]["contigs"] = ["1", "2"]
            config["analysis"]["wgs"]["scatter_count"] = 2
            config["reference"]["genome"] = str(reference)
            config["reference"]["genome_dict"] = str(reference_dict)
            config["reference"]["wgs_calling_regions"] = str(source)
            config["reference"]["wgs_exclude_regions"] = None
            config["storage"]["tmp_dir"] = str(root / "tmp")
            config["storage"]["local_scratch"] = str(root / "scratch")
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            subprocess.run(
                [
                    str(Path(sys.executable).with_name("snakemake")),
                    "-n",
                    "--cores",
                    "1",
                    "--configfile",
                    str(config_path),
                    str(root / "tmp" / "analysis" / "mutect2_shards"),
                ],
                cwd=ROOT,
                env={**os.environ, "XDG_CACHE_HOME": str(root / "cache")},
                check=True,
                stdout=subprocess.DEVNULL,
            )

    # Build the GRCh38 WES reference-profile DAG against a dictionary-compatible fixture.
    if not (ROOT / "results").exists():
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lengths = PRIMARY_CONTIG_LENGTHS["grch38"]
            reference = root / "genome.fa"
            reference.write_text(">chr1\nA\n")
            Path(f"{reference}.fai").write_text(
                "".join(f"{name}\t{length}\t0\t0\t0\n" for name, length in lengths.items())
            )
            for extension in ("amb", "ann", "bwt", "pac", "sa", "alt"):
                Path(f"{reference}.64.{extension}").write_text("fixture\n")
            reference_dict = root / "genome.dict"
            reference_dict.write_text(
                "@HD\tVN:1.6\n" + "".join(
                    f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in lengths.items()
                )
            )
            bed = root / "capture.bed"
            bed.write_text("chr1\t100\t200\n")
            bed_gz = root / "capture.bed.gz"
            vcf = root / "gnomad.vcf"
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                f"##contig=<ID=chr1,length={lengths['chr1']}>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )
            vcf_gz = root / "gnomad.vcf.gz"
            bgzip = Path(sys.executable).with_name("bgzip")
            tabix = Path(sys.executable).with_name("tabix")
            with bed_gz.open("wb") as handle:
                subprocess.run([str(bgzip), "-c", str(bed)], stdout=handle, check=True)
            subprocess.run([str(tabix), "-f", "-p", "bed", str(bed_gz)], check=True)
            with vcf_gz.open("wb") as handle:
                subprocess.run([str(bgzip), "-c", str(vcf)], stdout=handle, check=True)
            subprocess.run([str(tabix), "-f", "-p", "vcf", str(vcf_gz)], check=True)

            config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
            config["reference"].update({
                "build": "grch38",
                "genome": str(reference),
                "genome_dict": str(reference_dict),
                "regions": str(bed_gz),
                "regions_source": "vendor",
                "regions_source_build": "grch38",
                "gnomad": str(vcf_gz),
                "contamination_sites": None,
            })
            config["chromosomes"] = list(lengths)
            config["mutational_signatures"]["reference_build"] = "grch38"
            config["hotspots"]["reference_build"] = "grch38"
            config["storage"]["tmp_dir"] = str(root / "tmp")
            config["storage"]["local_scratch"] = str(root / "scratch")
            config_path = root / "config.yaml"
            config_path.write_text(yaml.safe_dump(config, sort_keys=False))
            subprocess.run(
                [
                    str(Path(sys.executable).with_name("snakemake")),
                    "-n", "--cores", "1", "--configfile", str(config_path),
                    str(root / "tmp" / "analysis" / "reference_profile.json"),
                ],
                cwd=ROOT,
                env={**os.environ, "XDG_CACHE_HOME": str(root / "cache")},
                check=True,
                stdout=subprocess.DEVNULL,
            )


if __name__ == "__main__":
    main()
