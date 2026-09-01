#!/usr/bin/env python3
"""Validate matched-normal germline configuration and DAG construction."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import yaml
from annotation import checksum_path


TASK_ROOT = ROOT / "tmp" / "codex" / "test-germline"
SNAKEMAKE = Path(sys.executable).with_name("snakemake")
BGZIP = Path(sys.executable).with_name("bgzip")
TABIX = Path(sys.executable).with_name("tabix")


def bgzip_and_index(path: Path, preset: str) -> Path:
    compressed = Path(f"{path}.gz")
    with compressed.open("wb") as handle:
        subprocess.run([str(BGZIP), "-c", str(path)], stdout=handle, check=True)
    subprocess.run([str(TABIX), "-f", "-p", preset, str(compressed)], check=True)
    return compressed


def fixture(name: str, mode: str, input_mode: str) -> tuple[Path, Path]:
    root = TASK_ROOT / name
    root.mkdir(parents=True)
    reference = root / "genome.fa"
    reference.write_text(">1\nA\n")
    Path(f"{reference}.fai").write_text("1\t1000\t0\t0\t0\n")
    (root / "genome.dict").write_text("@HD\tVN:1.6\n@SQ\tSN:1\tLN:1000\n")
    for extension in ("amb", "ann", "bwt", "pac", "sa"):
        Path(f"{reference}.{extension}").write_text("fixture\n")

    bed = root / "regions.bed"
    bed.write_text("1\t0\t1000\n")
    regions = bgzip_and_index(bed, "bed")
    gnomad_vcf = root / "gnomad.vcf"
    gnomad_vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    gnomad = bgzip_and_index(gnomad_vcf, "vcf")
    intervals = root / "wgs.interval_list"
    intervals.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n"
        "1\t1\t1000\t+\tCALLABLE\n"
    )

    if input_mode == "cram":
        for sample in ("T1", "T2", "N"):
            (root / f"{sample}.cram").write_bytes(b"fixture")
            (root / f"{sample}.cram.crai").write_bytes(b"fixture")
        samples = {
            "samples": {
                sample: {"cram": str(root / f"{sample}.cram")}
                for sample in ("T1", "T2", "N")
            },
            "tumours": {"T1": "N", "T2": "N"},
        }
    elif input_mode == "fastq":
        for sample in ("T", "N"):
            for read in (1, 2):
                (root / f"{sample}_R{read}.fastq.gz").write_bytes(b"fixture")
        samples = {
            "samples": {
                sample: {
                    "fastq_1": str(root / f"{sample}_R1.fastq.gz"),
                    "fastq_2": str(root / f"{sample}_R2.fastq.gz"),
                }
                for sample in ("T", "N")
            },
            "tumours": {"T": "N"},
        }
    else:
        final_vcf = root / "T.final.vcf.gz"
        final_vcf.write_bytes(b"fixture")
        samples = {"samples": {"T": {"final_vcf": str(final_vcf)}}}

    samples_path = root / "samples.yaml"
    samples_path.write_text(yaml.safe_dump(samples, sort_keys=False))
    config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    config["analysis"]["type"] = mode
    config["analysis"]["wgs"]["contigs"] = ["1"]
    config["analysis"]["wgs"]["scatter_count"] = 1
    config["reference"].update({
        "genome": str(reference),
        "genome_dict": str(root / "genome.dict"),
        "regions": str(regions),
        "wgs_calling_regions": str(intervals),
        "wgs_exclude_regions": None,
        "gnomad": str(gnomad),
        "contamination_sites": None,
        "panel_of_normals": str(gnomad),
    })
    config["chromosomes"] = ["1"]
    config["storage"]["tmp_dir"] = str(root / "work")
    config["storage"]["local_scratch"] = str(root / "scratch")
    config_path = root / "config.yaml"
    config["run_management"] = {"samples_file": str(samples_path), "config_file": str(config_path)}
    if input_mode == "vcf":
        resources = {}
        for name in ("cache", "plugins", "clinvar", "gnomad", "revel", "cadd_snv", "cadd_indel", "alphamissense", "spliceai_snv", "spliceai_indel"):
            path = root / f"annotation-{name}"
            if name in ("cache", "plugins"):
                path.mkdir()
                (path / "release.txt").write_text("116\n")
            else:
                path.write_text(name)
            resource = {"build": "grch37", "version": "116", "path": str(path),
                        "expected_sha256": checksum_path(path), "source": "fixture", "licence": "fixture", "access_date": "2026-08-25"}
            if path.is_file():
                index = Path(f"{path}.tbi")
                index.write_text("index")
                resource["index"] = {"path": str(index), "expected_sha256": checksum_path(index)}
            resources[name] = resource
        config["annotation"].update(enabled=True, resources=resources)
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    return root, config_path


def dry_run(
    root: Path, config: Path, target: str | None = None, slurm_profile: bool = False
) -> str:
    (root / "python-tmp").mkdir(exist_ok=True)
    (root / "cache").mkdir(exist_ok=True)
    command = [str(SNAKEMAKE)]
    if target:
        command.append(target)
    command.extend(["-n", "-p", "--cores", "64" if slurm_profile else "1", "--configfile", str(config)])
    if slurm_profile:
        command.extend(["--profile", str(ROOT / "config" / "slurm")])
    completed = subprocess.run(
        command,
        cwd=ROOT,
        env={
            **os.environ,
            "TMPDIR": str(root / "python-tmp"),
            "XDG_CACHE_HOME": str(root / "cache"),
        },
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode:
        raise AssertionError(completed.stdout)
    return completed.stdout


def main() -> None:
    if TASK_ROOT.exists():
        shutil.rmtree(TASK_ROOT)
    TASK_ROOT.mkdir(parents=True)
    try:
        rules = (ROOT / "workflow" / "rules" / "germline.smk").read_text()
        config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
        assert config["germline"]["enabled"] is True
        assert config["germline"]["max_concurrent_haplotypecaller_shards"] > 0
        assert config["germline"]["hard_filters"] == {
            "snp": {
                "qd_min": 2.0, "qual_min": 30.0, "sor_max": 3.0,
                "fs_max": 60.0, "mq_min": 40.0, "mq_rank_sum_min": -12.5,
                "read_pos_rank_sum_min": -8.0,
            },
            "indel": {
                "qd_min": 2.0, "qual_min": 30.0, "fs_max": 200.0,
                "read_pos_rank_sum_min": -20.0,
            },
        }
        for filter_name in (
            "SNP_QD", "SNP_QUAL", "SNP_SOR", "SNP_FS", "SNP_MQ",
            "SNP_MQRankSum", "SNP_ReadPosRankSum", "INDEL_QD", "INDEL_QUAL",
            "INDEL_FS", "INDEL_ReadPosRankSum",
        ):
            assert filter_name in rules
        assert "--select-type-to-include MIXED" in rules
        assert "--exclude-non-variants true" in rules

        wes_cram, wes_cram_config = fixture("wes-cram", "wes", "cram")
        output = dry_run(
            wes_cram,
            wes_cram_config,
            "results/germline/N.haplotypecaller.filtered.vcf.gz",
        )
        for rule in (
            "germline_haplotypecaller_shard", "germline_merge_gvcf",
            "germline_genotype", "germline_filter",
        ):
            assert rule in output, output
        for expression in (
            "QD < 2.0", "QUAL < 30.0", "SOR > 3.0", "FS > 60.0",
            "MQ < 40.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0",
            "FS > 200.0", "ReadPosRankSum < -20.0",
        ):
            assert expression in output
        assert output.count("normal=N") >= 1

        wes_fastq, wes_fastq_config = fixture("wes-fastq", "wes", "fastq")
        output = dry_run(
            wes_fastq,
            wes_fastq_config,
            "results/germline/N.haplotypecaller.g.vcf.gz",
        )
        assert "mark_duplicates" in output
        assert "germline_haplotypecaller_shard" in output
        assert "sort_alignment_cram" not in output

        alignment_rules = (ROOT / "workflow" / "rules" / "alignment.smk").read_text()
        assert "rule sort_alignment_cram:" not in alignment_rules
        paired_rule = alignment_rules.split("rule bwa_mem_paired:", 1)[1].split("rule mark_duplicates:", 1)[0]
        assert "threads: 32" in paired_rule and "bwa mem -M -t 24" in paired_rule and "samtools sort -@ 8" in paired_rule
        mark_rule = alignment_rules.split("rule mark_duplicates:", 1)[1]
        assert 'gatk_java_options("44g")' in mark_rule and "threads: 1" in mark_rule
        slurm = yaml.safe_load((ROOT / "config" / "slurm" / "config.yaml").read_text())
        assert slurm["set-resources"]["mark_duplicates"]["mem_mb"] == 49152
        output = dry_run(wes_fastq, wes_fastq_config, "results/cram/T.sorted.dups.cram", slurm_profile=True)
        assert "rule bwa_mem_paired:" in output and "threads: 32" in output
        mark_job = output.split("rule mark_duplicates:", 1)[1]
        assert "mem_mb=49152" in mark_job

        wgs_cram, wgs_cram_config = fixture("wgs-cram", "wgs", "cram")
        output = dry_run(
            wgs_cram,
            wgs_cram_config,
            "results/germline/N.haplotypecaller.filtered.vcf.gz",
        )
        assert "prepare_wgs_territory" in output and "split_wgs_intervals" in output
        assert "germline_haplotypecaller_shard" in output
        assert "germline_filter" in output

        output = dry_run(
            wgs_cram,
            wgs_cram_config,
            "results/germline/N.haplotypecaller.filtered.vcf.gz",
            slurm_profile=True,
        )
        assert "haplotypecaller_shards=1" in output

        vcf_only, vcf_only_config = fixture("vcf-only", "wes", "vcf")
        output = dry_run(vcf_only, vcf_only_config)
        assert "annotate_somatic_vep" in output
        assert "aggregate_mutation_burden" in output
        assert "aggregate_recurrent_variants" in output
        assert "annotation_resource_report" in output
        assert "germline_haplotypecaller_shard" not in output
        assert "germline_genotype" not in output
        assert "germline_filter" not in output
        assert "alignment_summary" not in output
        assert "wgs_per_contig_coverage" not in output
        assert "exon_gene_coverage" not in output
        assert "somalier_extract" not in output

        wgs_vcf_only, wgs_vcf_only_config = fixture("wgs-vcf-only", "wgs", "vcf")
        output = dry_run(wgs_vcf_only, wgs_vcf_only_config)
        assert "annotate_somatic_vep" in output and "aggregate_mutation_burden" in output
        assert "prepare_wgs_territory" in output and "split_wgs_intervals" not in output
        assert "alignment_summary" not in output and "wgs_per_contig_coverage" not in output
    finally:
        shutil.rmtree(TASK_ROOT)


if __name__ == "__main__":
    main()
