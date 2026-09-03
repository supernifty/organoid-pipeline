import csv
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from allele_recount import count_observations  # noqa: E402
from benchmark_metrics import cosine_similarity, metrics  # noqa: E402
from caller_tiers import build_tiers  # noqa: E402
from sample_inputs import comparison_map, matched_normal_samples, validate_samples  # noqa: E402
from sbs96 import canonical_channel  # noqa: E402


def test_workflow_helper_scripts_are_repository_root_safe():
    workflow_sources = [ROOT / "Snakefile", *sorted((ROOT / "workflow/rules").glob("*.smk"))]
    unsafe = []
    pattern = re.compile(r"\bpython3?\s+scripts/[A-Za-z0-9_.-]+\.py")
    for path in workflow_sources:
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            if pattern.search(line):
                unsafe.append(f"{path.relative_to(ROOT)}:{line_number}: {line.strip()}")
    assert unsafe == []


def manifest():
    return {
        "samples": {
            "B": {"role": "baseline", "donor": "D", "lineage": "L", "cram": "B.cram"},
            "O1": {"role": "organoid", "donor": "D", "lineage": "L", "cram": "O1.cram"},
            "O2": {
                "role": "organoid",
                "donor": "D",
                "lineage": "L",
                "fastq_1": "1.fq.gz",
                "fastq_2": "2.fq.gz",
                "read_group": {"platform": "ILLUMINA", "library": "L", "unit": "U"},
            },
        },
        "comparisons": {"O1": {"baseline": "B"}, "O2": {"baseline": "B"}},
    }


def test_manifest_supports_shared_baseline():
    value = manifest()
    validate_samples(value, "fixture")
    assert comparison_map(value, "fixture") == {"O1": "B", "O2": "B"}
    assert matched_normal_samples(value, "fixture") == ["B"]


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda value: value["samples"]["O1"].update(donor="other"), "inconsistent donor"),
        (lambda value: value["samples"]["O1"].update(fastq_1="x"), "both 'fastq_1' and 'fastq_2'"),
        (lambda value: value["comparisons"].pop("O1"), "organoids without comparisons"),
        (lambda value: value["samples"]["O2"].pop("read_group"), "read_group metadata"),
        (
            lambda value: value["samples"]["O1"].update(crai="indexes/O1.crai"),
            "samtools-discoverable path",
        ),
    ],
)
def test_manifest_rejects_unsafe_models(change, message):
    value = manifest()
    change(value)
    with pytest.raises(ValueError, match=message):
        validate_samples(value, "fixture")


def write_vcf(path, rows):
    path.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "".join(rows)
    )


def test_exact_allele_caller_tiers(tmp_path):
    mutect2, strelka = tmp_path / "m.vcf", tmp_path / "s.vcf"
    write_vcf(
        mutect2,
        [
            "chr1\t2\t.\tA\tC\t.\tPASS\tM2=shared\n",
            "chr1\t3\t.\tG\tT\t.\tPASS\tM2=only\n",
        ],
    )
    write_vcf(strelka, ["chr1\t2\t.\tA\tC\t.\tPASS\t.\n", "chr1\t4\t.\tC\tT\t.\tPASS\t.\n"])
    _, m2, st, both, union, support = build_tiers(str(mutect2), [str(strelka)])
    assert len(m2) == 2 and len(st) == 2 and len(both) == 1 and len(union) == 3
    assert support[("chr1", 2, "A", "C")] == "Mutect2,Strelka2"
    assert support[("chr1", 3, "G", "T")] == "Mutect2"
    assert support[("chr1", 4, "C", "T")] == "Strelka2"
    assert union[("chr1", 2, "A", "C")][7] == "M2=shared"


def test_recount_snv_and_indel_strands():
    assert count_observations("A", "C", ".,Cc")[:4] == (1, 1, 1, 1)
    assert count_observations("A", "AT", ".+1T,+1t")[:4] == (0, 0, 1, 1)
    assert count_observations("AT", "A", ".-1T,-1t")[:4] == (0, 0, 1, 1)


def test_sbs96_canonicalization():
    assert canonical_channel("ACA", "C", "T") == "A[C>T]A"
    assert canonical_channel("TGT", "G", "A") == "A[C>T]A"


def test_benchmark_metrics():
    result = metrics({1, 2, 3}, {2, 3, 4}, 1_000_000_000)
    assert result["precision"] == pytest.approx(2 / 3)
    assert result["recall"] == pytest.approx(2 / 3)
    assert result["f1"] == pytest.approx(2 / 3)
    assert result["false_positives_per_callable_gb"] == 1
    assert cosine_similarity({"a": 1, "b": 1}, {"a": 2, "b": 2}) == pytest.approx(1)


def test_high_vaf_benchmark_does_not_count_lower_vaf_truth_as_false_positive(tmp_path):
    truth, calls = tmp_path / "truth.vcf", tmp_path / "calls.vcf"
    write_vcf(
        truth,
        [
            "chr1\t2\t.\tA\tC\t.\tPASS\tVAF=0.30\n",
            "chr1\t3\t.\tA\tG\t.\tPASS\tVAF=0.10\n",
        ],
    )
    write_vcf(
        calls,
        [
            "chr1\t2\t.\tA\tC\t.\tPASS\t.\n",
            "chr1\t3\t.\tA\tG\t.\tPASS\t.\n",
            "chr1\t4\t.\tA\tT\t.\tPASS\t.\n",
        ],
    )
    output = tmp_path / "metrics.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow/scripts/benchmark_metrics.py"),
            "--truth",
            str(truth),
            "--calls",
            str(calls),
            "--callable-bases",
            "100",
            "--output",
            str(output),
        ],
        check=True,
    )
    high_vaf = yaml.safe_load(output.read_text())["truth_vaf_ge_0.25"]
    assert high_vaf["tp"] == 1 and high_vaf["fp"] == 1


def test_population_recurrence_and_evidence_filters(tmp_path):
    sample_manifest = manifest()
    manifest_path = tmp_path / "samples.yaml"
    manifest_path.write_text(yaml.safe_dump(sample_manifest))
    o1, o2, population = tmp_path / "O1.vcf", tmp_path / "O2.vcf", tmp_path / "pop.vcf"
    variants = [
        "chr1\t2\t.\tA\tC\t.\tPASS\tCALLER_SUPPORT=Mutect2\n",
        "chr1\t3\t.\tA\tG\t.\tPASS\tCALLER_SUPPORT=Strelka2\n",
        "chr1\t4\t.\tA\tT\t.\tPASS\tCALLER_SUPPORT=Mutect2,Strelka2\n",
        "chr1\t5\t.\tC\tA\t.\tPASS\tCALLER_SUPPORT=Mutect2\n",
        "chr1\t6\t.\tC\tG\t.\tPASS\tCALLER_SUPPORT=Mutect2\n",
        "chr1\t7\t.\tC\tT\t.\tPASS\tCALLER_SUPPORT=Mutect2\n",
    ]
    write_vcf(o1, variants)
    write_vcf(o2, [variants[1]])
    write_vcf(population, ["chr1\t2\t.\tA\tC\t.\tPASS\tAF=0.0001\n"])
    counts_path = tmp_path / "counts.tsv"
    fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "sample",
        "depth",
        "ref_count",
        "alt_count",
        "vaf",
        "ref_forward",
        "ref_reverse",
        "alt_forward",
        "alt_reverse",
        "mean_baseq",
        "mean_mapq",
    ]
    with counts_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in variants:
            columns = line.split("\t")
            pos, ref, alt = int(columns[1]), columns[3], columns[4]
            later_alt = 1 if pos == 5 else 3
            baseline_depth = 2 if pos == 6 else 10
            baseline_alt = 2 if pos == 7 else 0
            for sample, depth, alt_count in (
                ("O1", 10, later_alt),
                ("B", baseline_depth, baseline_alt),
            ):
                writer.writerow(
                    {
                        "chrom": "chr1",
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "sample": sample,
                        "depth": depth,
                        "ref_count": depth - alt_count,
                        "alt_count": alt_count,
                        "vaf": alt_count / depth,
                        "ref_forward": 0,
                        "ref_reverse": 0,
                        "alt_forward": 0,
                        "alt_reverse": 0,
                        "mean_baseq": 30,
                        "mean_mapq": 60,
                    }
                )
    prefix = tmp_path / "catalog"
    command = [
        sys.executable,
        str(ROOT / "workflow/scripts/catalog_filter.py"),
        "--sample",
        "O1",
        "--manifest",
        str(manifest_path),
        "--union",
        f"O1={o1}",
        "--union",
        f"O2={o2}",
        "--counts",
        str(counts_path),
        "--population-vcf",
        str(population),
    ]
    for option in (
        "audit-vcf",
        "audit-tsv",
        "stringent-vcf",
        "stringent-tsv",
        "sensitivity-vcf",
        "sensitivity-tsv",
        "rejected-vcf",
        "rejected-tsv",
        "shared-lineage-vcf",
        "shared-lineage-tsv",
        "stage-counts",
    ):
        command.extend((f"--{option}", str(prefix) + f".{option}"))
    subprocess.run(command, check=True)
    stringent = (Path(str(prefix) + ".stringent-tsv")).read_text()
    sensitivity = (Path(str(prefix) + ".sensitivity-tsv")).read_text()
    shared = (Path(str(prefix) + ".shared-lineage-tsv")).read_text()
    assert "\t4\tA\tT\t" in stringent
    assert "\t2\tA\tC\t" not in stringent and "\t2\tA\tC\t" in sensitivity
    assert "\t3\tA\tG\t" in shared and "SHARED_LINEAGE" in shared
    rejected = (Path(str(prefix) + ".rejected-tsv")).read_text()
    assert "LOW_LATER_ALT" in rejected and "LOW_BASELINE_DEPTH" in rejected
    assert "BASELINE_ALT_EVIDENCE" in rejected


def test_trimming_absent_from_active_dag_and_dependencies():
    snakefile = (ROOT / "Snakefile").read_text()
    pixi = (ROOT / "pixi.toml").read_text()
    config = (ROOT / "config" / "config.yaml").read_text()
    assert "trimming.smk" not in snakefile
    assert "trimmomatic" not in pixi.lower()
    assert "trimmomatic" not in config.lower()


def test_apptainer_runtime_is_supported():
    snakefile = (ROOT / "Snakefile").read_text()
    pull_images = (ROOT / "scripts/pull_images.sh").read_text()
    assert 'runtime in {"apptainer", "singularity"}' in snakefile
    assert "apptainer|singularity)" in pull_images


def test_slurm_walltimes_are_placeable_for_low_depth_wgs():
    profile = yaml.safe_load((ROOT / "config/slurm/config.yaml").read_text())["set-resources"]
    for rule in (
        "mutect2_sample_pon_chromosome",
        "mutect2_pon_shard",
        "mutect2_chromosome",
    ):
        assert profile[rule]["runtime"] == 720
    for rule in ("bwa_mem", "bwa_mem_paired", "strelka_somatic", "germline_haplotypecaller_shard"):
        assert profile[rule]["runtime"] == 1440

    variant_rules = (ROOT / "workflow/rules/variant_calling.smk").read_text()
    for rule in (
        "mutect2_sample_pon_chromosome",
        "mutect2_pon_shard",
        "mutect2_chromosome",
    ):
        body = variant_rules.split(f"rule {rule}:", 1)[1].split("\nrule ", 1)[0]
        assert "runtime=720" in body
    strelka = variant_rules.split("rule strelka_somatic:", 1)[1].split("\nrule ", 1)[0]
    assert 'runtime=1440 if ANALYSIS_TYPE == "wgs" else 360' in strelka
    alignment = (ROOT / "workflow/rules/alignment.smk").read_text()
    assert "runtime=1440" in alignment.split("rule bwa_mem_paired:", 1)[1].split("\nrule ", 1)[0]
    germline = (ROOT / "workflow/rules/germline.smk").read_text()
    assert (
        "runtime=1440"
        in germline.split("rule germline_haplotypecaller_shard:", 1)[1].split("\nrule ", 1)[0]
    )


def test_memory_headroom_for_jvm_and_observed_samtools_workloads():
    profile = yaml.safe_load((ROOT / "config/slurm/config.yaml").read_text())["set-resources"]
    expected = {
        "fastqc": 4096,
        "samtools_alignment_qc": 8192,
        "split_wgs_intervals": 8192,
        "contamination_sites": 8192,
        "mutect2_orientation_model": 8192,
    }
    for rule, mem_mb in expected.items():
        assert profile[rule]["mem_mb"] == mem_mb

    rule_sources = {
        "fastqc": ROOT / "workflow/rules/fastqc.smk",
        "samtools_alignment_qc": ROOT / "workflow/rules/qc.smk",
        "split_wgs_intervals": ROOT / "workflow/rules/variant_calling.smk",
        "contamination_sites": ROOT / "workflow/rules/variant_calling.smk",
        "mutect2_orientation_model": ROOT / "workflow/rules/variant_calling.smk",
    }
    for rule, path in rule_sources.items():
        body = path.read_text().split(f"rule {rule}:", 1)[1].split("\nrule ", 1)[0]
        assert f"mem_mb={expected[rule]}" in body

    # FastQC uses 512 MB per requested thread; the three GATK rules use -Xmx4g.
    fastqc_body = rule_sources["fastqc"].read_text().split("rule fastqc:", 1)[1]
    assert "threads: 4" in fastqc_body
    assert profile["fastqc"]["mem_mb"] >= 4 * 512 + 2048
    for rule in ("split_wgs_intervals", "contamination_sites", "mutect2_orientation_model"):
        body = rule_sources[rule].read_text().split(f"rule {rule}:", 1)[1].split("\nrule ", 1)[0]
        assert 'gatk_java_options("4g")' in body
        assert profile[rule]["mem_mb"] >= 4096 + 4096


def test_alignment_inputs_are_ignored_by_git():
    ignore = (ROOT / ".gitignore").read_text().splitlines()
    assert "/data/" in ignore
    for pattern in ("*.bam", "*.bai", "*.cram", "*.crai"):
        assert pattern in ignore
