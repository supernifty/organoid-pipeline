#!/usr/bin/env python3
"""Validate mixed FASTQ/CRAM sample input handling."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from workflow.scripts.sample_inputs import (
    aligned_cram,
    aligned_cram_index,
    bam_sample_name,
    final_vcf,
    final_vcf_samples,
    missing_sample_files,
    matched_normal_samples,
    sample_has_cram,
    sample_has_final_vcf,
    sample_has_fastqs,
    validate_samples,
    vcf_only_mode,
)


SOURCE = "test samples"


def expect_error(samples, message):
    try:
        validate_samples(samples, SOURCE)
    except ValueError as error:
        if message not in str(error):
            raise AssertionError(f"expected {message!r} in {error!s}") from error
    else:
        raise AssertionError("expected sample validation to fail")


def main() -> None:
    samples = {
        "samples": {
            "T": {
                "cram": "/data/crams/T.sorted.dups.cram",
                "crai": "/data/crams/T.sorted.dups.cram.crai",
                "bam_sample": "T_HEADER",
            },
            "N": {
                "fastq_1": "data/N_R1.fastq.gz",
                "fastq_2": "data/N_R2.fastq.gz",
            },
            "C": {
                "cram": "/data/crams/C.sorted.dups.cram",
            },
        },
        "tumours": {"T": "N"},
    }

    validate_samples(samples, SOURCE)
    assert sample_has_cram(samples, "T", SOURCE)
    assert not sample_has_fastqs(samples, "T", SOURCE)
    assert sample_has_fastqs(samples, "N", SOURCE)
    assert not sample_has_cram(samples, "N", SOURCE)
    assert aligned_cram(samples, "T", SOURCE) == "/data/crams/T.sorted.dups.cram"
    assert aligned_cram_index(samples, "T", SOURCE) == "/data/crams/T.sorted.dups.cram.crai"
    assert bam_sample_name(samples, "T", SOURCE) == "T_HEADER"
    assert bam_sample_name(samples, "N", SOURCE) == "N"
    assert aligned_cram(samples, "N", SOURCE) == "results/cram/N.sorted.dups.cram"
    assert aligned_cram_index(samples, "N", SOURCE) == "results/cram/N.sorted.dups.cram.crai"
    assert aligned_cram_index(samples, "C", SOURCE) == "/data/crams/C.sorted.dups.cram.crai"
    assert missing_sample_files(samples, SOURCE) == [
        ("C", "cram", "/data/crams/C.sorted.dups.cram"),
        ("C", "crai(default)", "/data/crams/C.sorted.dups.cram.crai"),
        ("N", "fastq_1", "data/N_R1.fastq.gz"),
        ("N", "fastq_2", "data/N_R2.fastq.gz"),
        ("T", "cram", "/data/crams/T.sorted.dups.cram"),
        ("T", "crai", "/data/crams/T.sorted.dups.cram.crai"),
    ]
    assert matched_normal_samples(samples, SOURCE) == ["N"]

    shared_normal = {
        "samples": {
            "T1": {"cram": "T1.cram"},
            "T2": {"cram": "T2.cram"},
            "N": {"cram": "N.cram"},
            "UNPAIRED": {"cram": "UNPAIRED.cram"},
        },
        "tumours": {"T2": "N", "T1": "N"},
    }
    assert matched_normal_samples(shared_normal, SOURCE) == ["N"]

    expect_error({"samples": {"S": {"fastq_1": "S_R1.fastq.gz"}}}, "both 'fastq_1' and 'fastq_2'")
    expect_error({"samples": {"S": {}}}, "either 'final_vcf', 'cram', or both 'fastq_1' and 'fastq_2'")
    expect_error({"samples": {"S": {"cram": "S.cram", "bam_sample": 1}}}, "key 'bam_sample' must be a string")
    expect_error({"samples": {"N": {"cram": "N.cram"}}, "tumours": {"T": "N"}}, "Tumour sample 'T'")
    expect_error({"samples": {"T": {"cram": "T.cram"}}, "tumours": {"T": "N"}}, "Normal sample 'N'")

    vcf_samples = {"samples": {"T": {"final_vcf": "variants/T.vcf.gz"}}}
    validate_samples(vcf_samples, SOURCE)
    assert sample_has_final_vcf(vcf_samples, "T", SOURCE)
    assert not sample_has_cram(vcf_samples, "T", SOURCE)
    assert not sample_has_fastqs(vcf_samples, "T", SOURCE)
    assert final_vcf(vcf_samples, "T", SOURCE) == "variants/T.vcf.gz"
    assert final_vcf_samples(vcf_samples, SOURCE) == ["T"]
    assert vcf_only_mode(vcf_samples, SOURCE)
    assert matched_normal_samples(vcf_samples, SOURCE) == []
    assert missing_sample_files(vcf_samples, SOURCE) == [("T", "final_vcf", "variants/T.vcf.gz")]

    expect_error(
        {"samples": {"T": {"final_vcf": "T.vcf", "cram": "T.cram"}}},
        "must not mix 'final_vcf'",
    )
    expect_error(
        {"samples": {"T": {"final_vcf": "T.vcf"}}, "tumours": {"T": "N"}},
        "must not define 'tumours'",
    )
    expect_error(
        {"samples": {"T": {"final_vcf": "T.vcf"}, "N": {"cram": "N.cram"}}},
        "must not mix final_vcf samples",
    )


if __name__ == "__main__":
    main()
