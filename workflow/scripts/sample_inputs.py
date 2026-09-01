"""Validation and deterministic resolution of organoid cohort inputs."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
ROLES = {"baseline", "organoid"}
READ_GROUP_FIELDS = ("platform", "library", "unit")


def _table(manifest: dict[str, Any], key: str, source: str) -> dict[str, Any]:
    value = manifest.get(key)
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{source} must define a non-empty '{key}' mapping")
    return value


def comparison_map(manifest: dict[str, Any], source: str) -> dict[str, str]:
    comparisons = _table(manifest, "comparisons", source)
    result = {}
    for organoid, values in comparisons.items():
        if not isinstance(values, dict) or set(values) != {"baseline"}:
            raise ValueError(f"Comparison '{organoid}' must contain exactly one 'baseline' field")
        baseline = values["baseline"]
        if not isinstance(baseline, str) or not baseline:
            raise ValueError(f"Comparison '{organoid}' baseline must be a sample identifier")
        result[organoid] = baseline
    return result


def validate_samples(manifest: dict[str, Any], source: str, check_readable: bool = False) -> None:
    samples = _table(manifest, "samples", source)
    comparisons = comparison_map(manifest, source)
    for sample, values in samples.items():
        if not isinstance(sample, str) or not SAFE_ID.fullmatch(sample):
            raise ValueError(
                f"Unsafe sample identifier {sample!r}; use letters, digits, '.', '_' or '-'"
            )
        if not isinstance(values, dict):
            raise ValueError(f"Sample '{sample}' must be a mapping")
        if values.get("role") not in ROLES:
            raise ValueError(f"Sample '{sample}' role must be one of: baseline, organoid")
        for key in ("donor", "lineage"):
            if not isinstance(values.get(key), str) or not values[key]:
                raise ValueError(f"Sample '{sample}' must define non-empty '{key}' metadata")
        fastq_keys = {"fastq_1", "fastq_2"} & set(values)
        has_fastqs = fastq_keys == {"fastq_1", "fastq_2"}
        has_cram = "cram" in values
        if fastq_keys and not has_fastqs:
            raise ValueError(f"Sample '{sample}' must define both 'fastq_1' and 'fastq_2'")
        if has_fastqs == has_cram:
            raise ValueError(
                f"Sample '{sample}' must define exactly one input type: paired FASTQs or CRAM"
            )
        for key in ("fastq_1", "fastq_2", "cram", "crai", "bam_sample"):
            if key in values and (not isinstance(values[key], str) or not values[key]):
                raise ValueError(f"Sample '{sample}' key '{key}' must be a non-empty string")
        if has_cram and "crai" in values:
            cram = Path(values["cram"])
            recognized_indexes = {Path(f"{cram}.crai"), cram.with_suffix(".crai")}
            if Path(values["crai"]) not in recognized_indexes:
                raise ValueError(
                    f"Sample '{sample}' CRAI must use a standard samtools-discoverable path: "
                    f"'{cram}.crai' or '{cram.with_suffix('.crai')}'"
                )
        if has_fastqs:
            read_group_data = values.get("read_group")
            if not isinstance(read_group_data, dict):
                raise ValueError(f"FASTQ sample '{sample}' must define read_group metadata")
            missing = [key for key in READ_GROUP_FIELDS if not read_group_data.get(key)]
            if missing:
                raise ValueError(f"FASTQ sample '{sample}' read_group lacks: {', '.join(missing)}")
        if check_readable:
            for key in ("fastq_1", "fastq_2", "cram"):
                if key in values and not os.access(values[key], os.R_OK):
                    raise ValueError(f"Sample '{sample}' {key} is not readable: {values[key]}")
            if has_cram and not os.access(aligned_cram_index(manifest, sample, source), os.R_OK):
                raise ValueError(f"Sample '{sample}' CRAM index is not readable")

    organoids = {name for name, values in samples.items() if values["role"] == "organoid"}
    baselines = {name for name, values in samples.items() if values["role"] == "baseline"}
    if set(comparisons) != organoids:
        missing = sorted(organoids - set(comparisons))
        extra = sorted(set(comparisons) - organoids)
        details = []
        if missing:
            details.append("organoids without comparisons: " + ", ".join(missing))
        if extra:
            details.append("non-organoid comparisons: " + ", ".join(extra))
        raise ValueError("; ".join(details))
    for organoid, baseline in comparisons.items():
        if baseline not in samples:
            raise ValueError(f"Baseline sample '{baseline}' for '{organoid}' is undefined")
        if baseline not in baselines:
            raise ValueError(f"Comparison baseline '{baseline}' does not have role: baseline")
        for key in ("donor", "lineage"):
            if samples[organoid][key] != samples[baseline][key]:
                raise ValueError(
                    f"Comparison '{organoid}' and baseline '{baseline}' have inconsistent {key}"
                )


def sample_config(samples: dict[str, Any], sample: str, source: str) -> dict[str, Any]:
    try:
        return samples["samples"][sample]
    except KeyError:
        raise ValueError(f"Sample '{sample}' is not defined in {source}") from None


def sample_has_fastqs(samples, sample, source):
    values = sample_config(samples, sample, source)
    return "fastq_1" in values and "fastq_2" in values


def sample_has_cram(samples, sample, source):
    return "cram" in sample_config(samples, sample, source)


def sample_has_final_vcf(samples, sample, source):
    return False


def final_vcf(samples, sample, source):
    raise ValueError("Final-VCF-only mode is not supported by the organoid workflow")


def final_vcf_samples(samples, source):
    return []


def vcf_only_mode(samples, source):
    return False


def matched_normal_samples(samples, source):
    return sorted(set(comparison_map(samples, source).values()))


def organoid_samples(samples, source):
    return sorted(comparison_map(samples, source))


def aligned_cram(samples, sample, source):
    values = sample_config(samples, sample, source)
    return values.get("cram", f"results/cram/{sample}.sorted.dups.cram")


def aligned_cram_index(samples, sample, source):
    values = sample_config(samples, sample, source)
    if "cram" in values:
        return values.get("crai", f"{values['cram']}.crai")
    return f"{aligned_cram(samples, sample, source)}.crai"


def bam_sample_name(samples, sample, source):
    return sample_config(samples, sample, source).get("bam_sample", sample)


def read_group(samples, sample, source):
    return dict(sample_config(samples, sample, source).get("read_group", {}))


def missing_sample_files(samples, source):
    missing = []
    for sample in sorted(samples.get("samples", {})):
        values = sample_config(samples, sample, source)
        for key in ("fastq_1", "fastq_2", "cram"):
            path = values.get(key)
            if path and not Path(path).is_file():
                missing.append((sample, key, path))
        if "cram" in values:
            index = aligned_cram_index(samples, sample, source)
            if not Path(index).is_file():
                missing.append((sample, "crai" if "crai" in values else "crai(default)", index))
    return missing
