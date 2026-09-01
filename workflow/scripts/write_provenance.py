#!/usr/bin/env python3
"""Write a config-only provenance table for the somatic pipeline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from config_utils import deep_merge

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - pipeline targets Python 3.11+
    tomllib = None


Row = tuple[str, str, str, str]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_yaml(path: Path) -> Any:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def read_toml(path: Path) -> dict[str, Any]:
    if tomllib is None:
        raise RuntimeError("tomllib is required; run with Python 3.11 or newer")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def value_to_tsv(value: Any) -> str:
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def add(rows: list[Row], key: str, value: Any, source: str, note: str = "") -> None:
    rows.append((key, value_to_tsv(value), source, note))


def add_mapping(
    rows: list[Row],
    prefix: str,
    mapping: dict[str, Any],
    source: str,
    notes: dict[str, str] | None = None,
) -> None:
    notes = notes or {}
    for key in sorted(mapping):
        add(rows, f"{prefix}.{key}", mapping[key], source, notes.get(key, ""))


def add_nested_mapping(rows: list[Row], prefix: str, value: Any, source: str) -> None:
    if isinstance(value, dict):
        for key in sorted(value):
            add_nested_mapping(rows, f"{prefix}.{key}", value[key], source)
    else:
        add(rows, prefix, value, source)


def sample_summary(samples: dict[str, Any]) -> dict[str, Any]:
    sample_ids = sorted(samples.get("samples", {}))
    comparisons = {sample: values["baseline"] for sample, values in sorted(samples.get("comparisons", {}).items())}
    organoid_ids = sorted(comparisons)
    baseline_ids = sorted(set(comparisons.values()))
    return {
        "ids": sample_ids,
        "organoid_ids": organoid_ids,
        "baseline_ids": baseline_ids,
        "organoid_baseline_comparisons": comparisons,
        "roles": {sample: samples["samples"][sample].get("role") for sample in sample_ids},
        "count": len(sample_ids),
        "organoid_count": len(organoid_ids),
        "baseline_count": len(baseline_ids),
    }


def add_dependency_rows(rows: list[Row], pixi: dict[str, Any], source: str) -> None:
    sections = {
        "dependency": pixi.get("dependencies", {}),
        "pypi_dependency": pixi.get("pypi-dependencies", {}),
    }
    for section_name, dependencies in sections.items():
        for name in sorted(dependencies):
            value = dependencies[name]
            note = f"declared {section_name.replace('_', ' ')}"
            if value == "*":
                note += "; exact resolved version not captured in config-only mode"
            add(rows, f"software.{name}.declared_dependency", value, source, note)


def build_rows(
    config: dict[str, Any],
    samples: dict[str, Any],
    pixi: dict[str, Any],
    slurm: dict[str, Any] | None,
    config_path: Path,
    samples_path: Path,
    pixi_path: Path,
    slurm_path: Path | None,
    run_path: Path,
    config_overlay_path: Path | None = None,
) -> list[Row]:
    config_source = str(config_path)
    if config_overlay_path:
        config_source = f"{config_source} + {config_overlay_path}"
    samples_source = str(samples_path)
    rows: list[Row] = []

    add(rows, "provenance.generated_at_utc", datetime.now(timezone.utc).isoformat(), "derived")
    add(rows, "provenance.version_capture_mode", "config_only", "derived")
    add(rows, "provenance.upstream_commit", "a533612", "repository specification")
    if config_overlay_path:
        add(rows, "provenance.config_overlay", str(config_overlay_path), "derived")
    add(rows, "run.path", str(run_path), "derived", "Website can derive batch identity from this path")
    add(rows, "samples.source_file", samples_source, "derived")

    add_mapping(rows, "pipeline", config.get("pipeline", {}), config_source)

    for section in ("analysis", "execution", "reference", "singularity", "storage", "output", "cluster"):
        if section in config:
            add_nested_mapping(rows, section, config[section], config_source)

    for scalar_key in ("container_runtime", "cram_version", "pixi_bin"):
        if scalar_key in config:
            add(rows, scalar_key, config[scalar_key], config_source)

    if "chromosomes" in config:
        add(rows, "chromosomes", config["chromosomes"], config_source)

    sample_rows = sample_summary(samples)
    for key in sorted(sample_rows):
        add(rows, f"samples.{key}", sample_rows[key], samples_source)

    for section in (
        "gatk",
        "picard",
        "strelka",
        "mutect2",
        "germline",
        "vt",
        "filtering",
        "mutational_signatures",
        "hotspots",
        "annotation",
        "coverage",
        "somalier",
    ):
        if section in config:
            add_nested_mapping(rows, section, config[section], config_source)

    for tool in ("gatk", "strelka", "picard"):
        tool_config = config.get(tool, {})
        for key in ("version", "image"):
            if key in tool_config:
                add(rows, f"software.{tool}.{key}", tool_config[key], config_source, "configured value")

    add_dependency_rows(rows, pixi, str(pixi_path))

    if slurm and slurm_path:
        add_nested_mapping(rows, "slurm_profile", slurm, str(slurm_path))

    reference = config.get("reference", {})
    checksum_paths = {
        "reference.fai_sha256": Path(f"{reference.get('genome', '')}.fai"),
        "reference.dictionary_sha256": Path(reference.get("genome_dict", "")),
        "reference.territory_sha256": Path(reference.get("wgs_calling_regions", "")),
        "reference.gnomad_sha256": Path(reference.get("gnomad", "")),
        "reference.population_vcf_sha256": Path(reference.get("population_vcf", "")),
        "reference.contamination_sites_sha256": Path(reference.get("contamination_sites", "")),
    }
    metadata = reference.get("regions_metadata")
    if metadata:
        checksum_paths["reference.capture_metadata_sha256"] = Path(metadata)
    for key, path in checksum_paths.items():
        if str(path) and path.is_file():
            add(rows, key, sha256(path), "derived")

    return rows


def write_rows(rows: list[Row], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["key", "value", "source", "note"])
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--config-overlay", type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--pixi", required=True, type=Path)
    parser.add_argument("--slurm-config", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    slurm = None
    if args.slurm_config and args.slurm_config.exists():
        slurm = read_yaml(args.slurm_config)

    resolved_config = read_yaml(args.config)
    if args.config_overlay:
        resolved_config = deep_merge(resolved_config, read_yaml(args.config_overlay))

    rows = build_rows(
        config=resolved_config,
        samples=read_yaml(args.samples),
        pixi=read_toml(args.pixi),
        slurm=slurm,
        config_path=args.config,
        samples_path=args.samples,
        pixi_path=args.pixi,
        slurm_path=args.slurm_config if slurm else None,
        run_path=Path.cwd().resolve(),
        config_overlay_path=args.config_overlay,
    )
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
