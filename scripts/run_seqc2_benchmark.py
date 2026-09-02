#!/usr/bin/env python3
"""Plan or explicitly execute the storage-intensive SEQC2 benchmark matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def absolute_paths(config):
    """Resolve pipeline resource/image paths before Snakemake changes directory."""
    for section, keys in {
        "reference": (
            "genome",
            "genome_dict",
            "wgs_calling_regions",
            "wgs_exclude_regions",
            "gnomad",
            "population_vcf",
            "contamination_sites",
            "panel_of_normals",
            "problematic_regions",
            "low_mappability_regions",
            "repeat_regions",
        ),
        "gatk": ("singularity_image",),
        "picard": ("singularity_image",),
        "strelka": ("singularity_image",),
    }.items():
        for key in keys:
            value = config.get(section, {}).get(key)
            if value and not Path(value).is_absolute():
                config[section][key] = str((ROOT / value).resolve())
    return config


def require_file(value, label):
    path = Path(value or "")
    if not value or not path.is_file():
        raise ValueError(f"{label} must name a readable local file")
    return path.resolve()


def verify_checksum(path, expected, label):
    if not expected or len(expected) != 64:
        raise ValueError(f"{label}_sha256 must be a reviewed SHA-256 digest")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest().lower() != expected.lower():
        raise ValueError(f"Checksum mismatch for {label}: {path}")


def command_text(command):
    return " ".join(str(value) for value in command)


def add_cosines(matrix_path, metric_paths):
    lines = Path(matrix_path).read_text().splitlines()
    names = lines[0].split("\t")[1:]
    vectors = {name: [] for name in names}
    for line in lines[1:]:
        values = line.split("\t")[1:]
        for name, value in zip(names, values, strict=True):
            vectors[name].append(float(value))
    truth = vectors["truth"]
    for name, path in metric_paths.items():
        vector = vectors[name]
        dot = sum(left * right for left, right in zip(truth, vector, strict=True))
        denominator = math.sqrt(sum(value * value for value in truth)) * math.sqrt(
            sum(value * value for value in vector)
        )
        payload = json.loads(Path(path).read_text())
        payload["sbs96_cosine_similarity"] = dot / denominator if denominator else 0.0
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--cores", type=int, default=8)
    args = parser.parse_args()
    settings = yaml.safe_load(args.config.read_text())
    if not settings.get("enabled"):
        raise ValueError("Set enabled: true after reviewing storage, inputs, and checksums")
    tumour = require_file(settings.get("tumour_alignment"), "tumour_alignment")
    baseline = require_file(settings.get("baseline_alignment"), "baseline_alignment")
    reference = require_file(settings.get("reference"), "reference")
    territory = require_file(settings.get("high_confidence_regions"), "high_confidence_regions")
    truth = require_file(settings.get("truth_snv"), "truth_snv")
    verify_checksum(tumour, settings.get("tumour_alignment_sha256"), "tumour_alignment")
    verify_checksum(baseline, settings.get("baseline_alignment_sha256"), "baseline_alignment")
    verify_checksum(reference, settings.get("reference_sha256"), "reference")
    verify_checksum(truth, settings.get("truth_snv_sha256"), "truth_snv")
    verify_checksum(
        territory, settings.get("high_confidence_regions_sha256"), "high_confidence_regions"
    )
    pipeline_config_path = require_file(settings.get("pipeline_config"), "pipeline_config")
    pipeline_config = absolute_paths(yaml.safe_load(pipeline_config_path.read_text()))
    reference_dict = require_file(
        pipeline_config.get("reference", {}).get("genome_dict"),
        "pipeline reference.genome_dict",
    )
    output_root = Path(settings.get("output_directory", "benchmarks/seqc2/runs")).resolve()
    seeds = [int(value) for value in settings["seeds"]]
    baseline_depths = [float(value) for value in settings["baseline_target_depths"]]
    if settings.get("normal_normal_negative_control", True) and 6.0 not in baseline_depths:
        raise ValueError("normal-normal control requires baseline_target_depths to include 6")
    tumour_target = float(settings.get("tumour_target_depth", 6))
    tumour_depth = float(settings["tumour_input_depth"])
    baseline_depth = float(settings["baseline_input_depth"])
    estimated = sum(tumour.stat().st_size * tumour_target / tumour_depth for _ in seeds)
    estimated += sum(
        baseline.stat().st_size * depth / baseline_depth for _ in seeds for depth in baseline_depths
    )
    if settings.get("normal_normal_negative_control", True):
        estimated += sum(baseline.stat().st_size * 6 / baseline_depth for _ in seeds)
    free = shutil.disk_usage(output_root.parent if output_root.parent.exists() else ROOT).free
    print(f"Estimated downsampled alignment storage: {estimated / 2**30:.1f} GiB")
    print(f"Available storage: {free / 2**30:.1f} GiB")
    if estimated * 1.5 > free:
        raise ValueError("Insufficient free space after applying a 1.5× working-space margin")

    commands = []
    cases = []
    for seed in seeds:
        tumour_out = output_root / "downsampled" / f"HCC1395.seed{seed}.{tumour_target:g}x.cram"
        tumour_command = [
            sys.executable,
            str(ROOT / "scripts/downsample_alignment.py"),
            "--input",
            str(tumour),
            "--reference",
            str(reference),
            "--reference-dict",
            str(reference_dict),
            "--territory",
            str(territory),
            "--input-depth",
            str(tumour_depth),
            "--target-depth",
            str(tumour_target),
            "--seed",
            str(seed),
            "--output",
            str(tumour_out),
            "--report",
            f"{tumour_out}.json",
        ]
        commands.append(tumour_command)
        for target in baseline_depths:
            baseline_out = output_root / "downsampled" / f"HCC1395BL.seed{seed}.{target:g}x.cram"
            commands.append(
                [
                    sys.executable,
                    str(ROOT / "scripts/downsample_alignment.py"),
                    "--input",
                    str(baseline),
                    "--reference",
                    str(reference),
                    "--reference-dict",
                    str(reference_dict),
                    "--territory",
                    str(territory),
                    "--input-depth",
                    str(baseline_depth),
                    "--target-depth",
                    str(target),
                    "--seed",
                    str(seed),
                    "--output",
                    str(baseline_out),
                    "--report",
                    f"{baseline_out}.json",
                ]
            )
            cases.append(
                (
                    f"seed{seed}.tumour{tumour_target:g}x.baseline{target:g}x",
                    tumour_out,
                    baseline_out,
                )
            )
        if settings.get("normal_normal_negative_control", True):
            control_out = output_root / "downsampled" / f"HCC1395BL.control.seed{seed}.6x.cram"
            commands.append(
                [
                    sys.executable,
                    str(ROOT / "scripts/downsample_alignment.py"),
                    "--input",
                    str(baseline),
                    "--reference",
                    str(reference),
                    "--reference-dict",
                    str(reference_dict),
                    "--territory",
                    str(territory),
                    "--input-depth",
                    str(baseline_depth),
                    "--target-depth",
                    "6",
                    "--seed",
                    str(seed + 100000),
                    "--output",
                    str(control_out),
                    "--report",
                    f"{control_out}.json",
                ]
            )
            baseline_six = output_root / "downsampled" / f"HCC1395BL.seed{seed}.6x.cram"
            cases.append((f"seed{seed}.normal-normal.6x", control_out, baseline_six))

    for command in commands:
        print(command_text(command))
        if args.execute and not Path(command[command.index("--output") + 1]).exists():
            subprocess.run(command, check=True)
    for name, tumour_cram, baseline_cram in cases:
        case = output_root / name
        manifest = {
            "samples": {
                "HCC1395BL": {
                    "role": "baseline",
                    "donor": "SEQC2",
                    "lineage": "HCC1395",
                    "cram": str(baseline_cram),
                    "crai": f"{baseline_cram}.crai",
                },
                "HCC1395": {
                    "role": "organoid",
                    "donor": "SEQC2",
                    "lineage": "HCC1395",
                    "cram": str(tumour_cram),
                    "crai": f"{tumour_cram}.crai",
                },
            },
            "comparisons": {"HCC1395": {"baseline": "HCC1395BL"}},
        }
        config = yaml.safe_load(yaml.safe_dump(pipeline_config))
        config["reference"]["wgs_calling_regions"] = str(territory)
        config["run_management"] = {
            "samples_file": str((case / "samples.yaml").resolve()),
            "config_file": str((case / "config.yaml").resolve()),
        }
        run = [
            "pixi",
            "run",
            "snakemake",
            "--snakefile",
            str(ROOT / "Snakefile"),
            "--directory",
            str(case),
            "--cores",
            str(args.cores),
            "--configfile",
            str((case / "config.yaml").resolve()),
        ]
        print(command_text(run))
        if args.execute:
            case.mkdir(parents=True, exist_ok=True)
            (case / "samples.yaml").write_text(yaml.safe_dump(manifest, sort_keys=False))
            (case / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
            subprocess.run(run, cwd=ROOT, check=True)
            callable_bases = sum(
                int(row.split()[2]) - int(row.split()[1])
                for row in territory.read_text().splitlines()
                if row and not row.startswith("#")
            )
            tier_paths = {
                "mutect2": "results/callers/HCC1395.mutect2.pass.vcf.gz",
                "strelka2": "results/callers/HCC1395.strelka.pass.vcf.gz",
                "intersection": "results/callers/HCC1395.intersection.vcf.gz",
                "stringent": "results/catalogs/HCC1395.stringent.vcf.gz",
            }
            metric_paths = {}
            for tier, relative in tier_paths.items():
                metric_path = case / f"{tier}.metrics.json"
                metric_paths[tier] = metric_path
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT / "workflow/scripts/benchmark_metrics.py"),
                        "--truth",
                        str(truth),
                        "--calls",
                        str(case / relative),
                        "--callable-bases",
                        str(callable_bases),
                        "--output",
                        str(metric_path),
                        "--false-positive-vcf",
                        str(case / f"{tier}.false-positive.vcf"),
                    ],
                    check=True,
                )
            matrix_path = case / "benchmark.sbs96.tsv"
            matrix_command = [
                sys.executable,
                str(ROOT / "workflow/scripts/sbs96.py"),
                "--reference",
                str(reference),
                "--catalog",
                f"truth={truth}",
            ]
            for tier, relative in tier_paths.items():
                matrix_command.extend(("--catalog", f"{tier}={case / relative}"))
                matrix_command.extend(
                    ("--catalog", f"{tier}_false_positive={case / f'{tier}.false-positive.vcf'}")
                )
            matrix_command.extend(("--output", str(matrix_path)))
            subprocess.run(matrix_command, check=True)
            add_cosines(matrix_path, metric_paths)
    summary = {
        "cases": [name for name, _, _ in cases],
        "estimated_bytes": estimated,
        "executed": args.execute,
        "normal_normal_negative_control": settings.get("normal_normal_negative_control", True),
    }
    if args.execute:
        (output_root / "benchmark_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
