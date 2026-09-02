#!/usr/bin/env python3
"""Summarize an isolated run-manager batch without conflating wall and elapsed time."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import yaml

VERSION = "1.0.0"
ALIASES = {
    "wall_seconds": ("s", "seconds", "wall_seconds", "walltime_seconds", "runtime_seconds"),
    "cpu_seconds": ("cpu_time", "cpu_seconds", "cpu_time_seconds"),
    "max_rss_mb": ("max_rss", "max_rss_mb", "max_memory_mb"),
    "io_in_mb": ("io_in", "io_in_mb", "input_mb"),
    "io_out_mb": ("io_out", "io_out_mb", "output_mb"),
}
GROUP_PATTERNS = (
    ("preparation", ("downsample", "prepare_fastq")),
    ("FastQC", ("fastqc",)),
    ("BWA/sort", (".bwa", "bwa_", "/bwa")),
    ("MarkDuplicates", ("mark_duplicates", "duplicates")),
    ("alignment QC", ("samtools", "coverage", "alignment_metrics", "insert_size", "preflight")),
    ("Mutect2", ("mutect2", "contamination")),
    ("Strelka2", ("strelka",)),
    (
        "recounting/filtering",
        ("caller_tier", "cohort_", "catalog", "filter", "intersect", "normalize"),
    ),
    ("SBS96", ("sbs96", "signature")),
    ("aggregate QC", ("multiqc", "aggregate", "final_variant_counts", "provenance", "manifest")),
)


class ReportError(RuntimeError):
    pass


def timestamp(value: str | None) -> dt.datetime | None:
    if not value:
        return None
    try:
        return dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None


def number(row: dict[str, str], aliases: tuple[str, ...]) -> float | None:
    for alias in aliases:
        value = row.get(alias)
        if value not in (None, "", "NA", "N/A", "."):
            try:
                result = float(value)
                return result if math.isfinite(result) else None
            except ValueError:
                continue
    return None


def benchmark_identity(path: Path, benchmark_root: Path) -> tuple[str, str | None]:
    relative = path.relative_to(benchmark_root)
    rule = relative.parent.as_posix() if relative.parent != Path(".") else relative.stem
    stem = relative.name.removesuffix(".tsv")
    sample = stem.split(".", 1)[0] if "." in stem else None
    return rule, sample


def group_for(rule: str, path: str) -> str:
    value = f"{rule}/{path}".lower()
    for group, patterns in GROUP_PATTERNS:
        if any(pattern.lower() in value for pattern in patterns):
            return group
    return "other"


def read_benchmarks(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not root.exists():
        return rows
    for path in sorted(root.rglob("*.tsv")):
        try:
            with path.open(newline="") as handle:
                parsed = list(csv.DictReader(handle, delimiter="\t"))
        except (OSError, csv.Error, UnicodeDecodeError):
            continue
        rule, sample = benchmark_identity(path, root)
        for index, row in enumerate(parsed, 1):
            item = {
                "rule": rule,
                "sample": sample,
                "group": group_for(rule, str(path)),
                "benchmark_path": str(path),
                "row": index,
                "completion_state": "complete",
            }
            item.update({key: number(row, aliases) for key, aliases in ALIASES.items()})
            rows.append(item)
    return rows


def preparation_rows(paths: list[Path], reports: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for path, report in zip(paths, reports, strict=True):
        start = timestamp(report.get("started_at"))
        end = timestamp(report.get("completed_at"))
        pairs = report.get("input_pairs")
        bases = report.get("input_sequenced_bases")
        rows.append(
            {
                "rule": "downsample_fastq_pair",
                "sample": report.get("sample") or report.get("reuse_identity", {}).get("sample"),
                "group": "preparation",
                "benchmark_path": str(path),
                "row": 1,
                "completion_state": report.get("status", "unavailable"),
                "wall_seconds": (end - start).total_seconds() if start and end else None,
                "cpu_seconds": None,
                "max_rss_mb": None,
                "io_in_mb": sum(
                    item.get("size", 0)
                    for item in report.get("reuse_identity", {}).get("inputs", [])
                )
                / 1024**2
                if report.get("reuse_identity", {}).get("inputs")
                else None,
                "io_out_mb": sum(item.get("size", 0) for item in report.get("outputs", []))
                / 1024**2
                if report.get("outputs")
                else None,
                "mean_input_read_length": bases / (2 * pairs) if pairs and bases else None,
            }
        )
    return rows


def aggregate(rows: Iterable[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "jobs": 0,
            "summed_job_wall_seconds": 0.0,
            "cpu_seconds": 0.0,
            "max_rss_mb": None,
            "io_in_mb": 0.0,
            "io_out_mb": 0.0,
        }
    )
    availability: dict[str, dict[str, bool]] = defaultdict(lambda: defaultdict(bool))
    for row in rows:
        name = row.get(key) or "unavailable"
        item = grouped[name]
        item["jobs"] += 1
        for metric in ("wall_seconds", "cpu_seconds", "io_in_mb", "io_out_mb"):
            if row.get(metric) is not None:
                target = "summed_job_wall_seconds" if metric == "wall_seconds" else metric
                item[target] += row[metric]
                availability[name][target] = True
        if row.get("max_rss_mb") is not None:
            current = item["max_rss_mb"]
            item["max_rss_mb"] = (
                row["max_rss_mb"] if current is None else max(current, row["max_rss_mb"])
            )
    for name, item in grouped.items():
        for metric in ("summed_job_wall_seconds", "cpu_seconds", "io_in_mb", "io_out_mb"):
            if not availability[name][metric]:
                item[metric] = None
        item[key] = name
    return [grouped[name] for name in sorted(grouped)]


def load_json(path: Path, required: bool = True) -> dict[str, Any]:
    if not path.is_file():
        if required:
            raise ReportError(f"required JSON file is missing: {path}")
        return {}
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(f"could not read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ReportError(f"expected a JSON object: {path}")
    return value


def elapsed_metrics(batch: dict[str, Any]) -> dict[str, Any]:
    launches = batch.get("launches", [])
    starts = [timestamp(item.get("submitted_at") or item.get("prepared_at")) for item in launches]
    ends = [
        timestamp(item.get("complete_at") or item.get("partial_at") or item.get("failed_at"))
        for item in launches
    ]
    starts = [value for value in starts if value]
    ends = [value for value in ends if value]
    elapsed = (max(ends) - min(starts)).total_seconds() if starts and ends else None
    calibration = [item for item in launches if item.get("targets")]
    full = [item for item in launches if item.get("default_dag")]

    def span(items: list[dict[str, Any]]) -> float | None:
        values = []
        for item in items:
            start = timestamp(item.get("submitted_at") or item.get("prepared_at"))
            end = timestamp(
                item.get("complete_at") or item.get("partial_at") or item.get("failed_at")
            )
            if start and end:
                values.append((end - start).total_seconds())
        return sum(values) if values else None

    return {
        "experienced_controller_elapsed_seconds": elapsed,
        "calibration_elapsed_seconds": span(calibration),
        "full_dag_elapsed_seconds": span(full),
        "scheduler_wait_seconds": None,
    }


def read_depths(results: Path) -> dict[str, float]:
    depths = {}
    for path in sorted((results / "qc/coverage").glob("*.wgs_coverage_mqc.tsv")):
        lines = [
            line for line in path.read_text().splitlines() if line and not line.startswith("#")
        ]
        if len(lines) < 2:
            continue
        row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=False))
        try:
            depths[row["sample"]] = float(row["mean_autosomal_depth"])
        except (KeyError, ValueError):
            continue
    return depths


def scheduler_accounting(job_ids: list[str]) -> list[dict[str, str]]:
    if not job_ids:
        return []
    command = [
        "sacct",
        "-P",
        "-n",
        "-j",
        ",".join(job_ids),
        "-o",
        "JobIDRaw,Submit,Start,End,State,ElapsedRaw,TotalCPU,MaxRSS",
    ]
    try:
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
    except OSError as exc:
        raise ReportError(f"could not execute sacct: {exc}") from exc
    if completed.returncode:
        raise ReportError(completed.stderr.strip() or "sacct failed")
    fields = (
        "job_id",
        "submit",
        "start",
        "end",
        "state",
        "elapsed_seconds",
        "total_cpu",
        "max_rss",
    )
    return [
        dict(zip(fields, line.split("|")[: len(fields)], strict=False))
        for line in completed.stdout.splitlines()
        if line.strip()
    ]


def optional_scheduler_accounting(
    requested: bool, job_ids: list[str]
) -> tuple[list[dict[str, str]], str | None]:
    if not requested:
        return [], None
    try:
        return scheduler_accounting(job_ids), None
    except ReportError as exc:
        return [], str(exc)


def scheduler_wait_seconds(rows: list[dict[str, str]]) -> float | None:
    waits = []
    for row in rows:
        submitted = timestamp(row.get("submit"))
        started = timestamp(row.get("start"))
        if submitted and started and started >= submitted:
            waits.append((started - submitted).total_seconds())
    return sum(waits) if waits else None


def output_storage(results: Path) -> dict[str, Any]:
    files = [path for path in results.rglob("*") if path.is_file()]
    categories = {}
    for directory in ("cram", "variants", "callers", "catalogs", "signatures", "qc", "aggregate"):
        selected = [path for path in files if path.relative_to(results).parts[0] == directory]
        categories[directory] = {
            "file_count": len(selected),
            "total_bytes": sum(path.stat().st_size for path in selected),
        }
    return {
        "file_count": len(files),
        "total_bytes": sum(path.stat().st_size for path in files),
        "categories": categories,
    }


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "rule",
        "sample",
        "group",
        "row",
        "completion_state",
        "wall_seconds",
        "cpu_seconds",
        "max_rss_mb",
        "io_in_mb",
        "io_out_mb",
        "benchmark_path",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "unavailable" if row.get(key) is None else row.get(key) for key in fields}
            )


def markdown(payload: dict[str, Any]) -> str:
    timing = payload["timing"]
    lines = [
        "# Run performance report",
        "",
        f"Batch: `{payload['batch']}`",
        "",
        "## Timing",
        "",
        "| Metric | Seconds |",
        "|---|---:|",
    ]
    labels = (
        ("experienced_controller_elapsed_seconds", "Experienced queue-inclusive elapsed"),
        ("summed_job_wall_seconds", "Summed job wall time"),
        ("summed_cpu_seconds", "Summed CPU time"),
        ("scheduler_wait_seconds", "Scheduler waiting"),
    )
    for key, label in labels:
        value = timing.get(key)
        lines.append(f"| {label} | {value if value is not None else 'unavailable'} |")
    lines.extend(
        [
            "",
            "## Grouped compute",
            "",
            "| Group | Jobs | Wall s | CPU s | Peak RSS MB |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for item in payload["grouped_costs"]:
        values = [
            item.get(key) if item.get(key) is not None else "unavailable"
            for key in ("group", "jobs", "summed_job_wall_seconds", "cpu_seconds", "max_rss_mb")
        ]
        lines.append("| " + " | ".join(map(str, values)) + " |")
    lines.extend(
        [
            "",
            "## Coverage and storage",
            "",
            f"Achieved mean autosomal depths: `{json.dumps(payload['achieved_depths'], sort_keys=True)}`  ",
            f"Final results storage: {payload['storage']['total_bytes']} bytes across {payload['storage']['file_count']} files.",
            "",
            "Cluster, queue, filesystem, sample, and concurrency choices affect these results. "
            "Unavailable metrics are not estimated.",
            "",
        ]
    )
    return "\n".join(lines)


def resolve_batch(value: Path) -> Path:
    if value.is_dir():
        return value.resolve()
    candidate = Path("runs") / value
    if candidate.is_dir():
        return candidate.resolve()
    raise ReportError(f"batch directory does not exist: {value}")


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--batch", required=True, type=Path)
    value.add_argument("--preparation-report", action="append", required=True, type=Path)
    value.add_argument("--calibration-history", required=True, type=Path)
    value.add_argument("--output-prefix", required=True, type=Path)
    value.add_argument(
        "--sacct", action="store_true", help="explicitly enrich controller scheduler timing"
    )
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    directory = resolve_batch(args.batch)
    batch = load_json(directory / "batch.json")
    preparation = [load_json(path) for path in args.preparation_report]
    calibration = load_json(args.calibration_history)
    config_path = directory / "config/current/config.yaml"
    config = yaml.safe_load(config_path.read_text()) if config_path.is_file() else {}
    manifest = load_json(directory / "results/analysis_manifest.json", required=False)
    rows = preparation_rows(args.preparation_report, preparation) + read_benchmarks(
        directory / "results/benchmarks"
    )
    timing = elapsed_metrics(batch)
    timing["summed_job_wall_seconds"] = (
        sum(row["wall_seconds"] for row in rows if row["wall_seconds"] is not None)
        if any(row["wall_seconds"] is not None for row in rows)
        else None
    )
    timing["summed_cpu_seconds"] = (
        sum(row["cpu_seconds"] for row in rows if row["cpu_seconds"] is not None)
        if any(row["cpu_seconds"] is not None for row in rows)
        else None
    )
    job_ids = sorted(
        {
            str(job)
            for launch in batch.get("launches", [])
            for job in launch.get("scheduler_ids", [])
        }
    )
    scheduler, scheduler_error = optional_scheduler_accounting(args.sacct, job_ids)
    if scheduler:
        timing["scheduler_wait_seconds"] = scheduler_wait_seconds(scheduler)
    territory_bases = manifest.get("territory_bases")
    wgs_config = config.get("analysis", {}).get("wgs", {})
    configured_shards = wgs_config.get("scatter_count")
    shard_count = configured_shards
    if shard_count is None and territory_bases and wgs_config.get("target_bases_per_shard"):
        shard_count = math.ceil(territory_bases / wgs_config["target_bases_per_shard"])
    payload = {
        "schema_version": 1,
        "tool_version": VERSION,
        "batch": batch.get("batch", directory.name),
        "batch_directory": str(directory),
        "git": next(
            (
                launch.get("git")
                for launch in reversed(batch.get("launches", []))
                if launch.get("git")
            ),
            None,
        ),
        "reference": {
            "build": config.get("reference", {}).get("build"),
            "identity": manifest.get("reference"),
            "callable_bases": next(
                (
                    item.get("reuse_identity", {}).get("callable_bases")
                    for item in preparation
                    if item
                ),
                None,
            ),
        },
        "preparation": preparation,
        "calibration": calibration,
        "achieved_depths": read_depths(directory / "results"),
        "configuration": {
            "containers": {
                key: config.get(key) for key in ("gatk", "strelka", "vep") if key in config
            },
            "mutect2_shards": shard_count,
            "mutect2_concurrency": wgs_config.get("max_concurrent_mutect2_shards"),
            "scheduler_friendly_runtime_limits_minutes": {
                "bwa_mem_paired": 1440,
                "mutect2_chromosome": 720,
                "strelka_somatic": 1440,
            },
        },
        "controller_launches": batch.get("launches", []),
        "timing": timing,
        "benchmarks": rows,
        "per_rule": aggregate(rows, "rule"),
        "grouped_costs": aggregate(rows, "group"),
        "maximum_observed_rss_mb": max(
            (row["max_rss_mb"] for row in rows if row["max_rss_mb"] is not None), default=None
        ),
        "scheduler_accounting": scheduler,
        "scheduler_accounting_error": scheduler_error,
        "storage": output_storage(directory / "results"),
        "limitations": "Operational benchmark only; cluster, queue, filesystem, sample, and concurrency dependent; no truth set.",
    }
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    json_path = args.output_prefix.with_suffix(".json")
    tsv_path = args.output_prefix.with_suffix(".tsv")
    md_path = args.output_prefix.with_suffix(".md")
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_tsv(tsv_path, rows)
    md_path.write_text(markdown(payload))
    print(
        json.dumps(
            {"json": str(json_path), "tsv": str(tsv_path), "markdown": str(md_path)}, sort_keys=True
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReportError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
