from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import summarize_run_performance as summary


def write_benchmark(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def test_benchmark_aliases_repeated_rows_and_aggregation(tmp_path: Path) -> None:
    root = tmp_path / "benchmarks"
    write_benchmark(
        root / "alignment/CASE.bwa.tsv",
        ["s", "max_rss", "io_in", "io_out", "cpu_time"],
        [[10, 100, 1, 2, 20], [12, 110, 3, 4, 24]],
    )
    write_benchmark(root / "mutect2/CASE.0001.tsv", ["wall_seconds", "max_rss_mb"], [[5, 200]])
    rows = summary.read_benchmarks(root)
    assert len(rows) == 3
    assert rows[0]["group"] == "BWA/sort"
    assert rows[0]["sample"] == "CASE"
    grouped = {item["group"]: item for item in summary.aggregate(rows, "group")}
    assert grouped["BWA/sort"]["summed_job_wall_seconds"] == 22
    assert grouped["BWA/sort"]["cpu_seconds"] == 44
    assert grouped["BWA/sort"]["max_rss_mb"] == 110
    assert grouped["Mutect2"]["cpu_seconds"] is None


def test_preparation_rows_include_elapsed_read_length_and_io(tmp_path: Path) -> None:
    path = tmp_path / "prep.json"
    report = {
        "sample": "CASE",
        "status": "complete",
        "started_at": "2026-01-01T00:00:00Z",
        "completed_at": "2026-01-01T00:00:10Z",
        "input_pairs": 2,
        "input_sequenced_bases": 400,
        "reuse_identity": {"inputs": [{"size": 1024**2}, {"size": 2 * 1024**2}]},
        "outputs": [{"size": 512}, {"size": 512}],
    }
    row = summary.preparation_rows([path], [report])[0]
    assert row["wall_seconds"] == 10
    assert row["mean_input_read_length"] == 100
    assert row["io_in_mb"] == 3


def test_elapsed_is_not_summed_parallel_wall_time() -> None:
    batch = {
        "launches": [
            {
                "prepared_at": "2026-01-01T00:00:00Z",
                "partial_at": "2026-01-01T00:01:40Z",
                "targets": ["coverage"],
            },
            {
                "prepared_at": "2026-01-01T00:02:00Z",
                "complete_at": "2026-01-01T00:05:00Z",
                "default_dag": True,
            },
        ]
    }
    values = summary.elapsed_metrics(batch)
    assert values["experienced_controller_elapsed_seconds"] == 300
    assert values["calibration_elapsed_seconds"] == 100
    assert values["full_dag_elapsed_seconds"] == 180


def test_scheduler_wait_is_only_derived_from_valid_timestamps() -> None:
    rows = [
        {"submit": "2026-01-01T00:00:00Z", "start": "2026-01-01T00:00:30Z"},
        {"submit": "unknown", "start": "2026-01-01T00:00:40Z"},
    ]
    assert summary.scheduler_wait_seconds(rows) == 30
    assert summary.scheduler_wait_seconds([]) is None


def test_optional_scheduler_failure_does_not_abort(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail(_job_ids: list[str]) -> list[dict[str, str]]:
        raise summary.ReportError("accounting unavailable")

    monkeypatch.setattr(summary, "scheduler_accounting", fail)
    rows, error = summary.optional_scheduler_accounting(True, ["123"])
    assert rows == []
    assert error == "accounting unavailable"


def test_complete_report_with_missing_metrics(tmp_path: Path) -> None:
    batch = tmp_path / "runs/benchmark"
    (batch / "config/current").mkdir(parents=True)
    (batch / "results/qc/coverage").mkdir(parents=True)
    (batch / "batch.json").write_text(
        json.dumps(
            {
                "batch": "benchmark",
                "launches": [
                    {
                        "prepared_at": "2026-01-01T00:00:00Z",
                        "complete_at": "2026-01-01T00:01:00Z",
                        "default_dag": True,
                        "git": {"revision": "abc", "dirty": False},
                    }
                ],
            }
        )
    )
    (batch / "config/current/config.yaml").write_text(
        yaml.safe_dump(
            {
                "reference": {"build": "grch38"},
                "analysis": {"wgs": {"max_concurrent_mutect2_shards": 8}},
            }
        )
    )
    (batch / "results/analysis_manifest.json").write_text(
        json.dumps({"reference": {"sha256": "xyz"}})
    )
    (batch / "results/qc/coverage/CASE.wgs_coverage_mqc.tsv").write_text(
        "# custom\nsample\tmean_autosomal_depth\nCASE\t6.1\n"
    )
    write_benchmark(
        batch / "results/benchmarks/qc/CASE.samtools.tsv", ["seconds", "max_memory_mb"], [[7, 50]]
    )
    prep_paths = []
    for sample in ("CASE", "NORMAL"):
        path = tmp_path / f"{sample}.json"
        path.write_text(json.dumps({"reuse_identity": {"callable_bases": 100}, "sample": sample}))
        prep_paths.append(path)
    history = tmp_path / "calibration.json"
    history.write_text(json.dumps({"attempts": [{"sample": "CASE", "depth": 6.1}]}))
    prefix = tmp_path / "report/performance"
    argv = [
        "--batch",
        str(batch),
        "--preparation-report",
        str(prep_paths[0]),
        "--preparation-report",
        str(prep_paths[1]),
        "--calibration-history",
        str(history),
        "--output-prefix",
        str(prefix),
    ]
    assert summary.main(argv) == 0
    payload = json.loads(prefix.with_suffix(".json").read_text())
    assert payload["timing"]["experienced_controller_elapsed_seconds"] == 60
    assert payload["timing"]["summed_job_wall_seconds"] == 7
    assert payload["timing"]["summed_cpu_seconds"] is None
    assert payload["achieved_depths"] == {"CASE": 6.1}
    assert payload["scheduler_accounting"] == []
    assert "Cluster, queue, filesystem" in prefix.with_suffix(".md").read_text()
