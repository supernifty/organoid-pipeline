from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import run_manager as manager


def write_yaml(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=True))


@pytest.fixture
def managed_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    root = tmp_path / "repo"
    (root / "config").mkdir(parents=True)
    genome = root / "reference/genome.fa"
    regions = root / "reference/regions.bed.gz"
    genome.parent.mkdir()
    genome.write_text(">1\nA\n")
    regions.write_bytes(b"regions")
    write_yaml(
        root / "config/config.yaml",
        {
            "analysis": {"type": "wes"},
            "reference": {"build": "grch37", "genome": "reference/genome.fa", "regions": "reference/regions.bed.gz", "panel_of_normals": None},
            "mutect2": {"create_panel_of_normals": True, "generated_panel_of_normals": "results/variants/mutect2.pon.vcf.gz"},
            "gatk": {"version": "4.4.0.0"},
            "germline": {
                "enabled": True,
                "max_concurrent_haplotypecaller_shards": 16,
                "hard_filters": {"snp": {"qd_min": 2.0}, "indel": {"qd_min": 2.0}},
            },
            "filtering": {"af_threshold": 0.1},
        },
    )
    fastq = root / "inputs/N.fastq.gz"
    fastq.parent.mkdir()
    fastq.write_bytes(b"reads")
    samples = root / "cohort.yaml"
    write_yaml(samples, {"samples": {"T": {"fastq_1": str(fastq)}, "N": {"fastq_1": str(fastq)}}, "tumours": {"T": "N"}})
    monkeypatch.setattr(manager, "repository_root", lambda: root)
    return root, samples


def prep(batch: str, samples: Path | None, resume: bool = False) -> dict:
    return manager.locked_prepare(argparse.Namespace(batch=batch, samples=str(samples) if samples else None, resume=resume, dry_run=False, target=[], command=["test"]))


def test_batch_name_and_create_resume_history(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    with pytest.raises(manager.RunError, match="must match"):
        prep("bad/name", samples)
    first = prep("batch-1", samples)
    assert Path(first["run_dir"]).is_dir()
    record = manager.load_record("batch-1")
    assert record["state"] == "created"
    assert record["samples_source"] == str(samples.resolve())
    with pytest.raises(manager.RunError, match="already exists"):
        prep("batch-1", samples)
    second = prep("batch-1", None, resume=True)
    assert second["launch"] != first["launch"]
    assert len(list((root / "runs/batch-1/config/history").iterdir())) == 2


def test_two_batches_have_fully_isolated_work_trees(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    prep("parallel-a", samples)
    prep("parallel-b", samples)
    for name in ("parallel-a", "parallel-b"):
        directory = root / "runs" / name
        for child in ("results", "tmp", "log", ".snakemake", "config/current", "config/history"):
            assert (directory / child).is_dir()
    assert (root / "runs/parallel-a/.snakemake").resolve() != (root / "runs/parallel-b/.snakemake").resolve()


def test_protected_identity_rejected_but_scientific_change_reconciles(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    prep("science", samples)
    result = root / "runs/science/results/variants/T.intersect.vcf.gz"
    result.parent.mkdir(parents=True, exist_ok=True)
    result.write_bytes(b"old")
    config = yaml.safe_load((root / "config/config.yaml").read_text())
    config["filtering"]["af_threshold"] = 0.2
    write_yaml(root / "config/config.yaml", config)
    prep("science", None, resume=True)
    assert not result.exists()
    assert (root / "runs/science/history/results/0001/variants/T.intersect.vcf.gz").is_file()
    config["reference"]["build"] = "grch38"
    write_yaml(root / "config/config.yaml", config)
    with pytest.raises(manager.RunError, match="reference identity"):
        prep("science", None, resume=True)


def test_sample_reconciliation_hardlinks_safe_outputs_only(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    prep("cohort", samples)
    results = root / "runs/cohort/results"
    cram = results / "cram/N.sorted.dups.cram"
    aggregate = results / "aggregate/qc_summary.html"
    cram.parent.mkdir(parents=True)
    aggregate.parent.mkdir(parents=True)
    cram.write_bytes(b"cram")
    aggregate.write_text("aggregate")
    values = yaml.safe_load(samples.read_text())
    values["samples"]["U"] = values["samples"]["N"].copy()
    write_yaml(samples, values)
    prep("cohort", None, resume=True)
    archived = root / "runs/cohort/history/results/0001"
    assert cram.is_file()
    assert os.stat(cram).st_ino == os.stat(archived / "cram/N.sorted.dups.cram").st_ino
    assert not aggregate.exists()
    assert (archived / "aggregate/qc_summary.html").is_file()


def test_pon_cache_publish_hit_and_corruption(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    prep("pon-a", samples)
    record = manager.load_record("pon-a")
    vcf = root / "runs/pon-a/results/variants/mutect2.pon.vcf.gz"
    vcf.parent.mkdir(parents=True)
    vcf.write_bytes(b"vcf")
    Path(f"{vcf}.tbi").write_bytes(b"index")
    manager.publish_pon("pon-a", record)
    entry = manager.cache_entry(record["pon_fingerprint"])
    assert manager.validate_cache(entry, record["pon_fingerprint"])
    second_samples = root / "second-cohort.yaml"
    values = yaml.safe_load(samples.read_text())
    values["samples"]["U"] = values["samples"].pop("T")
    values["tumours"] = {"U": "N"}
    write_yaml(second_samples, values)
    prep("pon-b", second_samples)
    assert manager.load_record("pon-b")["pon_fingerprint"] == record["pon_fingerprint"]
    effective = yaml.safe_load((root / "runs/pon-b/config/current/config.yaml").read_text())
    assert effective["reference"]["panel_of_normals"] == str(entry / "mutect2.pon.vcf.gz")
    (entry / "mutect2.pon.vcf.gz").write_bytes(b"corrupt")
    with pytest.raises(manager.RunError, match="corrupt"):
        manager.validate_cache(entry, record["pon_fingerprint"])


def test_germline_reuse_and_filter_invalidation(managed_repo: tuple[Path, Path]) -> None:
    root, samples = managed_repo
    prep("germline", samples)
    results = root / "runs/germline/results/germline"
    results.mkdir(parents=True)
    gvcf = results / "N.haplotypecaller.g.vcf.gz"
    gvcf_tbi = results / "N.haplotypecaller.g.vcf.gz.tbi"
    filtered = results / "N.haplotypecaller.filtered.vcf.gz"
    filtered_tbi = results / "N.haplotypecaller.filtered.vcf.gz.tbi"
    for path in (gvcf, gvcf_tbi, filtered, filtered_tbi):
        path.write_bytes(path.name.encode())

    config = yaml.safe_load((root / "config/config.yaml").read_text())
    config["germline"]["hard_filters"]["snp"]["qd_min"] = 3.0
    write_yaml(root / "config/config.yaml", config)
    prep("germline", None, resume=True)

    archived = root / "runs/germline/history/results/0001/germline"
    assert gvcf.is_file() and gvcf_tbi.is_file()
    assert os.stat(gvcf).st_ino == os.stat(archived / gvcf.name).st_ino
    assert not filtered.exists() and not filtered_tbi.exists()

    filtered.write_bytes(b"refiltered")
    filtered_tbi.write_bytes(b"refiltered-index")
    config["gatk"]["version"] = "4.5.0.0"
    write_yaml(root / "config/config.yaml", config)
    prep("germline", None, resume=True)
    assert not gvcf.exists() and not gvcf_tbi.exists()
    assert not filtered.exists() and not filtered_tbi.exists()


def test_new_qc_and_annotation_fingerprints_and_reuse_policy() -> None:
    base = {
        "analysis": {"type": "wes"}, "reference": {"build": "grch37", "genome": "ref.fa"},
        "annotation": {"enabled": True, "version": 116},
        "coverage": {"exon_enabled": True, "tumour_complete_coverage_depth": 20},
        "somalier": {"enabled": True, "minimum_depth": 20},
    }
    changed = yaml.safe_load(yaml.safe_dump(base))
    changed["coverage"]["tumour_complete_coverage_depth"] = 50
    old_fp, new_fp = manager.category_fingerprints(base), manager.category_fingerprints(changed)
    assert old_fp["coverage"] != new_fp["coverage"]
    assert old_fp["annotation"] == new_fp["annotation"]
    unchanged = {"T", "N"}
    assert manager.reusable(Path("qc/coverage/T.exon_coverage.tsv"), unchanged, {"T"}, True, old_fp, old_fp)
    assert manager.reusable(Path("qc/somalier/extracted/T.somalier"), unchanged, {"T"}, True, old_fp, old_fp)
    assert not manager.reusable(Path("qc/somalier/cohort.pairs.tsv"), unchanged, {"T"}, True, old_fp, old_fp)
    assert not manager.reusable(Path("aggregate/mutation_burden.tsv"), unchanged, {"T"}, True, old_fp, old_fp)
    assert manager.reusable(Path("annotations/somatic/T.intersect.annotated.vcf.gz"), unchanged, {"T"}, True, old_fp, old_fp)

def test_controller_states_and_target_never_completes(managed_repo: tuple[Path, Path]) -> None:
    _, samples = managed_repo
    prepared = manager.locked_prepare(argparse.Namespace(batch="target", samples=str(samples), resume=False, dry_run=False, target=["one.vcf"], command=["test"]))
    args = argparse.Namespace(batch="target", launch=prepared["launch"], command=["/bin/sh", "-c", "true"])
    assert manager.controller(args) == 0
    assert manager.load_record("target")["state"] == "partial"

    failed = manager.locked_prepare(argparse.Namespace(batch="failure", samples=str(samples), resume=False, dry_run=False, target=[], command=["test"]))
    args = argparse.Namespace(batch="failure", launch=failed["launch"], command=["/bin/sh", "-c", "exit 7"])
    assert manager.controller(args) == 7
    assert manager.load_record("failure")["state"] == "failed"


def test_stale_recovery_is_explicit(managed_repo: tuple[Path, Path]) -> None:
    _, samples = managed_repo
    prepared = prep("stale", samples)
    manager.transition("stale", prepared["launch"], "submitted", "123")
    with pytest.raises(manager.RunError, match="recover it explicitly"):
        prep("stale", None, resume=True)
    manager.recover("stale")
    assert manager.load_record("stale")["state"] == "failed"
    prep("stale", None, resume=True)


def test_mocked_scheduler_reconciliation_never_infers_completion(managed_repo: tuple[Path, Path]) -> None:
    _, samples = managed_repo
    prepared = prep("scheduler", samples)
    manager.reconcile_scheduler("scheduler", "PENDING")
    assert manager.load_record("scheduler")["state"] == "submitted"
    manager.reconcile_scheduler("scheduler", "RUNNING")
    assert manager.load_record("scheduler")["state"] == "running"
    manager.reconcile_scheduler("scheduler", "COMPLETED")
    record = manager.load_record("scheduler")
    assert record["state"] == "failed"
    assert "did not finalize" in record["last_error"]
