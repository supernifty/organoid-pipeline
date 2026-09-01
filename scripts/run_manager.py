#!/usr/bin/env python3
"""Persistent, concurrency-safe run management for the somatic pipeline."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any, Iterable

import yaml


SCHEMA_VERSION = 1
BATCH_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
ACTIVE_STATES = {"submitted", "running"}
CATEGORIES = (
    "alignment", "strelka", "mutect2", "germline_calling", "germline_filtering",
    "filtering", "signatures", "annotation", "coverage", "somalier", "qc",
)
PATH_KEYS = {
    "genome", "genome_dict", "regions", "regions_metadata", "wgs_calling_regions",
    "wgs_exclude_regions", "gnomad", "contamination_sites", "panel_of_normals",
    "image", "singularity_image", "sbs_definition", "id_definition", "dbs_definition",
    "resource", "pixi_bin", "path", "sites_vcf", "sites_vcf_index", "labels",
    "reference_somalier_dir",
}
SAMPLE_PATH_KEYS = {"fastq_1", "fastq_2", "cram", "crai", "final_vcf"}


class RunError(RuntimeError):
    pass


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def read_yaml(path: Path) -> dict[str, Any]:
    try:
        value = yaml.safe_load(path.read_text()) or {}
    except OSError as exc:
        raise RunError(f"cannot read YAML {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RunError(f"YAML root must be a mapping: {path}")
    return value


def deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_yaml(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w") as handle:
        yaml.safe_dump(value, handle, sort_keys=True)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def file_identity(path_value: str | None) -> Any:
    if not path_value:
        return None
    path = Path(path_value)
    value: dict[str, Any] = {"path": str(path)}
    try:
        stat = path.stat()
        value.update(size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    except OSError:
        value["missing"] = True
    return value


def resolve_paths(value: Any, base: Path, keys: set[str], parent_key: str = "") -> Any:
    if isinstance(value, dict):
        return {key: resolve_paths(item, base, keys, key) for key, item in value.items()}
    if isinstance(value, list):
        return [resolve_paths(item, base, keys, parent_key) for item in value]
    if parent_key in keys and isinstance(value, str) and value:
        path = Path(value).expanduser()
        return str(path.resolve() if path.is_absolute() else (base / path).resolve())
    return value


def repository_root() -> Path:
    return Path(__file__).resolve().parent.parent


def run_dir(batch: str) -> Path:
    validate_batch(batch)
    return repository_root() / "runs" / batch


def validate_batch(batch: str) -> None:
    if not BATCH_RE.fullmatch(batch):
        raise RunError("batch must match [A-Za-z0-9][A-Za-z0-9._-]*")


def load_record(batch: str) -> dict[str, Any]:
    path = run_dir(batch) / "batch.json"
    if not path.is_file():
        raise RunError(f"batch does not exist: {batch}")
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RunError(f"invalid batch record {path}: {exc}") from exc


def save_record(batch: str, record: dict[str, Any]) -> None:
    record["updated_at"] = utcnow()
    atomic_json(run_dir(batch) / "batch.json", record)


def git_provenance(root: Path, history: Path) -> dict[str, Any]:
    def git(*args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(["git", *args], cwd=root, text=True, capture_output=True, check=False)

    head = git("rev-parse", "HEAD")
    status = git("status", "--porcelain")
    diff = git("diff", "--binary", "HEAD", "--")
    diff_path = history / "tracked.diff"
    diff_path.write_text(diff.stdout if diff.returncode == 0 else "")
    return {
        "revision": head.stdout.strip() if head.returncode == 0 else None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
        "tracked_diff": str(diff_path),
    }


def protected_identity(config: dict[str, Any]) -> dict[str, Any]:
    analysis = config.get("analysis", {})
    reference = config.get("reference", {})
    analysis_type = analysis.get("type", "wes")
    territory_keys = ["regions"] if analysis_type == "wes" else ["wgs_calling_regions", "wgs_exclude_regions"]
    return {
        "analysis_type": analysis_type,
        "reference": {"build": reference.get("build"), "genome": file_identity(reference.get("genome"))},
        "callable_territory": {key: file_identity(reference.get(key)) for key in territory_keys},
    }


def category_fingerprints(config: dict[str, Any]) -> dict[str, str]:
    reference = config.get("reference", {})
    analysis = config.get("analysis", {})
    territory = {key: reference.get(key) for key in ("regions", "wgs_calling_regions", "wgs_exclude_regions")}
    values = {
        "alignment": {
            "cram_version": config.get("cram_version"), "trimmomatic": config.get("trimmomatic"),
            "skip_trimming": config.get("skip_trimming"), "picard": config.get("picard"),
            "reference": {key: reference.get(key) for key in ("build", "genome", "bwa_index_suffix")},
        },
        "strelka": {"strelka": config.get("strelka"), "analysis_type": analysis.get("type"),
                    "reference": {"genome": reference.get("genome"), **territory}},
        "mutect2": {
            "mutect2": config.get("mutect2"), "gatk": config.get("gatk"), "analysis": analysis,
            "reference": {key: reference.get(key) for key in
                          ("build", "genome", "genome_dict", "gnomad", "contamination_sites", "panel_of_normals")}
                         | territory,
        },
        "germline_calling": {
            "enabled": config.get("germline", {}).get("enabled", True),
            "gatk": config.get("gatk"), "analysis": analysis,
            "reference": {key: reference.get(key) for key in ("build", "genome", "genome_dict")} | territory,
        },
        "germline_filtering": {
            "hard_filters": config.get("germline", {}).get("hard_filters"),
            "gatk": config.get("gatk"),
        },
        "filtering": {"filtering": config.get("filtering"), "territory": territory},
        "signatures": {"mutational_signatures": config.get("mutational_signatures"), "vt": config.get("vt"),
                       "reference": {"build": reference.get("build"), "genome": reference.get("genome")}},
        "annotation": {"annotation": config.get("annotation"), "reference": {"build": reference.get("build"), "genome": reference.get("genome")}},
        "coverage": {"coverage": config.get("coverage"), "analysis_type": analysis.get("type"),
                     "reference": {"build": reference.get("build"), "genome": reference.get("genome")}},
        "somalier": {"somalier": config.get("somalier"), "reference": {"build": reference.get("build"), "genome": reference.get("genome")}},
        "qc": {"qc": config.get("qc"), "hotspots": config.get("hotspots"), "analysis_type": analysis.get("type"),
               "reference": {"genome": reference.get("genome")}},
    }
    return {category: digest(values[category]) for category in CATEGORIES}


def sample_semantics(samples: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(samples)
    for values in result.get("samples", {}).values():
        if isinstance(values, dict):
            for key in SAMPLE_PATH_KEYS:
                if key in values:
                    values[key] = file_identity(values[key])
    return result


def normal_ids(samples: dict[str, Any]) -> list[str]:
    all_samples = set(samples.get("samples", {}))
    tumours = set(samples.get("tumours", {}))
    paired_normals = set(samples.get("tumours", {}).values())
    return sorted((all_samples - tumours) | paired_normals)


def pon_fingerprint(config: dict[str, Any], samples: dict[str, Any], identity: dict[str, Any]) -> str:
    sample_map = samples.get("samples", {})
    normals = {}
    for normal in normal_ids(samples):
        values = copy.deepcopy(sample_map.get(normal, {}))
        inputs = {key: file_identity(values.get(key)) for key in SAMPLE_PATH_KEYS if values.get(key)}
        normals[normal] = {"bam_sample": values.get("bam_sample", normal), "inputs": inputs}
    value = {
        "analysis": identity["analysis_type"],
        "reference": identity["reference"],
        "territory": identity["callable_territory"],
        "normals": normals,
        "alignment": {key: config.get(key) for key in ("cram_version", "trimmomatic", "skip_trimming")},
        "tools": {"mutect2": config.get("mutect2"), "gatk": config.get("gatk")},
    }
    return digest(value)


def cache_entry(fingerprint: str) -> Path:
    return repository_root() / "cache" / "pon" / fingerprint


def sha256_file(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def validate_cache(path: Path, fingerprint: str) -> bool:
    if not path.exists():
        return False
    manifest_path = path / "manifest.json"
    vcf = path / "mutect2.pon.vcf.gz"
    tbi = path / "mutect2.pon.vcf.gz.tbi"
    try:
        manifest = json.loads(manifest_path.read_text())
        valid = (
            manifest.get("fingerprint") == fingerprint
            and vcf.is_file() and tbi.is_file()
            and manifest.get("vcf", {}).get("size") == vcf.stat().st_size
            and manifest.get("vcf", {}).get("sha256") == sha256_file(vcf)
            and manifest.get("tbi", {}).get("size") == tbi.stat().st_size
            and manifest.get("tbi", {}).get("sha256") == sha256_file(tbi)
        )
    except (OSError, json.JSONDecodeError):
        valid = False
    if not valid:
        raise RunError(f"corrupt or mismatched PoN cache entry: {path}")
    return True


def publish_pon(batch: str, record: dict[str, Any]) -> None:
    fingerprint = record.get("pon_fingerprint")
    if not fingerprint:
        return
    source = run_dir(batch) / "results/variants/mutect2.pon.vcf.gz"
    source_tbi = Path(f"{source}.tbi")
    if not source.is_file() or not source_tbi.is_file():
        return
    parent = cache_entry(fingerprint).parent
    parent.mkdir(parents=True, exist_ok=True)
    lock_path = parent / f".{fingerprint}.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        destination = cache_entry(fingerprint)
        if destination.exists():
            validate_cache(destination, fingerprint)
            return
        required = source.stat().st_size + source_tbi.stat().st_size
        if shutil.disk_usage(parent).free < required * 2:
            raise RunError(f"insufficient capacity to publish {required} bytes to {parent}")
        staging = parent / f".{fingerprint}.{uuid.uuid4().hex}.staging"
        staging.mkdir()
        try:
            for original, name in ((source, "mutect2.pon.vcf.gz"), (source_tbi, "mutect2.pon.vcf.gz.tbi")):
                target = staging / name
                try:
                    os.link(original, target)
                except OSError:
                    shutil.copy2(original, target)
            manifest = {
                "fingerprint": fingerprint,
                "created_at": utcnow(),
                "vcf": {"size": (staging / "mutect2.pon.vcf.gz").stat().st_size,
                        "sha256": sha256_file(staging / "mutect2.pon.vcf.gz")},
                "tbi": {"size": (staging / "mutect2.pon.vcf.gz.tbi").stat().st_size,
                        "sha256": sha256_file(staging / "mutect2.pon.vcf.gz.tbi")},
            }
            atomic_json(staging / "manifest.json", manifest)
            validate_cache(staging, fingerprint)
            os.replace(staging, destination)
        finally:
            if staging.exists():
                shutil.rmtree(staging)


def unchanged_samples(old: dict[str, Any], new: dict[str, Any]) -> set[str]:
    old_values = sample_semantics(old).get("samples", {})
    new_values = sample_semantics(new).get("samples", {})
    return {name for name in old_values.keys() & new_values.keys() if old_values[name] == new_values[name]}


def reusable(relative: Path, unchanged: set[str], unchanged_pairs: set[str], same_normals: bool,
             old_fp: dict[str, str], new_fp: dict[str, str]) -> bool:
    parts = relative.parts
    if not parts or parts[0] == "aggregate" or relative.name.endswith(".done"):
        return False
    text = relative.as_posix()
    if parts[0] in {"bam", "cram", "trimmed"} or text.startswith("qc/fastqc/") or text.startswith("qc/metrics/"):
        category = "qc" if text.startswith("qc/") else "alignment"
        return old_fp.get("alignment") == new_fp.get("alignment") and (
            category != "qc" or old_fp.get("qc") == new_fp.get("qc")
        ) and any(name in relative.name or (len(parts) > 2 and parts[2] == name) for name in unchanged)
    if text.startswith("qc/coverage/"):
        return all(old_fp.get(key) == new_fp.get(key) for key in ("alignment", "coverage")) and any(name in relative.name for name in unchanged)
    if text.startswith("qc/somalier/extracted/"):
        return all(old_fp.get(key) == new_fp.get(key) for key in ("alignment", "somalier")) and any(name in relative.name for name in unchanged)
    if text.startswith("qc/somalier/"):
        return False
    if parts[0] == "variants":
        tumour = relative.name.split(".", 1)[0]
        if tumour not in unchanged_pairs:
            return False
        if ".strelka." in relative.name:
            return old_fp.get("strelka") == new_fp.get("strelka") and old_fp.get("alignment") == new_fp.get("alignment")
        if ".mutect2." in relative.name:
            return same_normals and old_fp.get("mutect2") == new_fp.get("mutect2") and old_fp.get("alignment") == new_fp.get("alignment")
        if ".intersect." in relative.name:
            return same_normals and all(old_fp.get(key) == new_fp.get(key) for key in ("alignment", "strelka", "mutect2", "filtering"))
        return False
    if parts[0] == "germline":
        normal = next(
            (name for name in sorted(unchanged, key=len, reverse=True) if relative.name.startswith(f"{name}.")),
            None,
        )
        categories = ["alignment", "germline_calling"]
        if ".filtered." in relative.name:
            categories.append("germline_filtering")
        return normal in unchanged and all(old_fp.get(key) == new_fp.get(key) for key in categories)
    if parts[0] == "annotations":
        sample = next((name for name in sorted(unchanged, key=len, reverse=True) if relative.name.startswith(f"{name}.")), None)
        upstream = ("alignment", "strelka", "mutect2", "filtering") if parts[1] == "somatic" else ("alignment", "germline_calling", "germline_filtering")
        return sample in unchanged and all(old_fp.get(key) == new_fp.get(key) for key in (*upstream, "annotation"))
    if parts[0] == "signatures":
        return same_normals and all(old_fp.get(key) == new_fp.get(key) for key in
                                    ("alignment", "strelka", "mutect2", "filtering", "signatures")) \
            and any(name in relative.name for name in unchanged_pairs)
    return False


def reconcile(directory: Path, old_samples: dict[str, Any], new_samples: dict[str, Any], old_fp: dict[str, str], new_fp: dict[str, str]) -> str | None:
    if digest(sample_semantics(old_samples)) == digest(sample_semantics(new_samples)) and old_fp == new_fp:
        return None
    results = directory / "results"
    history_root = directory / "history/results"
    history_root.mkdir(parents=True, exist_ok=True)
    revision = f"{len(list(history_root.iterdir())) + 1:04d}"
    archived = history_root / revision
    if results.exists():
        os.replace(results, archived)
    else:
        archived.mkdir()
    results.mkdir()
    unchanged = unchanged_samples(old_samples, new_samples)
    old_pairs = old_samples.get("tumours", {})
    new_pairs = new_samples.get("tumours", {})
    pairs = {tumour for tumour in old_pairs.keys() & new_pairs.keys()
             if old_pairs[tumour] == new_pairs[tumour] and tumour in unchanged and new_pairs[tumour] in unchanged}
    old_normals = normal_ids(old_samples)
    new_normals = normal_ids(new_samples)
    same_normals = old_normals == new_normals
    for source in archived.rglob("*"):
        if source.is_symlink():
            continue
        if source.is_file():
            relative = source.relative_to(archived)
            if reusable(relative, unchanged, pairs, same_normals, old_fp, new_fp):
                target = results / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                os.link(source, target)
    temporary = directory / "tmp"
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir()
    return revision


def prepare(args: argparse.Namespace) -> dict[str, Any]:
    root = repository_root()
    directory = run_dir(args.batch)
    exists = (directory / "batch.json").exists()
    if exists and not args.resume:
        raise RunError(f"batch already exists; use --resume: {args.batch}")
    if not exists and args.resume:
        raise RunError(f"cannot resume missing batch: {args.batch}")
    old_record = load_record(args.batch) if exists else None
    if old_record and old_record.get("state") in ACTIVE_STATES:
        raise RunError(f"batch is {old_record['state']}; recover it explicitly if stale")
    samples_source = Path(args.samples).expanduser().resolve() if args.samples else None
    if samples_source is None and old_record:
        samples_source = Path(old_record["samples_source"])
    if samples_source is None:
        raise RunError("--samples is required when creating a batch")
    if not samples_source.is_file():
        raise RunError(f"samples file does not exist: {samples_source}")

    base = read_yaml(root / "config/config.yaml")
    local = root / "config/config.local.yaml"
    effective = deep_merge(base, read_yaml(local)) if local.exists() else base
    effective = resolve_paths(effective, root, PATH_KEYS)
    samples = resolve_paths(read_yaml(samples_source), samples_source.parent, SAMPLE_PATH_KEYS)
    identity = protected_identity(effective)
    fingerprints = category_fingerprints(effective)
    if old_record and old_record.get("identity") != identity:
        raise RunError("analysis mode, reference identity, or callable territory changed within the batch")

    directory.mkdir(parents=True, exist_ok=True)
    for path in ("config/current", "config/history", "results", "tmp", "log", "history/results", ".snakemake"):
        (directory / path).mkdir(parents=True, exist_ok=True)
    launch_id = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    history = directory / "config/history" / launch_id
    history.mkdir()
    old_samples = read_yaml(directory / "config/current/samples.yaml") if old_record else {}
    revision = reconcile(directory, old_samples, samples, old_record.get("fingerprints", {}), fingerprints) if old_record else None

    pon_fp = pon_fingerprint(effective, samples, identity)
    entry = cache_entry(pon_fp)
    if entry.exists():
        validate_cache(entry, pon_fp)
        effective.setdefault("reference", {})["panel_of_normals"] = str(entry / "mutect2.pon.vcf.gz")
    effective["run_management"] = {
        "batch": args.batch,
        "config_file": str(directory / "config/current/config.yaml"),
        "samples_file": str(directory / "config/current/samples.yaml"),
        "repository": str(root),
        "fingerprints": fingerprints,
        "pon_fingerprint": pon_fp,
    }
    atomic_yaml(history / "config.yaml", effective)
    atomic_yaml(history / "samples.yaml", samples)
    atomic_yaml(directory / "config/current/config.yaml", effective)
    atomic_yaml(directory / "config/current/samples.yaml", samples)
    provenance = git_provenance(root, history)
    for immutable in (history / "config.yaml", history / "samples.yaml", history / "tracked.diff"):
        immutable.chmod(0o444)
    targets = list(args.target or [])
    launch = {
        "id": launch_id, "prepared_at": utcnow(), "command": args.command or [], "targets": targets,
        "default_dag": not targets, "dry_run": args.dry_run, "samples_hash": digest(sample_semantics(samples)),
        "config_hash": digest(effective), "fingerprints": fingerprints, "git": provenance,
        "state": "prepared", "scheduler_ids": [], "results_revision": revision,
    }
    atomic_json(history / "launch.json", launch)
    if old_record:
        record = old_record
    else:
        record = {"schema_version": SCHEMA_VERSION, "batch": args.batch, "created_at": utcnow(), "launches": []}
    record.update(state="created", samples_source=str(samples_source), identity=identity,
                  samples_hash=launch["samples_hash"], config_hash=launch["config_hash"],
                  fingerprints=fingerprints, pon_fingerprint=pon_fp, current_launch=launch_id,
                  targets=targets, last_error=None, active_controller=None)
    record.setdefault("launches", []).append(launch)
    save_record(args.batch, record)
    return {"batch": args.batch, "run_dir": str(directory), "launch": launch_id,
            "config": str(directory / "config/current/config.yaml"), "default_dag": not targets}


def locked_prepare(args: argparse.Namespace) -> dict[str, Any]:
    validate_batch(args.batch)
    lock_root = repository_root() / "runs"
    lock_root.mkdir(parents=True, exist_ok=True)
    with (lock_root / f".{args.batch}.prepare.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return prepare(args)


def find_launch(record: dict[str, Any], launch_id: str) -> dict[str, Any]:
    for launch in record.get("launches", []):
        if launch.get("id") == launch_id:
            return launch
    raise RunError(f"unknown launch {launch_id}")


def transition(batch: str, launch_id: str, state: str, job_id: str | None = None, error: str | None = None) -> None:
    record = load_record(batch)
    launch = find_launch(record, launch_id)
    state_is_late_submission = state == "submitted" and record.get("state") in {"running", "partial", "failed", "complete"}
    if not state_is_late_submission:
        record["state"] = state
        launch["state"] = state
        launch[f"{state}_at"] = utcnow()
    if job_id:
        launch.setdefault("scheduler_ids", []).append(job_id)
        record.setdefault("scheduler_ids", []).append(job_id)
    if error:
        record["last_error"] = error
        launch["error"] = error
    save_record(batch, record)
    atomic_json(run_dir(batch) / "config/history" / launch_id / "launch.json", launch)


def controller(args: argparse.Namespace) -> int:
    directory = run_dir(args.batch)
    lock_path = directory / "controller.lock"
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RunError("another controller is active for this batch") from exc
        record = load_record(args.batch)
        if record.get("current_launch") != args.launch:
            raise RunError("launch is no longer current")
        record["active_controller"] = {"pid": os.getpid(), "host": os.uname().nodename, "started_at": utcnow()}
        save_record(args.batch, record)
        transition(args.batch, args.launch, "running")
        try:
            completed = subprocess.run(args.command, cwd=directory, check=False)
        except OSError as exc:
            record = load_record(args.batch)
            record["active_controller"] = None
            save_record(args.batch, record)
            transition(args.batch, args.launch, "failed", error=f"could not start controller command: {exc}")
            raise RunError(f"could not start controller command: {exc}") from exc
        record = load_record(args.batch)
        record["active_controller"] = None
        save_record(args.batch, record)
        if completed.returncode == 0:
            try:
                publish_pon(args.batch, record)
                current_launch = find_launch(record, args.launch)
                final_state = "complete" if current_launch.get("default_dag") and not current_launch.get("dry_run") else "partial"
                transition(args.batch, args.launch, final_state)
            except Exception as exc:
                transition(args.batch, args.launch, "failed", error=str(exc))
                raise
        else:
            transition(args.batch, args.launch, "failed", error=f"controller command exited {completed.returncode}")
        return completed.returncode


def recover(batch: str) -> None:
    directory = run_dir(batch)
    record = load_record(batch)
    with (directory / "controller.lock").open("a+") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RunError("controller lock is still held; recovery refused") from exc
        if record.get("state") not in ACTIVE_STATES:
            raise RunError(f"batch state is not stale-active: {record.get('state')}")
        record.setdefault("recoveries", []).append({"at": utcnow(), "previous_state": record["state"]})
        record["active_controller"] = None
        record["last_error"] = "explicit stale-controller recovery"
        record["state"] = "failed"
        save_record(batch, record)


def reconcile_scheduler(batch: str, scheduler_state: str) -> None:
    """Reconcile an observed scheduler state without ever inferring DAG success."""
    record = load_record(batch)
    state = scheduler_state.strip().upper().split("+", 1)[0]
    record.setdefault("scheduler_observations", []).append({"at": utcnow(), "state": state})
    launch = find_launch(record, record["current_launch"])
    if state in {"PENDING", "CONFIGURING", "REQUEUED", "RESIZING"}:
        if record.get("state") == "created":
            record["state"] = "submitted"
            launch["state"] = "submitted"
    elif state in {"RUNNING", "COMPLETING", "STAGE_OUT"}:
        if record.get("state") not in {"complete", "partial", "failed"}:
            record["state"] = "running"
            launch["state"] = "running"
    elif state in {"COMPLETED"}:
        if record.get("state") not in {"complete", "partial"}:
            record["state"] = "failed"
            launch["state"] = "failed"
            record["last_error"] = "scheduler completed but controller did not finalize the DAG"
    elif state in {"BOOT_FAIL", "CANCELLED", "DEADLINE", "FAILED", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "REVOKED", "TIMEOUT"}:
        if record.get("state") not in {"complete", "partial"}:
            record["state"] = "failed"
            launch["state"] = "failed"
            record["last_error"] = f"scheduler reported {state}"
    else:
        record["last_error"] = f"unrecognized scheduler state observed: {state}"
    save_record(batch, record)
    atomic_json(run_dir(batch) / "config/history" / launch["id"] / "launch.json", launch)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--batch", required=True)
    prep.add_argument("--samples")
    prep.add_argument("--resume", action="store_true")
    prep.add_argument("--dry-run", action="store_true")
    prep.add_argument("--target", action="append")
    prep.add_argument("--command", nargs=argparse.REMAINDER)
    ctl = sub.add_parser("controller")
    ctl.add_argument("--batch", required=True)
    ctl.add_argument("--launch", required=True)
    ctl.add_argument("command", nargs=argparse.REMAINDER)
    trans = sub.add_parser("transition")
    trans.add_argument("--batch", required=True)
    trans.add_argument("--launch", required=True)
    trans.add_argument("--state", required=True, choices=["created", "submitted", "running", "partial", "failed", "complete"])
    trans.add_argument("--job-id")
    trans.add_argument("--error")
    status = sub.add_parser("status")
    status.add_argument("--batch", required=True)
    sub.add_parser("list")
    recovery = sub.add_parser("recover")
    recovery.add_argument("--batch", required=True)
    reconcile = sub.add_parser("reconcile-slurm")
    reconcile.add_argument("--batch", required=True)
    reconcile.add_argument("--state", required=True, help="state reported by squeue or sacct")
    publish = sub.add_parser("publish-pon")
    publish.add_argument("--batch", required=True)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.action == "prepare":
            print(json.dumps(locked_prepare(args), sort_keys=True))
        elif args.action == "controller":
            if not args.command:
                raise RunError("controller command is required")
            return controller(args)
        elif args.action == "transition":
            transition(args.batch, args.launch, args.state, args.job_id, args.error)
        elif args.action == "status":
            print(json.dumps(load_record(args.batch), sort_keys=True, indent=2))
        elif args.action == "list":
            root = repository_root() / "runs"
            records = []
            for path in sorted(root.glob("*/batch.json")) if root.exists() else []:
                try:
                    record = json.loads(path.read_text())
                    records.append({key: record.get(key) for key in ("batch", "state", "created_at", "updated_at", "current_launch")})
                except (OSError, json.JSONDecodeError):
                    records.append({"batch": path.parent.name, "state": "invalid"})
            print(json.dumps(records, sort_keys=True, indent=2))
        elif args.action == "recover":
            recover(args.batch)
        elif args.action == "reconcile-slurm":
            reconcile_scheduler(args.batch, args.state)
        elif args.action == "publish-pon":
            record = load_record(args.batch)
            publish_pon(args.batch, record)
        return 0
    except RunError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
