#!/usr/bin/env python3
"""Annotation configuration validation and deterministic CSQ helpers."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any


REQUIRED_RESOURCES = (
    "cache", "plugins", "clinvar", "gnomad", "revel", "cadd_snv", "cadd_indel",
    "alphamissense", "spliceai_snv", "spliceai_indel",
)
REQUIRED_METADATA = ("build", "version", "path", "expected_sha256", "source", "licence", "access_date")
DEFAULT_PICK_ORDER = (
    "mane_select", "mane_plus_clinical", "canonical", "appris", "tsl",
    "biotype", "ccds", "rank", "length",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checksum_path(path: Path) -> str:
    if path.is_file():
        return sha256(path)
    digest = hashlib.sha256()
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        digest.update(child.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        digest.update(sha256(child).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def validate_annotation_config(config: dict[str, Any], reference_build: str) -> list[dict[str, Any]]:
    annotation = config.get("annotation", {})
    if not annotation.get("enabled", False):
        return []
    if annotation.get("version") != 116 or annotation.get("cache_version") != 116 or annotation.get("plugin_release") != 116:
        raise ValueError("annotation requires VEP, cache, and plugin release 116")
    image = annotation.get("docker_image", "")
    if not re.fullmatch(r"ensemblorg/ensembl-vep@sha256:[0-9a-f]{64}", str(image)):
        raise ValueError("annotation.docker_image must pin the official Ensembl VEP image by full digest")
    pick_order = annotation.get("pick_order")
    if not isinstance(pick_order, list) or not pick_order or len(pick_order) != len(set(pick_order)):
        raise ValueError("annotation.pick_order must be a non-empty list without duplicates")
    resources = annotation.get("resources")
    if not isinstance(resources, dict):
        raise ValueError("annotation.resources must be a mapping")
    missing = [name for name in REQUIRED_RESOURCES if name not in resources]
    if missing:
        raise ValueError("annotation resources missing: " + ", ".join(missing))
    validated = []
    for name in sorted(resources):
        resource = resources[name]
        if not isinstance(resource, dict):
            raise ValueError(f"annotation resource {name!r} must be a mapping")
        absent = [key for key in REQUIRED_METADATA if not resource.get(key)]
        if absent:
            raise ValueError(f"annotation resource {name!r} missing metadata: {', '.join(absent)}")
        if str(resource["build"]).lower() != reference_build:
            raise ValueError(f"annotation resource {name!r} build does not match {reference_build}")
        path = Path(resource["path"])
        if not path.exists():
            raise ValueError(f"annotation resource {name!r} does not exist: {path}")
        checksum = checksum_path(path)
        if checksum != str(resource["expected_sha256"]).lower():
            raise ValueError(f"annotation resource {name!r} checksum mismatch")
        index = resource.get("index")
        if index is not None:
            if not isinstance(index, dict) or not index.get("path") or not index.get("expected_sha256"):
                raise ValueError(f"annotation resource {name!r} index requires path and expected_sha256")
            index_path = Path(index["path"])
            if not index_path.is_file():
                raise ValueError(f"annotation resource {name!r} index does not exist: {index_path}")
            if sha256(index_path) != str(index["expected_sha256"]).lower():
                raise ValueError(f"annotation resource {name!r} index checksum mismatch")
            if path.is_file() and index_path.stat().st_mtime_ns < path.stat().st_mtime_ns:
                raise ValueError(f"annotation resource {name!r} index is stale")
        size = path.stat().st_size if path.is_file() else sum(child.stat().st_size for child in path.rglob("*") if child.is_file())
        validated.append({"name": name, **resource, "size": size, "checksum": checksum, "validation_status": "valid"})
    return validated


def parse_csq_fields(header_line: str) -> list[str]:
    match = re.search(r"Format:\s*([^\"]+)", header_line)
    if not match:
        raise ValueError("CSQ header lacks a Format declaration")
    fields = [field.strip().rstrip(">") for field in match.group(1).split("|")]
    if not fields or any(not field for field in fields):
        raise ValueError("invalid CSQ Format declaration")
    return fields


def parse_csq(value: str | None, fields: list[str]) -> list[dict[str, str]]:
    if not value or value == ".":
        return []
    consequences = []
    for item in value.split(","):
        values = item.split("|")
        values.extend([""] * (len(fields) - len(values)))
        consequences.append(dict(zip(fields, values[:len(fields)])))
    return consequences


def select_pick(consequences: list[dict[str, str]], pick_order: list[str] | tuple[str, ...] = DEFAULT_PICK_ORDER) -> dict[str, str]:
    if not consequences:
        return {}
    flagged = [item for item in consequences if item.get("PICK") == "1"]
    if len(flagged) == 1:
        return flagged[0]
    key_map = {
        "mane_select": "MANE_SELECT", "mane_plus_clinical": "MANE_PLUS_CLINICAL",
        "canonical": "CANONICAL", "appris": "APPRIS", "tsl": "TSL",
        "biotype": "BIOTYPE", "ccds": "CCDS", "rank": "EXON", "length": "Feature",
    }
    def rank(item: dict[str, str]) -> tuple[Any, ...]:
        values = []
        for name in pick_order:
            value = item.get(key_map.get(name, name.upper()), "")
            values.append((0 if value else 1, value))
        values.extend((item.get("Gene", ""), item.get("Feature", ""), item.get("Consequence", "")))
        return tuple(values)
    return min(consequences, key=rank)
