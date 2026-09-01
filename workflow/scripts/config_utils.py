#!/usr/bin/env python3
"""Helpers for loading the tracked config with an optional local overlay."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


def deep_merge(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    """Return a recursive merge matching Snakemake's configfile semantics."""
    merged = deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = deep_merge(dict(merged[key]), value)
        else:
            merged[key] = deepcopy(value)
    return merged


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"Configuration file must contain a mapping: {path}")
    return value


def load_merged_config(base: Path, overlay: Path | None = None) -> dict[str, Any]:
    config = read_yaml(base)
    if overlay and overlay.exists():
        config = deep_merge(config, read_yaml(overlay))
    return config
