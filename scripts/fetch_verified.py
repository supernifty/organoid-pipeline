#!/usr/bin/env python3
"""Explicitly download one reviewed resource after size/capacity checks."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import urllib.request
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if len(args.sha256) != 64:
        raise ValueError("A reviewed SHA-256 checksum is required")
    request = urllib.request.Request(args.url, method="HEAD")
    with urllib.request.urlopen(request) as response:
        size = int(response.headers.get("Content-Length", 0))
    destination = args.output.resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(destination.parent).free
    print(f"Remote size: {size / 2**30:.2f} GiB; available: {free / 2**30:.2f} GiB")
    if size and size * 1.2 > free:
        raise ValueError("Insufficient space with a 1.2× download margin")
    if not args.execute:
        print("Plan only; add --execute after reviewing URL, size, destination, and checksum.")
        return
    temporary = destination.with_suffix(destination.suffix + ".partial")
    try:
        with urllib.request.urlopen(args.url) as source, temporary.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        observed = sha256(temporary)
        if observed.lower() != args.sha256.lower():
            raise ValueError(f"SHA-256 mismatch: expected {args.sha256}, observed {observed}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
