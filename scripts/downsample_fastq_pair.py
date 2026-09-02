#!/usr/bin/env python3
"""Plan or execute deterministic, pair-preserving FASTQ downsampling."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO, Callable, Iterator, TextIO

VERSION = "1.1.0"
GIB = 1024**3
ALGORITHM = {
    "version": 1,
    "digest": "blake2b",
    "digest_bytes": 8,
    "seed_encoding": "UTF-8 decimal followed by NUL",
    "canonical_name": "first header token without leading @ or terminal /1 or /2",
    "threshold": "unsigned big-endian digest < floor(fraction * 2^64)",
}


class FastqError(RuntimeError):
    """Raised for invalid inputs or unsafe execution conditions."""


@dataclass(frozen=True)
class Scan:
    pairs: int
    bases: int


@dataclass
class ProgressReporter:
    sample: str
    phase: str
    total_bytes: int
    interval_seconds: float = 60.0
    stream: TextIO = sys.stderr
    last_reported: float = field(default_factory=time.monotonic)

    def start(self) -> None:
        print(
            f"progress sample={self.sample} phase={self.phase} pairs=0 "
            f"compressed_bytes=0/{self.total_bytes} percent=0.00",
            file=self.stream,
            flush=True,
        )

    def update(self, pairs: int, r1_bytes: int, r2_bytes: int, force: bool = False) -> None:
        now = time.monotonic()
        if not force and now - self.last_reported < self.interval_seconds:
            return
        consumed = min(self.total_bytes, r1_bytes + r2_bytes)
        percent = 100.0 * consumed / self.total_bytes if self.total_bytes else 100.0
        print(
            f"progress sample={self.sample} phase={self.phase} pairs={pairs} "
            f"compressed_bytes={consumed}/{self.total_bytes} percent={percent:.2f}",
            file=self.stream,
            flush=True,
        )
        self.last_reported = now

    def finish(self, pairs: int) -> None:
        self.update(pairs, self.total_bytes, 0, force=True)


def utcnow() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def identity(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }


def sequence_bytes(line: bytes) -> bytes:
    return line.rstrip(b"\r\n")


def canonical_name(header: bytes) -> bytes:
    if not header.startswith(b"@"):
        raise FastqError("FASTQ header does not start with @")
    token = header[1:].rstrip(b"\r\n").split(None, 1)[0]
    if not token:
        raise FastqError("FASTQ header has an empty read name")
    if token.endswith((b"/1", b"/2")):
        token = token[:-2]
    return token


def records(handle: BinaryIO, label: str) -> Iterator[tuple[bytes, bytes, int]]:
    record_number = 0
    while True:
        header = handle.readline()
        if not header:
            return
        record_number += 1
        sequence = handle.readline()
        separator = handle.readline()
        quality = handle.readline()
        if not sequence or not separator or not quality:
            raise FastqError(f"{label} is truncated at record {record_number}")
        if not header.startswith(b"@"):
            raise FastqError(f"{label} record {record_number} header does not start with @")
        if not separator.startswith(b"+"):
            raise FastqError(f"{label} record {record_number} separator does not start with +")
        sequence_value = sequence_bytes(sequence)
        quality_value = sequence_bytes(quality)
        if not sequence_value:
            raise FastqError(f"{label} record {record_number} has an empty sequence")
        if len(sequence_value) != len(quality_value):
            raise FastqError(
                f"{label} record {record_number} sequence/quality lengths differ "
                f"({len(sequence_value)} != {len(quality_value)})"
            )
        yield header + sequence + separator + quality, canonical_name(header), len(sequence_value)


def compressed_position(handle: BinaryIO) -> int:
    """Return the underlying compressed-stream position for a gzip handle."""
    raw = getattr(handle, "fileobj", None)
    return int(raw.tell()) if raw is not None else 0


def paired_records(
    r1: Path,
    r2: Path,
    progress: Callable[[int, int, int], None] | None = None,
) -> Iterator[tuple[bytes, bytes, bytes, int]]:
    try:
        with gzip.open(r1, "rb") as one, gzip.open(r2, "rb") as two:
            left = records(one, str(r1))
            right = records(two, str(r2))
            pair_number = 0
            while True:
                a = next(left, None)
                b = next(right, None)
                if a is None and b is None:
                    return
                pair_number += 1
                if a is None or b is None:
                    raise FastqError(
                        f"FASTQ mates have unequal record counts at pair {pair_number}"
                    )
                record1, name1, bases1 = a
                record2, name2, bases2 = b
                if name1 != name2:
                    raise FastqError(
                        f"FASTQ mates are unsynchronized at pair {pair_number}: "
                        f"{name1!r} != {name2!r}"
                    )
                if progress is not None:
                    progress(pair_number, compressed_position(one), compressed_position(two))
                yield record1, record2, name1, bases1 + bases2
    except (gzip.BadGzipFile, EOFError, OSError) as exc:
        raise FastqError(f"could not read gzip FASTQs: {exc}") from exc


def scan_pair(r1: Path, r2: Path, progress: ProgressReporter | None = None) -> Scan:
    pairs = bases = 0
    if progress is not None:
        progress.start()
    callback = progress.update if progress is not None else None
    for _record1, _record2, _name, pair_bases in paired_records(r1, r2, callback):
        pairs += 1
        bases += pair_bases
    if pairs == 0:
        raise FastqError("FASTQ pair contains no records")
    if progress is not None:
        progress.finish(pairs)
    return Scan(pairs=pairs, bases=bases)


def territory_bases(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    total = 0
    try:
        with opener(path, "rt") as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith(("#", "@")):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    raise FastqError(
                        f"malformed territory row {number}: expected at least 3 columns"
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise FastqError(f"malformed territory coordinates at row {number}") from exc
                # Picard interval lists are 1-based inclusive; BED is 0-based half-open.
                length = (
                    end
                    - start
                    + (1 if path.name.endswith((".interval_list", ".interval_list.gz")) else 0)
                )
                if length <= 0:
                    raise FastqError(f"non-positive territory interval at row {number}")
                total += length
    except (gzip.BadGzipFile, EOFError, OSError) as exc:
        raise FastqError(f"could not read territory {path}: {exc}") from exc
    if total <= 0:
        raise FastqError(f"territory contains no callable bases: {path}")
    return total


def selected(name: bytes, seed: int, fraction: float) -> bool:
    threshold = math.floor(fraction * 2**64)
    digest = hashlib.blake2b(str(seed).encode() + b"\0" + name, digest_size=8).digest()
    return int.from_bytes(digest, "big") < threshold


def calibrated_fraction(
    old_fraction: float, observed_depth: float, target_depth: float = 6.0
) -> float:
    if not 0 < old_fraction <= 1 or observed_depth <= 0 or target_depth <= 0:
        raise FastqError(
            "calibration fractions and depths must be positive, with fraction no greater than 1"
        )
    value = old_fraction * target_depth / observed_depth
    if value > 1:
        raise FastqError(f"calibration requires fraction {value:.8g}, which exceeds 1")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def capacity(output1: Path, output2: Path, estimate: int) -> tuple[int, int]:
    parents = set()
    for output in (output1, output2):
        parent = output.parent.resolve()
        while not parent.exists():
            parent = parent.parent
        parents.add(parent)
    available = min(shutil.disk_usage(parent).free for parent in parents)
    required = max(math.ceil(estimate * 1.5), estimate + GIB)
    return available, required


def matching_report(report: Path, expected: dict, output1: Path, output2: Path) -> dict | None:
    if not report.is_file() or not output1.is_file() or not output2.is_file():
        return None
    try:
        payload = json.loads(report.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if payload.get("status") != "complete" or payload.get("reuse_identity") != expected:
        return None
    outputs = payload.get("outputs", [])
    if len(outputs) != 2:
        return None
    for path, metadata in zip((output1, output2), outputs, strict=True):
        if metadata.get("path") != str(path.resolve()) or metadata.get("sha256") != sha256(path):
            return None
        try:
            with gzip.open(path, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
        except (gzip.BadGzipFile, EOFError, OSError):
            return None
    return payload


def validated_scan(
    report: Path,
    r1: Path,
    r2: Path,
    sample: str,
    seed: int,
    target_depth: float,
    territory: Path,
    callable_bases: int,
    fraction_override: float | None,
) -> Scan:
    """Reuse source totals only when a prior plan exactly matches this request."""
    try:
        payload = json.loads(report.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise FastqError(f"could not read validated plan {report}: {exc}") from exc
    expected = {
        "inputs": [identity(r1), identity(r2)],
        "sample": sample,
        "seed": seed,
        "target_depth": target_depth,
        "territory": identity(territory),
        "callable_bases": callable_bases,
        "selection_algorithm": ALGORITHM,
        "calibrated_fraction_override": fraction_override,
    }
    mismatches = [key for key, value in expected.items() if payload.get(key) != value]
    if mismatches:
        raise FastqError(
            f"validated plan does not match current inputs or parameters: {', '.join(mismatches)}"
        )
    pairs = payload.get("input_pairs")
    bases = payload.get("input_sequenced_bases")
    if not isinstance(pairs, int) or pairs <= 0 or not isinstance(bases, int) or bases <= 0:
        raise FastqError("validated plan has invalid source totals")
    calculated = target_depth * callable_bases / bases
    applied = fraction_override if fraction_override is not None else calculated
    if (
        payload.get("calculated_initial_fraction") != calculated
        or payload.get("applied_fraction") != applied
    ):
        raise FastqError("validated plan has inconsistent sampling fractions")
    return Scan(pairs=pairs, bases=bases)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def pigz_writer(path: Path, threads: int) -> subprocess.Popen:
    handle = path.open("wb")
    try:
        process = subprocess.Popen(
            ["pigz", "-n", "-p", str(threads), "-c"], stdin=subprocess.PIPE, stdout=handle
        )
    except OSError as exc:
        handle.close()
        raise FastqError(f"could not start pigz: {exc}") from exc
    process._output_handle = handle  # type: ignore[attr-defined]
    return process


def close_writer(process: subprocess.Popen) -> None:
    assert process.stdin is not None
    process.stdin.close()
    returncode = process.wait()
    process._output_handle.close()  # type: ignore[attr-defined]
    if returncode:
        raise FastqError(f"pigz exited with status {returncode}")


def execute(
    r1: Path,
    r2: Path,
    output1: Path,
    output2: Path,
    fraction: float,
    seed: int,
    threads: int,
    sample: str,
    progress_seconds: float,
) -> Scan:
    token = uuid.uuid4().hex
    temporary1 = output1.with_name(f".{output1.name}.{token}.tmp")
    temporary2 = output2.with_name(f".{output2.name}.{token}.tmp")
    writer1 = writer2 = None
    output_pairs = output_bases = processed_pairs = 0
    try:
        writer1 = pigz_writer(temporary1, threads // 2)
        writer2 = pigz_writer(temporary2, threads - threads // 2)
        assert writer1.stdin is not None and writer2.stdin is not None
        source_progress = ProgressReporter(
            sample,
            "selection",
            r1.stat().st_size + r2.stat().st_size,
            progress_seconds,
        )
        source_progress.start()
        for record1, record2, name, pair_bases in paired_records(r1, r2, source_progress.update):
            processed_pairs += 1
            if selected(name, seed, fraction):
                writer1.stdin.write(record1)
                writer2.stdin.write(record2)
                output_pairs += 1
                output_bases += pair_bases
        source_progress.finish(processed_pairs)
        close_writer(writer1)
        writer1 = None
        close_writer(writer2)
        writer2 = None
        validated = scan_pair(
            temporary1,
            temporary2,
            ProgressReporter(
                sample,
                "output_validation",
                temporary1.stat().st_size + temporary2.stat().st_size,
                progress_seconds,
            ),
        )
        if validated != Scan(output_pairs, output_bases):
            raise FastqError("temporary output validation totals do not match written totals")
        output1.parent.mkdir(parents=True, exist_ok=True)
        output2.parent.mkdir(parents=True, exist_ok=True)
        os.replace(temporary1, output1)
        os.replace(temporary2, output2)
        return validated
    finally:
        for process in (writer1, writer2):
            if process is not None:
                if process.stdin and not process.stdin.closed:
                    process.stdin.close()
                if process.poll() is None:
                    process.terminate()
                process.wait()
                if not process._output_handle.closed:  # type: ignore[attr-defined]
                    process._output_handle.close()  # type: ignore[attr-defined]
        temporary1.unlink(missing_ok=True)
        temporary2.unlink(missing_ok=True)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--r1", required=True, type=Path)
    value.add_argument("--r2", required=True, type=Path)
    value.add_argument("--output-r1", required=True, type=Path)
    value.add_argument("--output-r2", required=True, type=Path)
    value.add_argument("--sample", required=True)
    value.add_argument("--target-depth", required=True, type=float)
    value.add_argument("--territory", required=True, type=Path)
    value.add_argument("--seed", required=True, type=int)
    value.add_argument("--threads", type=int, default=2)
    value.add_argument("--report", required=True, type=Path)
    value.add_argument("--sampling-fraction", type=float)
    value.add_argument(
        "--validated-plan",
        type=Path,
        help="reuse source totals from an exactly matching prior plan during execution",
    )
    value.add_argument(
        "--progress-seconds",
        type=float,
        default=60.0,
        help="minimum interval between progress messages (default: 60)",
    )
    value.add_argument("--execute", action="store_true")
    return value


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    started = utcnow()
    for path, label in ((args.r1, "R1"), (args.r2, "R2"), (args.territory, "territory")):
        if not path.is_file() or path.stat().st_size == 0:
            raise FastqError(f"{label} is missing or empty: {path}")
    if args.output_r1.resolve() == args.output_r2.resolve():
        raise FastqError("R1 and R2 output paths must differ")
    if args.target_depth <= 0:
        raise FastqError("target depth must be positive")
    if args.threads < 2:
        raise FastqError("at least two threads are required")
    if args.progress_seconds < 0:
        raise FastqError("progress interval must not be negative")
    if args.validated_plan is not None and not args.execute:
        raise FastqError("--validated-plan requires --execute")
    sample = args.sample.strip()
    if not sample:
        raise FastqError("sample label must not be empty")
    callable_bases = territory_bases(args.territory)
    if args.validated_plan is None:
        scan = scan_pair(
            args.r1,
            args.r2,
            ProgressReporter(
                sample,
                "validation_count",
                args.r1.stat().st_size + args.r2.stat().st_size,
                args.progress_seconds,
            ),
        )
        validated_plan_reused = False
    else:
        scan = validated_scan(
            args.validated_plan,
            args.r1,
            args.r2,
            sample,
            args.seed,
            args.target_depth,
            args.territory,
            callable_bases,
            args.sampling_fraction,
        )
        validated_plan_reused = True
        print(
            f"progress sample={sample} phase=validation_count status=reused pairs={scan.pairs}",
            file=sys.stderr,
            flush=True,
        )
    initial_fraction = args.target_depth * callable_bases / scan.bases
    if initial_fraction > 1:
        raise FastqError(
            f"source cannot supply target depth: required fraction {initial_fraction:.8g} exceeds 1"
        )
    fraction = args.sampling_fraction if args.sampling_fraction is not None else initial_fraction
    if not 0 < fraction <= 1:
        raise FastqError("applied sampling fraction must be greater than 0 and no greater than 1")
    estimate = math.ceil((args.r1.stat().st_size + args.r2.stat().st_size) * fraction)
    available, required = capacity(args.output_r1, args.output_r2, estimate)
    reuse_identity = {
        "inputs": [identity(args.r1), identity(args.r2)],
        "sample": sample,
        "seed": args.seed,
        "target_depth": args.target_depth,
        "territory": identity(args.territory),
        "callable_bases": callable_bases,
        "selection_algorithm": ALGORITHM,
        "calculated_initial_fraction": initial_fraction,
        "calibrated_fraction_override": args.sampling_fraction,
        "applied_fraction": fraction,
    }
    plan = {
        **reuse_identity,
        "input_pairs": scan.pairs,
        "input_sequenced_bases": scan.bases,
        "estimated_compressed_output_bytes": estimate,
        "required_free_bytes": required,
        "available_free_bytes": available,
        "validated_plan_reused": validated_plan_reused,
        "mode": "execute" if args.execute else "plan",
    }
    if not args.execute:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0
    if available < required:
        raise FastqError(
            f"insufficient output capacity: {available} bytes available, {required} required"
        )
    args.output_r1.parent.mkdir(parents=True, exist_ok=True)
    args.output_r2.parent.mkdir(parents=True, exist_ok=True)
    existing = matching_report(args.report, reuse_identity, args.output_r1, args.output_r2)
    if existing is not None:
        existing["reused"] = True
        print(json.dumps(existing, indent=2, sort_keys=True))
        return 0
    output_scan = execute(
        args.r1,
        args.r2,
        args.output_r1,
        args.output_r2,
        fraction,
        args.seed,
        args.threads,
        sample,
        args.progress_seconds,
    )
    payload = {
        "schema_version": 1,
        "tool": "downsample_fastq_pair.py",
        "tool_version": VERSION,
        "command": sys.argv if argv is None else ["downsample_fastq_pair.py", *argv],
        "started_at": started,
        "completed_at": utcnow(),
        "status": "complete",
        "reused": False,
        "reuse_identity": reuse_identity,
        **plan,
        "output_pairs": output_scan.pairs,
        "output_sequenced_bases": output_scan.bases,
        "nominal_output_depth": output_scan.bases / callable_bases,
        "outputs": [
            {**identity(path), "sha256": sha256(path)} for path in (args.output_r1, args.output_r2)
        ],
    }
    atomic_json(args.report, payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FastqError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
