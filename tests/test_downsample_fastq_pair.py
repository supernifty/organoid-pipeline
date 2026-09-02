from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import downsample_fastq_pair as fastq


def write_fastq(path: Path, names: list[str], mate: int, newline: bytes = b"\n") -> None:
    with gzip.GzipFile(filename="", mode="wb", fileobj=path.open("wb"), mtime=0) as handle:
        for index, name in enumerate(names):
            sequence = (b"ACGT" if index % 2 == 0 else b"TGCA") + bytes(str(mate), "ascii")
            handle.write(
                b"@"
                + name.encode()
                + f"/{mate} metadata".encode()
                + newline
                + sequence
                + newline
                + b"+"
                + newline
                + b"I" * len(sequence)
                + newline
            )


def fixture(tmp_path: Path, count: int = 100) -> tuple[Path, Path, Path]:
    names = [f"read-{number}" for number in range(count)]
    r1, r2, territory = (
        tmp_path / "source_R1.fastq.gz",
        tmp_path / "source_R2.fastq.gz",
        tmp_path / "territory.bed",
    )
    write_fastq(r1, names, 1)
    write_fastq(r2, names, 2)
    territory.write_text("chr1\t0\t100\n")
    return r1, r2, territory


def arguments(tmp_path: Path, r1: Path, r2: Path, territory: Path) -> list[str]:
    return [
        "--r1",
        str(r1),
        "--r2",
        str(r2),
        "--output-r1",
        str(tmp_path / "out" / "R1.fastq.gz"),
        "--output-r2",
        str(tmp_path / "out" / "R2.fastq.gz"),
        "--sample",
        "CASE",
        "--target-depth",
        "1",
        "--territory",
        str(territory),
        "--seed",
        "42",
        "--threads",
        "2",
        "--report",
        str(tmp_path / "out" / "report.json"),
    ]


def test_scan_preserves_crlf_and_normalizes_mate_suffixes(tmp_path: Path) -> None:
    r1, r2, _ = fixture(tmp_path, 3)
    write_fastq(r1, ["a", "b"], 1, b"\r\n")
    write_fastq(r2, ["a", "b"], 2, b"\r\n")
    assert fastq.scan_pair(r1, r2) == fastq.Scan(pairs=2, bases=20)
    records = list(fastq.paired_records(r1, r2))
    assert records[0][0].endswith(b"\r\n")


def test_selection_is_deterministic_nested_and_seeded() -> None:
    names = [f"read-{number}".encode() for number in range(1000)]
    low = {name for name in names if fastq.selected(name, 7, 0.2)}
    high = {name for name in names if fastq.selected(name, 7, 0.6)}
    other = {name for name in names if fastq.selected(name, 8, 0.2)}
    assert low < high
    assert low != other
    assert all(fastq.selected(name, 7, 1.0) for name in names)


def test_calibration_arithmetic_and_fraction_gate() -> None:
    assert fastq.calibrated_fraction(0.4, 4.0) == pytest.approx(0.6)
    with pytest.raises(fastq.FastqError, match="exceeds 1"):
        fastq.calibrated_fraction(0.9, 4.0)


@pytest.mark.parametrize("fault", ["mismatch", "unequal", "truncated", "corrupt"])
def test_invalid_fastqs_fail(tmp_path: Path, fault: str) -> None:
    r1, r2, _ = fixture(tmp_path, 3)
    if fault == "mismatch":
        write_fastq(r2, ["read-0", "wrong", "read-2"], 2)
    elif fault == "unequal":
        write_fastq(r2, ["read-0", "read-1"], 2)
    elif fault == "truncated":
        with gzip.open(r2, "ab") as handle:
            handle.write(b"@partial\nAC\n")
    else:
        r2.write_bytes(b"not gzip")
    with pytest.raises(fastq.FastqError):
        fastq.scan_pair(r1, r2)


def test_plan_writes_nothing_and_reports_capacity(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    r1, r2, territory = fixture(tmp_path)
    argv = arguments(tmp_path, r1, r2, territory)
    assert fastq.main(argv) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "plan"
    assert payload["input_pairs"] == 100
    assert (
        payload["required_free_bytes"] >= payload["estimated_compressed_output_bytes"] + fastq.GIB
    )
    assert not (tmp_path / "out").exists()


def test_insufficient_source_and_capacity_fail(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    r1, r2, territory = fixture(tmp_path, 1)
    territory.write_text("chr1\t0\t100000\n")
    with pytest.raises(fastq.FastqError, match="cannot supply"):
        fastq.main(arguments(tmp_path, r1, r2, territory))
    territory.write_text("chr1\t0\t1\n")
    monkeypatch.setattr(fastq, "capacity", lambda *_: (0, 1))
    with pytest.raises(fastq.FastqError, match="insufficient output capacity"):
        fastq.main([*arguments(tmp_path, r1, r2, territory), "--execute"])


def test_execute_atomic_reproducible_and_restart_safe(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    r1, r2, territory = fixture(tmp_path)
    argv = [*arguments(tmp_path, r1, r2, territory), "--sampling-fraction", "0.4", "--execute"]
    monkeypatch.setattr(fastq, "GIB", 0)
    assert fastq.main(argv) == 0
    report_path = tmp_path / "out" / "report.json"
    report = json.loads(report_path.read_text())
    assert report["status"] == "complete"
    assert (
        report["output_pairs"]
        == fastq.scan_pair(tmp_path / "out/R1.fastq.gz", tmp_path / "out/R2.fastq.gz").pairs
    )
    retained = {
        name: (record1, record2)
        for record1, record2, name, _bases in fastq.paired_records(
            tmp_path / "out/R1.fastq.gz", tmp_path / "out/R2.fastq.gz"
        )
    }
    sources = {
        name: (record1, record2) for record1, record2, name, _bases in fastq.paired_records(r1, r2)
    }
    assert retained
    assert all(records == sources[name] for name, records in retained.items())
    before = [(tmp_path / f"out/R{mate}.fastq.gz").read_bytes() for mate in (1, 2)]
    assert fastq.main(argv) == 0
    assert before == [(tmp_path / f"out/R{mate}.fastq.gz").read_bytes() for mate in (1, 2)]
    assert json.loads(report_path.read_text())["reused"] is False

    changed = argv.copy()
    changed[changed.index("0.4")] = "0.2"
    assert fastq.main(changed) == 0
    assert json.loads(report_path.read_text())["reuse_identity"]["applied_fraction"] == 0.2
    assert not list((tmp_path / "out").glob(".*.tmp"))


def test_failure_does_not_publish_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    r1, r2, territory = fixture(tmp_path)
    monkeypatch.setattr(fastq, "GIB", 0)

    def fail(*_args, **_kwargs):
        raise fastq.FastqError("injected")

    monkeypatch.setattr(fastq, "pigz_writer", fail)
    with pytest.raises(fastq.FastqError, match="injected"):
        fastq.main([*arguments(tmp_path, r1, r2, territory), "--execute"])
    assert not (tmp_path / "out/R1.fastq.gz").exists()
    assert not (tmp_path / "out/R2.fastq.gz").exists()
    assert not (tmp_path / "out/report.json").exists()
