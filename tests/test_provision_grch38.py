"""Tests for explicit, restart-safe GRCh38 resource provisioning."""

from __future__ import annotations

import base64
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "provision_grch38", ROOT / "scripts" / "provision_grch38.py"
)
provision_grch38 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = provision_grch38
SPEC.loader.exec_module(provision_grch38)


def resource(payload: bytes) -> provision_grch38.Resource:
    checksum = base64.b64encode(hashlib.md5(payload).digest()).decode()  # noqa: S324
    return provision_grch38.Resource(
        "resource.dat", "example", "resource.dat", "1", len(payload), checksum
    )


class Response(io.BytesIO):
    def __init__(self, payload: bytes, status: int, content_range: str = ""):
        super().__init__(payload)
        self.status = status
        self.headers = {"Content-Range": content_range}


def test_plan_is_non_mutating(tmp_path, monkeypatch, capsys):
    payload = b"resource"
    monkeypatch.setattr(provision_grch38, "RESOURCES", (resource(payload),))
    destination = tmp_path / "not-created"
    provision_grch38.provision(destination, execute=False)
    assert not destination.exists()
    output = capsys.readouterr().out
    assert "Plan only" in output
    assert "population_vcf" in output


def test_capacity_rejection(tmp_path, monkeypatch):
    monkeypatch.setattr(
        provision_grch38.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=10),
    )
    with pytest.raises(ValueError, match="Insufficient capacity"):
        provision_grch38.check_capacity(tmp_path, required=100)


def test_checksum_rejection(tmp_path):
    output = tmp_path / "resource.dat"
    output.write_bytes(b"wrong")
    expected = resource(b"right")
    expected = provision_grch38.Resource(
        expected.filename,
        expected.bucket,
        expected.object_name,
        expected.generation,
        output.stat().st_size,
        expected.md5_base64,
    )
    with pytest.raises(ValueError, match="wrong checksum"):
        provision_grch38.validate_complete(output, expected)


def test_partial_download_resumes_and_publishes_atomically(tmp_path, monkeypatch):
    payload = b"0123456789"
    expected = resource(payload)
    output = tmp_path / expected.filename
    partial = tmp_path / f"{expected.filename}.partial"
    partial.write_bytes(payload[:4])

    def open_response(request):
        assert request.headers["Range"] == "bytes=4-"
        return Response(payload[4:], 206, "bytes 4-9/10")

    monkeypatch.setattr(provision_grch38.urllib.request, "urlopen", open_response)
    assert provision_grch38.download(expected, output) == "downloaded"
    assert output.read_bytes() == payload
    assert not partial.exists()


def test_complete_partial_is_recovered_without_network(tmp_path, monkeypatch):
    payload = b"complete"
    expected = resource(payload)
    output = tmp_path / expected.filename
    partial = tmp_path / f"{expected.filename}.partial"
    partial.write_bytes(payload)

    def unexpected_request(_request):
        raise AssertionError("a complete partial file must not be downloaded again")

    monkeypatch.setattr(provision_grch38.urllib.request, "urlopen", unexpected_request)
    assert provision_grch38.download(expected, output) == "recovered"
    assert output.read_bytes() == payload
    assert not partial.exists()


def test_execute_writes_manifest_and_overlay(tmp_path, monkeypatch):
    payload = b"resource"
    expected = resource(payload)
    monkeypatch.setattr(provision_grch38, "RESOURCES", (expected,))

    def fake_download(_resource, output):
        output.write_bytes(payload)
        return "downloaded"

    def fake_interval(_dictionary, output):
        output.write_text("fixture interval\n")

    monkeypatch.setattr(provision_grch38, "download", fake_download)
    monkeypatch.setattr(provision_grch38, "create_interval_list", fake_interval)
    monkeypatch.setattr(provision_grch38, "validate_resources", lambda _destination: None)
    destination = tmp_path / "references"
    provision_grch38.provision(destination, execute=True)

    manifest = json.loads((destination / "resource-manifest.json").read_text())
    assert manifest["resources"][0]["generation"] == "1"
    assert manifest["resources"][0]["status"] == "downloaded"
    overlay = destination / "organoid-pipeline.reference.yaml"
    assert str(destination.resolve()) in overlay.read_text()
