#!/usr/bin/env python3
"""Plan or provision the pinned Broad GRCh38 resources used by the pipeline."""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from analysis_mode import (  # noqa: E402
    PRIMARY_CONTIG_LENGTHS,
    open_text,
    read_fai,
    read_interval_list,
    read_sequence_dictionary,
    read_vcf_dictionary,
    validate_intervals,
)


@dataclass(frozen=True)
class Resource:
    filename: str
    bucket: str
    object_name: str
    generation: str
    size: int
    md5_base64: str

    @property
    def url(self) -> str:
        encoded = urllib.parse.quote(self.object_name, safe="/")
        return (
            f"https://storage.googleapis.com/{self.bucket}/{encoded}?generation={self.generation}"
        )


REFERENCE_BUCKET = "gcp-public-data--broad-references"
REFERENCE_PREFIX = "hg38/v0"
SOMATIC_BUCKET = "gatk-best-practices"
SOMATIC_PREFIX = "somatic-hg38"

RESOURCES = (
    Resource(
        "Homo_sapiens_assembly38.fasta",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta",
        "1575676516681666",
        3_249_912_778,
        "f/E0lT3MqMiZdFO7uAtrXg==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.fai",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.fai",
        "1575676516578189",
        160_928,
        "92NxsRNzSlbN4ja8A3LeCg==",
    ),
    Resource(
        "Homo_sapiens_assembly38.dict",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.dict",
        "1575676516486325",
        581_712,
        "OITGLrDlP6kkWe2b/xM65g==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.alt",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.alt",
        "1575676516489805",
        487_553,
        "sH5lqkQlvDZRQXVvXJgyjA==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.amb",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.amb",
        "1575676516504704",
        20_199,
        "5NxP23NYGY4IRxBlmVIKqQ==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.ann",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.ann",
        "1575676516518309",
        455_474,
        "r2Ee0LuUh/sbpKoafnrSHA==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.bwt",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.bwt",
        "1575676516703507",
        3_217_347_004,
        "fwyNz8hrfCzj46VBGNaPvQ==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.pac",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.pac",
        "1575676516584835",
        804_336_731,
        "F4hip5sEOi+XTvEOOHfvhg==",
    ),
    Resource(
        "Homo_sapiens_assembly38.fasta.64.sa",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly38.fasta.64.sa",
        "1575676516653931",
        1_608_673_512,
        "kaXV7TmG24p0eC5fRRnrXw==",
    ),
    Resource(
        "af-only-gnomad.hg38.vcf.gz",
        SOMATIC_BUCKET,
        f"{SOMATIC_PREFIX}/af-only-gnomad.hg38.vcf.gz",
        "1503507110526458",
        3_184_275_189,
        "pCCb5/tLWlqNO3eBMst0AQ==",
    ),
    Resource(
        "af-only-gnomad.hg38.vcf.gz.tbi",
        SOMATIC_BUCKET,
        f"{SOMATIC_PREFIX}/af-only-gnomad.hg38.vcf.gz.tbi",
        "1503506571246428",
        2_443_190,
        "p+/MsVGfBGwZzfnyhVnXRw==",
    ),
    Resource(
        "small_exac_common_3.hg38.vcf.gz",
        SOMATIC_BUCKET,
        f"{SOMATIC_PREFIX}/small_exac_common_3.hg38.vcf.gz",
        "1503507022893515",
        1_297_183,
        "THXBdVpFxk6K93hNt/3gCQ==",
    ),
    Resource(
        "small_exac_common_3.hg38.vcf.gz.tbi",
        SOMATIC_BUCKET,
        f"{SOMATIC_PREFIX}/small_exac_common_3.hg38.vcf.gz.tbi",
        "1503507007898573",
        242_095,
        "9lDR3aa9aMumXXfxMRR5hQ==",
    ),
)


def md5_base64(path: Path) -> str:
    digest = hashlib.md5()  # noqa: S324 - verifies the publisher's GCS checksum
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return base64.b64encode(digest.digest()).decode()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def human_size(size: int) -> str:
    return f"{size / 2**30:.2f} GiB"


def existing_parent(path: Path) -> Path:
    candidate = path.expanduser().absolute()
    while not candidate.exists():
        parent = candidate.parent
        if parent == candidate:
            raise ValueError(f"No existing parent for destination: {path}")
        candidate = parent
    return candidate


def validate_complete(path: Path, resource: Resource) -> bool:
    if not path.exists():
        return False
    if not path.is_file():
        raise ValueError(f"Resource destination is not a file: {path}")
    if path.stat().st_size != resource.size:
        raise ValueError(f"Existing resource has the wrong size: {path}; move it aside and rerun")
    if md5_base64(path) != resource.md5_base64:
        raise ValueError(
            f"Existing resource has the wrong checksum: {path}; move it aside and rerun"
        )
    return True


def remaining_bytes(destination: Path, resources=RESOURCES) -> int:
    remaining = 0
    for resource in resources:
        output = destination / resource.filename
        if output.exists():
            if validate_complete(output, resource):
                continue
        partial = output.with_name(f"{output.name}.partial")
        downloaded = partial.stat().st_size if partial.exists() else 0
        if downloaded > resource.size:
            raise ValueError(f"Partial download is larger than expected: {partial}")
        remaining += resource.size - downloaded
    return remaining


def check_capacity(destination: Path, required: int, margin: float = 1.2) -> int:
    available = shutil.disk_usage(existing_parent(destination)).free
    if required * margin > available:
        raise ValueError(
            f"Insufficient capacity: need {human_size(int(required * margin))} "
            f"including the {margin:.1f}x margin; {human_size(available)} available"
        )
    return available


def download(resource: Resource, output: Path) -> str:
    if validate_complete(output, resource):
        return "reused"
    partial = output.with_name(f"{output.name}.partial")
    offset = partial.stat().st_size if partial.exists() else 0
    if offset > resource.size:
        raise ValueError(f"Partial download is larger than expected: {partial}")
    if offset == resource.size:
        observed = md5_base64(partial)
        if observed != resource.md5_base64:
            raise ValueError(
                f"Checksum mismatch for {resource.filename}: "
                f"expected {resource.md5_base64}, observed {observed}"
            )
        os.replace(partial, output)
        return "recovered"

    request = urllib.request.Request(resource.url)
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    with urllib.request.urlopen(request) as response:  # noqa: S310 - pinned HTTPS URLs
        status = getattr(response, "status", None)
        append = offset > 0 and status == 206
        if append:
            content_range = response.headers.get("Content-Range", "")
            if not content_range.startswith(f"bytes {offset}-"):
                raise ValueError(f"Unexpected resume response for {resource.filename}")
        mode = "ab" if append else "wb"
        with partial.open(mode) as handle:
            downloaded = offset if append else 0
            started = time.monotonic()
            last_report = started
            while block := response.read(8 * 1024 * 1024):
                handle.write(block)
                downloaded += len(block)
                now = time.monotonic()
                if now - last_report >= 30 or downloaded == resource.size:
                    elapsed = max(now - started, 0.001)
                    transferred = downloaded - (offset if append else 0)
                    speed = transferred / elapsed / 2**20
                    percent = 100 * downloaded / resource.size
                    print(
                        f"progress: {resource.filename}: {human_size(downloaded)} / "
                        f"{human_size(resource.size)} ({percent:.1f}%), {speed:.1f} MiB/s",
                        flush=True,
                    )
                    last_report = now

    if partial.stat().st_size != resource.size:
        raise ValueError(
            f"Incomplete download for {resource.filename}: "
            f"{partial.stat().st_size} of {resource.size} bytes"
        )
    observed = md5_base64(partial)
    if observed != resource.md5_base64:
        raise ValueError(
            f"Checksum mismatch for {resource.filename}: "
            f"expected {resource.md5_base64}, observed {observed}"
        )
    os.replace(partial, output)
    return "downloaded"


def create_interval_list(dictionary: Path, output: Path) -> None:
    order, lengths = read_sequence_dictionary(dictionary)
    primary = PRIMARY_CONTIG_LENGTHS["grch38"]
    for contig, expected in primary.items():
        if lengths.get(contig) != expected:
            raise ValueError(
                f"Broad dictionary is incompatible at {contig}: "
                f"{lengths.get(contig)!r} != {expected}"
            )
    temporary = output.with_name(f"{output.name}.partial")
    with temporary.open("w") as handle:
        handle.write("@HD\tVN:1.6\tSO:coordinate\n")
        for contig in order:
            handle.write(f"@SQ\tSN:{contig}\tLN:{lengths[contig]}\n")
        for contig, length in primary.items():
            handle.write(f"{contig}\t1\t{length}\t+\tCALLABLE\n")
    os.replace(temporary, output)


def validate_resources(destination: Path) -> None:
    fai_order, fai_lengths = read_fai(destination / "Homo_sapiens_assembly38.fasta.fai")
    dictionary = destination / "Homo_sapiens_assembly38.dict"
    dictionary_order, dictionary_lengths = read_sequence_dictionary(dictionary)
    if (fai_order, fai_lengths) != (dictionary_order, dictionary_lengths):
        raise ValueError("Reference FAI and sequence dictionary do not match")
    primary = PRIMARY_CONTIG_LENGTHS["grch38"]
    if any(fai_lengths.get(contig) != length for contig, length in primary.items()):
        raise ValueError("Reference dictionary does not match canonical GRCh38")
    interval_dictionary, intervals = read_interval_list(
        destination / "wgs_calling_regions.interval_list"
    )
    if interval_dictionary != dictionary_lengths:
        raise ValueError("WGS interval-list dictionary does not match the reference")
    intervals = validate_intervals(intervals, fai_order, fai_lengths)
    if [contig for contig, _, _ in intervals] != list(primary):
        raise ValueError("WGS interval list does not contain each primary contig exactly once")
    for name in (
        "af-only-gnomad.hg38.vcf.gz",
        "small_exac_common_3.hg38.vcf.gz",
    ):
        vcf_dictionary = read_vcf_dictionary(destination / name)
        missing = [contig for contig in primary if contig not in vcf_dictionary]
        if missing:
            raise ValueError(f"{name} lacks primary contigs: {', '.join(missing)}")
    has_af = False
    with open_text(destination / "af-only-gnomad.hg38.vcf.gz") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=AF,"):
                has_af = True
            if line.startswith("#CHROM"):
                break
    if not has_af:
        raise ValueError("AF-only gnomAD does not declare the required INFO/AF field")


def config_text(destination: Path) -> str:
    destination = destination.resolve()
    value = lambda name: json.dumps(str(destination / name))  # noqa: E731
    return (
        "reference:\n"
        "  build: grch38\n"
        f"  genome: {value('Homo_sapiens_assembly38.fasta')}\n"
        f"  genome_dict: {value('Homo_sapiens_assembly38.dict')}\n"
        '  bwa_index_suffix: ".64"\n'
        f"  wgs_calling_regions: {value('wgs_calling_regions.interval_list')}\n"
        "  wgs_exclude_regions: null\n"
        f"  gnomad: {value('af-only-gnomad.hg38.vcf.gz')}\n"
        f"  population_vcf: {value('af-only-gnomad.hg38.vcf.gz')}\n"
        '  population_af_field: "AF"\n'
        f"  contamination_sites: {value('small_exac_common_3.hg38.vcf.gz')}\n"
        "  panel_of_normals: null\n"
        "  problematic_regions: null\n"
        "  low_mappability_regions: null\n"
        "  repeat_regions: null\n"
    )


def write_atomic(path: Path, content: str) -> None:
    temporary = path.with_name(f"{path.name}.partial")
    temporary.write_text(content)
    os.replace(temporary, path)


def provision(destination: Path, execute: bool) -> None:
    destination = destination.expanduser().absolute()
    required = remaining_bytes(destination)
    available = check_capacity(destination, required)
    print(f"Destination: {destination}")
    print(f"Pinned download total: {human_size(sum(item.size for item in RESOURCES))}")
    print(f"Remaining download: {human_size(required)}")
    print(f"Available capacity: {human_size(available)}")
    for resource in RESOURCES:
        print(f"  {resource.filename}\t{human_size(resource.size)}\t{resource.url}")
    print("\nConfiguration overlay:\n")
    print(config_text(destination), end="")
    if not execute:
        print("\nPlan only; rerun with --execute after reviewing storage and sources.")
        return

    destination.mkdir(parents=True, exist_ok=True)
    records = []
    for resource in RESOURCES:
        status = download(resource, destination / resource.filename)
        print(f"{status}: {resource.filename}")
        source_version = (
            "Broad GRCh38 v0" if resource.bucket == REFERENCE_BUCKET else "GATK somatic-hg38"
        )
        records.append(
            asdict(resource)
            | {
                "url": resource.url,
                "source_version": source_version,
                "accessed_at": datetime.now(timezone.utc).isoformat(),
                "status": status,
            }
        )

    interval_list = destination / "wgs_calling_regions.interval_list"
    create_interval_list(destination / "Homo_sapiens_assembly38.dict", interval_list)
    validate_resources(destination)
    overlay = destination / "organoid-pipeline.reference.yaml"
    write_atomic(overlay, config_text(destination))
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "bundle": "Broad GRCh38 v0 plus GATK somatic-hg38",
        "resources": records,
        "generated": [
            {
                "filename": interval_list.name,
                "source": "Homo_sapiens_assembly38.dict primary chr1-chr22,chrX,chrY",
                "size": interval_list.stat().st_size,
                "sha256": sha256(interval_list),
                "created_at": datetime.now(timezone.utc).isoformat(),
            },
            {
                "filename": overlay.name,
                "source": "provision_grch38.py",
                "size": overlay.stat().st_size,
                "sha256": sha256(overlay),
                "created_at": datetime.now(timezone.utc).isoformat(),
            },
        ],
    }
    write_atomic(
        destination / "resource-manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    print(f"\nProvisioning complete. Copy settings from: {overlay}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plan or download the pinned GRCh38 resources for this pipeline."
    )
    parser.add_argument(
        "--destination",
        required=True,
        type=Path,
        help="Absolute shared-filesystem destination visible from compute nodes",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Perform downloads; without this flag the command only prints a plan",
    )
    args = parser.parse_args()
    if not args.destination.expanduser().is_absolute():
        raise ValueError("--destination must be absolute")
    provision(args.destination, args.execute)


if __name__ == "__main__":
    main()
