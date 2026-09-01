#!/usr/bin/env python3
"""Plan or provision pinned Broad GRCh37 resources for native WGS analysis."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
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
from provision_grch38 import (  # noqa: E402
    Resource,
    check_capacity,
    download,
    human_size,
    remaining_bytes,
    sha256,
    write_atomic,
)

REFERENCE_BUCKET = "gcp-public-data--broad-references"
REFERENCE_PREFIX = "hg19/v0"
SOMATIC_BUCKET = "gatk-best-practices"
SOMATIC_PREFIX = "somatic-b37"
TRANSFORM_HEADROOM = 5 * 2**30

CORE_RESOURCES = (
    Resource(
        "Homo_sapiens_assembly19.fasta",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta",
        "1575676516091053",
        3_140_756_381,
        "iGuhVZOT91hywc9FnrV/LQ==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.fai",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.fai",
        "1575676515993281",
        2_780,
        "/cCrZ59kYdeJgN4qLpfo8w==",
    ),
    Resource(
        "Homo_sapiens_assembly19.dict",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.dict",
        "1575676515936444",
        14_811,
        "zTsd3gQnr9ihFbJzepZxPg==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.64.amb",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.64.amb",
        "1575676515917590",
        6_597,
        "t81JYCU/zcoYQtX4ZLoeiQ==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.64.ann",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.64.ann",
        "1575676515912674",
        6_901,
        "Dvt1S0f1oG5572HCsgyWWQ==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.64.bwt",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.64.bwt",
        "1575676516115149",
        3_101_976_644,
        "nkYeNZ53yg7E7TAksRw4wg==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.64.pac",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.64.pac",
        "1575676515970832",
        775_494_142,
        "1CzPI6ZDbjjj+7nOyoTLpg==",
    ),
    Resource(
        "Homo_sapiens_assembly19.fasta.64.sa",
        REFERENCE_BUCKET,
        f"{REFERENCE_PREFIX}/Homo_sapiens_assembly19.fasta.64.sa",
        "1575676516011652",
        1_550_988_336,
        "PVOTRACm5FCrncCjbOOzFQ==",
    ),
)

SOURCE_VCFS = (
    (
        Resource(
            "af-only-gnomad.raw.sites.vcf",
            SOMATIC_BUCKET,
            f"{SOMATIC_PREFIX}/af-only-gnomad.raw.sites.vcf",
            "1501865611361723",
            14_062_778_952,
            "LNqQauwGDll57AGAWLppIg==",
        ),
        "af-only-gnomad.grch37.vcf.gz",
    ),
    (
        Resource(
            "small_exac_common_3.vcf",
            SOMATIC_BUCKET,
            f"{SOMATIC_PREFIX}/small_exac_common_3.vcf",
            "1501865613293081",
            3_177_496,
            "s/rtvulcq5GVCti1DEs3Xw==",
        ),
        "small_exac_common_3.grch37.vcf.gz",
    ),
)


def load_manifest(destination: Path) -> dict:
    path = destination / "resource-manifest.json"
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def generated_record(manifest: dict, filename: str) -> dict | None:
    return next(
        (item for item in manifest.get("generated", []) if item.get("filename") == filename),
        None,
    )


def derived_ready(destination: Path, final_name: str, manifest: dict) -> bool:
    for filename in (final_name, f"{final_name}.tbi"):
        path = destination / filename
        record = generated_record(manifest, filename)
        if not path.is_file() or not record:
            return False
        if path.stat().st_size != record.get("size") or sha256(path) != record.get("sha256"):
            return False
    return True


def create_interval_list(dictionary: Path, output: Path) -> None:
    order, lengths = read_sequence_dictionary(dictionary)
    primary = PRIMARY_CONTIG_LENGTHS["grch37"]
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


def bgzip_with_reference_dictionary(
    source: Path,
    output: Path,
    reference_fai: Path,
    threads: int,
    bgzip: str = "bgzip",
    tabix: str = "tabix",
) -> None:
    order, lengths = read_fai(reference_fai)
    temporary = output.with_name(f".{output.name}.partial.vcf.gz")
    temporary_index = Path(f"{temporary}.tbi")
    temporary.unlink(missing_ok=True)
    temporary_index.unlink(missing_ok=True)
    total = source.stat().st_size
    started = time.monotonic()
    last_report = started
    consumed = 0
    saw_header = False
    with temporary.open("wb") as compressed:
        process = subprocess.Popen(
            [bgzip, "-@", str(threads), "-c"],
            stdin=subprocess.PIPE,
            stdout=compressed,
        )
        assert process.stdin is not None
        try:
            with source.open("rb") as uncompressed:
                for line in uncompressed:
                    consumed += len(line)
                    if line.startswith(b"##contig=<"):
                        continue
                    if line.startswith(b"#CHROM"):
                        for contig in order:
                            process.stdin.write(
                                f"##contig=<ID={contig},length={lengths[contig]}>\n".encode()
                            )
                        saw_header = True
                    process.stdin.write(line)
                    now = time.monotonic()
                    if now - last_report >= 30:
                        elapsed = max(now - started, 0.001)
                        print(
                            f"transform progress: {source.name}: {human_size(consumed)} / "
                            f"{human_size(total)} ({100 * consumed / total:.1f}%), "
                            f"{consumed / elapsed / 2**20:.1f} MiB/s",
                            flush=True,
                        )
                        last_report = now
            process.stdin.close()
            return_code = process.wait()
        except BaseException:
            process.kill()
            process.wait()
            raise
    if return_code != 0:
        raise RuntimeError(f"bgzip failed for {source} with exit code {return_code}")
    if not saw_header:
        raise ValueError(f"VCF has no #CHROM header: {source}")
    subprocess.run([bgzip, "-t", str(temporary)], check=True)
    subprocess.run([tabix, "-f", "-p", "vcf", str(temporary)], check=True)
    subprocess.run([tabix, "-l", str(temporary)], check=True, stdout=subprocess.DEVNULL)
    os.replace(temporary, output)
    os.replace(temporary_index, Path(f"{output}.tbi"))


def vcf_has_af(path: Path) -> bool:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=AF,"):
                return True
            if line.startswith("#CHROM"):
                return False
    return False


def validate_resources(destination: Path) -> None:
    fai_order, fai_lengths = read_fai(destination / "Homo_sapiens_assembly19.fasta.fai")
    dictionary_order, dictionary_lengths = read_sequence_dictionary(
        destination / "Homo_sapiens_assembly19.dict"
    )
    if (fai_order, fai_lengths) != (dictionary_order, dictionary_lengths):
        raise ValueError("Reference FAI and sequence dictionary do not match")
    primary = PRIMARY_CONTIG_LENGTHS["grch37"]
    if any(fai_lengths.get(contig) != length for contig, length in primary.items()):
        raise ValueError("Reference dictionary does not match canonical GRCh37")
    interval_dictionary, intervals = read_interval_list(
        destination / "wgs_calling_regions.interval_list"
    )
    if interval_dictionary != dictionary_lengths:
        raise ValueError("WGS interval-list dictionary does not match the reference")
    intervals = validate_intervals(intervals, fai_order, fai_lengths)
    if [contig for contig, _, _ in intervals] != list(primary):
        raise ValueError("WGS interval list does not contain each primary contig exactly once")
    for name in (
        "af-only-gnomad.grch37.vcf.gz",
        "small_exac_common_3.grch37.vcf.gz",
    ):
        vcf_dictionary = read_vcf_dictionary(destination / name)
        missing = [contig for contig in primary if contig not in vcf_dictionary]
        if missing:
            raise ValueError(f"{name} lacks primary contigs: {', '.join(missing)}")
    if not vcf_has_af(destination / "af-only-gnomad.grch37.vcf.gz"):
        raise ValueError("AF-only gnomAD does not declare the required INFO/AF field")


def config_text(destination: Path) -> str:
    destination = destination.resolve()

    def value(name: str) -> str:
        return json.dumps(str(destination / name))

    contigs = json.dumps(list(PRIMARY_CONTIG_LENGTHS["grch37"]))
    return (
        "analysis:\n"
        "  type: wgs\n"
        "  wgs:\n"
        f"    contigs: {contigs}\n"
        "reference:\n"
        "  build: grch37\n"
        f"  genome: {value('Homo_sapiens_assembly19.fasta')}\n"
        f"  genome_dict: {value('Homo_sapiens_assembly19.dict')}\n"
        '  bwa_index_suffix: ".64"\n'
        "  regions_source_build: grch37\n"
        f"  wgs_calling_regions: {value('wgs_calling_regions.interval_list')}\n"
        "  wgs_exclude_regions: null\n"
        f"  gnomad: {value('af-only-gnomad.grch37.vcf.gz')}\n"
        f"  population_vcf: {value('af-only-gnomad.grch37.vcf.gz')}\n"
        '  population_af_field: "AF"\n'
        f"  contamination_sites: {value('small_exac_common_3.grch37.vcf.gz')}\n"
        "  panel_of_normals: null\n"
        "  problematic_regions: null\n"
        "  low_mappability_regions: null\n"
        "  repeat_regions: null\n"
        "mutational_signatures:\n"
        "  reference_build: grch37\n"
        "hotspots:\n"
        "  reference_build: grch37\n"
        "  resource: config/hotspots.tsv\n"
        f"chromosomes: {contigs}\n"
    )


def generated_file_record(path: Path, source: str) -> dict:
    return {
        "filename": path.name,
        "source": source,
        "size": path.stat().st_size,
        "sha256": sha256(path),
        "created_at": datetime.now(timezone.utc).isoformat(),
    }


def provision(destination: Path, execute: bool, threads: int) -> None:
    destination = destination.expanduser().absolute()
    previous_manifest = load_manifest(destination)
    source_needed = [
        (resource, final_name)
        for resource, final_name in SOURCE_VCFS
        if not derived_ready(destination, final_name, previous_manifest)
    ]
    selected = CORE_RESOURCES + tuple(resource for resource, _ in source_needed)
    remaining = remaining_bytes(destination, selected)
    working = remaining + (TRANSFORM_HEADROOM if source_needed else 0)
    available = check_capacity(destination, working)
    remote_total = sum(item.size for item in CORE_RESOURCES) + sum(
        item.size for item, _ in SOURCE_VCFS
    )
    print(f"Destination: {destination}")
    print(f"Pinned remote objects: {human_size(remote_total)}")
    print(f"Remaining download: {human_size(remaining)}")
    print(
        "Temporary VCF transformation headroom: "
        f"{human_size(TRANSFORM_HEADROOM if source_needed else 0)}"
    )
    print(f"Available capacity: {human_size(available)}")
    for resource in selected:
        print(f"  {resource.filename}\t{human_size(resource.size)}\t{resource.url}")
    print("\nConfiguration overlay:\n")
    print(config_text(destination), end="")
    if not execute:
        print("\nPlan only; rerun with --execute after reviewing storage and sources.")
        return

    destination.mkdir(parents=True, exist_ok=True)
    previous_records = {item["filename"]: item for item in previous_manifest.get("resources", [])}
    records = dict(previous_records)
    for resource in selected:
        status = download(resource, destination / resource.filename)
        print(f"{status}: {resource.filename}")
        source_version = (
            "Broad GRCh37/hg19 v0" if resource.bucket == REFERENCE_BUCKET else "GATK somatic-b37"
        )
        records[resource.filename] = asdict(resource) | {
            "url": resource.url,
            "source_version": source_version,
            "accessed_at": datetime.now(timezone.utc).isoformat(),
            "status": status,
            "retained": resource in CORE_RESOURCES,
        }

    interval_list = destination / "wgs_calling_regions.interval_list"
    create_interval_list(destination / "Homo_sapiens_assembly19.dict", interval_list)
    for resource, final_name in source_needed:
        print(f"Converting and indexing: {resource.filename} -> {final_name}", flush=True)
        bgzip_with_reference_dictionary(
            destination / resource.filename,
            destination / final_name,
            destination / "Homo_sapiens_assembly19.fasta.fai",
            threads,
        )
    validate_resources(destination)

    overlay = destination / "organoid-pipeline.grch37.reference.yaml"
    write_atomic(overlay, config_text(destination))
    generated_paths = [interval_list, overlay]
    for _, final_name in SOURCE_VCFS:
        generated_paths.extend((destination / final_name, destination / f"{final_name}.tbi"))
    generated = [generated_file_record(path, "provision_grch37.py") for path in generated_paths]
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "bundle": "Broad GRCh37/hg19 v0 plus GATK somatic-b37",
        "resources": [records[name] for name in sorted(records)],
        "generated": generated,
    }
    write_atomic(
        destination / "resource-manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    for resource, _ in SOURCE_VCFS:
        source = destination / resource.filename
        if source.exists():
            source.unlink()
            print(f"removed verified uncompressed source after conversion: {source.name}")
    print(f"\nProvisioning complete. Copy settings from: {overlay}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plan or download pinned GRCh37 resources for native WGS analysis."
    )
    parser.add_argument(
        "--destination",
        required=True,
        type=Path,
        help="Absolute shared-filesystem destination visible from compute nodes",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Perform downloads; without this flag the command only prints a plan",
    )
    args = parser.parse_args()
    if not args.destination.expanduser().is_absolute():
        raise ValueError("--destination must be absolute")
    if args.threads <= 0:
        raise ValueError("--threads must be positive")
    provision(args.destination, args.execute, args.threads)


if __name__ == "__main__":
    main()
