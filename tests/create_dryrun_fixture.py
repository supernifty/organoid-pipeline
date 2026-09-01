#!/usr/bin/env python3
"""Create a tiny parse-only fixture for a complete GRCh37 or GRCh38 WGS DAG."""

import argparse
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("tmp/codex/dryrun"))
    parser.add_argument("--build", choices=("grch37", "grch38"), default="grch38")
    args = parser.parse_args()
    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    contig = "1" if args.build == "grch37" else "chr1"
    reference = root / "genome.fa"
    reference.write_text(f">{contig}\nA\n")
    Path(f"{reference}.fai").write_text(f"{contig}\t1000\t0\t0\t0\n")
    (root / "genome.dict").write_text(f"@HD\tVN:1.6\n@SQ\tSN:{contig}\tLN:1000\n")
    (root / "territory.interval_list").write_text(
        f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{contig}\tLN:1000\n{contig}\t1\t1000\t+\tCALLABLE\n"
    )
    extensions = ["amb", "ann", "bwt", "pac", "sa"]
    if args.build == "grch38":
        extensions.append("alt")
    for extension in extensions:
        Path(f"{reference}.64.{extension}").touch()
    for name in (
        "gnomad.vcf.gz",
        "gnomad.vcf.gz.tbi",
        "population.vcf.gz",
        "population.vcf.gz.tbi",
        "sites.vcf.gz",
        "sites.vcf.gz.tbi",
        "B.cram",
        "B.cram.crai",
        "O.cram",
        "O.cram.crai",
    ):
        (root / name).touch()
    samples = {
        "samples": {
            "B": {
                "role": "baseline",
                "donor": "D",
                "lineage": "L",
                "cram": str(root / "B.cram"),
                "crai": str(root / "B.cram.crai"),
            },
            "O": {
                "role": "organoid",
                "donor": "D",
                "lineage": "L",
                "cram": str(root / "O.cram"),
                "crai": str(root / "O.cram.crai"),
            },
        },
        "comparisons": {"O": {"baseline": "B"}},
    }
    (root / "samples.yaml").write_text(yaml.safe_dump(samples, sort_keys=False))
    config = yaml.safe_load(Path("config/config.yaml").read_text())
    config["container_runtime"] = "apptainer"
    config["analysis"]["wgs"]["contigs"] = [contig]
    config["analysis"]["wgs"]["scatter_count"] = 1
    config["chromosomes"] = [contig]
    config["reference"].update(
        {
            "build": args.build,
            "regions_source_build": args.build,
            "genome": str(reference),
            "genome_dict": str(root / "genome.dict"),
            "wgs_calling_regions": str(root / "territory.interval_list"),
            "gnomad": str(root / "gnomad.vcf.gz"),
            "population_vcf": str(root / "population.vcf.gz"),
            "contamination_sites": str(root / "sites.vcf.gz"),
        }
    )
    config["mutational_signatures"]["reference_build"] = args.build
    config["hotspots"]["reference_build"] = args.build
    config["run_management"] = {
        "samples_file": str(root / "samples.yaml"),
        "config_file": str(root / "config.yaml"),
    }
    config["storage"] = {
        "tmp_dir": str(root / "work"),
        "local_scratch": str(root / "scratch"),
        "pon_genomicsdb_scratch_mb": 1024,
    }
    path = root / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    print(path)


if __name__ == "__main__":
    main()
