# Release Notes

## Unreleased

- Added native GRCh37 WGS support for verified legacy alignments, including a
  plan-first, checksummed Broad/GATK resource provisioner and a site-neutral
  configuration overlay.
- Added fail-closed alignment/reference checks before deterministic
  downsampling, plus complete GRCh37 and GRCh38 DAG dry runs.
- Fixed source-BAM validation for current samtools `idxstats` syntax and made
  external-command failures report their captured stderr.
- Fixed the same `idxstats` incompatibility in production CRAM preflight and
  canonicalized helper-script paths across shared-filesystem aliases.
- Forced deterministic downsampling outputs to configured CRAM 3.0 for pinned
  GATK compatibility and made preflight reject incompatible cached CRAMs.
- Added resumable download progress reporting to both reference provisioners.
- Expanded the site-neutral SLURM runbook with locked setup, a one-pair/small-territory scheduler smoke run, driver preview, submission, monitoring, and recovery.
- Added first-class `container_runtime: apptainer` support and explicit binding of external reference directories.
- Stopped the cluster setup template from re-running `pixi install` inside every controller launch.
- Established the GRCh38 WGS organoid workflow from upstream commit `a533612`.
- Added explicit baseline/organoid roles, lineage-safe shared baselines, complete FASTQ read groups, and CRAM preflight validation.
- Removed trimming from the active DAG and dependency set; original reads are aligned directly.
- Disabled automatic biological-baseline PoN construction and added scattered F1R2 orientation modeling.
- Added exact caller tiers, cohort recounting, population/recurrence/baseline annotations, reversible reason-coded catalogs, and four-source SBS96 counts.
- Added Pixi/uv dependency ownership, locked Python tooling, unit tests, a complete synthetic DAG dry-run fixture, and explicit SEQC2 benchmark resources/downsampling support.

- Reworked alignment into directly coordinate-sorted CRAM streams: paired
  alignment uses 24 BWA and 8 sort threads, unpaired streams use 4 + 4, and the
  merged sorted CRAM goes directly to single-CPU MarkDuplicates with a 44 GB
  heap. Removed the separate whole-CRAM sort.
- Added the validating, restart-safe legacy BAM importer and checksummed import
  manifest without modifying legacy outputs.
- Added opt-in offline VEP 116 annotation with digest-pinned containers,
  fail-closed resource metadata, detailed indexed VCFs, PASS-only callable-base
  burden, normalized recurrence, and annotation-resource reporting.
- Replaced WGS per-base `samtools depth` streaming with bounded mosdepth output;
  added optional exon/gene 10×/20×/50×/100× coverage and separately configurable
  tumour/normal warnings.
- Added optional container-isolated Somalier extraction, shared-normal groups,
  relatedness, ancestry, and stable mismatch/swap flags.
- Extended MultiQC dependencies, provenance, batch fingerprints, and selective
  reuse. Cohort aggregates and cohort-level Somalier outputs are never reused
  after cohort changes.
- Per-sample annotated TSVs and genes-of-interest reports are intentionally not
  part of this release.

- Added default per-normal HaplotypeCaller gVCF generation and hard-filtered
  germline SNP/short-indel VCFs for alignment-backed runs.
- Added configurable SNP and INDEL/MIXED filters, Tabix sidecars, SLURM
  resources and concurrency control, provenance rows, and selective run-manager
  reuse.
- Final-somatic-VCF-only runs continue to skip all alignment and germline
  calling work.

## 0.3.0

Release `0.3.0` focuses on cluster-run robustness and safer restart behavior.

### Highlights

- Hardened CRAM, VCF, QC, and aggregate output publication with explicit sidecar
  outputs, temp-file publishing, and post-write validation.
- Declared VCF index sidecars as Snakemake outputs so missing or failed indexes
  trigger reruns.
- Added gzip, tabix, CRAM quickcheck, and non-empty output checks before
  publishing key outputs.
- Made cluster wrapper launches pass `--rerun-incomplete` explicitly.

### Verification

Recommended release checks:

```bash
pixi install --locked
bash tests/test_pipeline.sh
./scripts/check_sample_files.py config/samples.yaml
pixi run snakemake -n
./scripts/run_cluster.sh -n
./scripts/run_driver_slurm.sh --driver-dry-run results/signatures/signatures.done
```

The Snakemake dry run requires local reference and sample inputs. The smoke test
can run without full production data, but Docker-specific checks are skipped if
Docker is unavailable.

### Known Limitations

- `--rerun-incomplete` only covers outputs Snakemake knows are incomplete; the
  validation checks added in this release are the main protection against
  syntactically present but truncated outputs.
- Cluster validation still depends on site-specific `scripts/setup.sh`,
  container availability, filesystem latency, and scheduler resource limits.

## 0.2.0

Release `0.2.0` builds on the `v0.1` pipeline release and focuses on making the
workflow easier to run on collaborator data, especially when samples have
already been aligned.

### Highlights

- Added support for starting samples from existing coordinate-sorted,
  duplicate-marked CRAMs, with optional `.crai` paths.
- Added optional `bam_sample` entries for CRAMs whose header `SM` sample names
  differ from workflow sample keys.
- Allowed FASTQ-backed and CRAM-backed samples in the same `config/samples.yaml`.
- Added sample-path validation with `scripts/check_sample_files.py`.
- Added hotspot aggregation output at `results/aggregate/hotspots.tsv`.
- Documented clean-run preparation so previous `results/`, temporary files, and
  Snakemake metadata do not affect a new run.
- Expanded release validation guidance, including sample-file checks and SLURM
  dry-run previews.

### Expected Inputs

- `config/samples.yaml` copied from `config/samples.example.yaml` and edited for
  local tumour-normal pairs.
- For FASTQ-backed samples, one R1/R2 FASTQ pair per sample. Multi-lane data
  should still be concatenated per sample/read before running this version.
- For CRAM-backed samples, duplicate-marked CRAMs aligned to the configured
  reference, plus either explicit `crai` entries or default `{cram}.crai`
  indexes. Set `bam_sample` when the CRAM header `SM` name differs from the
  sample key.
- GRCh37/hg19 reference FASTA, FASTA index, sequence dictionary, BWA index,
  target regions BED, and GATK b37 gnomAD resource files at the configured
  paths.
- Singularity SIF images or Docker access for GATK and Strelka.

### Verification

Recommended release checks:

```bash
pixi install --locked
bash tests/test_pipeline.sh
./scripts/check_sample_files.py config/samples.yaml
pixi run snakemake -n
```

For SLURM deployments, also preview cluster submission:

```bash
./scripts/run_cluster.sh -n
./scripts/run_driver_slurm.sh --driver-dry-run results/signatures/signatures.done
```

The Snakemake dry run requires local reference and sample inputs. The smoke test
can run without full production data, but Docker-specific checks are skipped if
Docker is unavailable.

### Known Limitations

- The pipeline is configured for GRCh37/hg19/b37 resources.
- Sample configuration supports one paired FASTQ set per FASTQ-backed sample.
- CRAM-backed inputs are assumed to be coordinate-sorted, duplicate-marked,
  indexed, and aligned to the configured reference genome.
- Reference genome, target regions, test data, and container images are external
  assets and are not committed to the repository.
- Generated panel-of-normals output can be reused on reruns, but large cohorts
  should consider providing a curated external PON.
- The provenance table captures configured dependency constraints, not every
  resolved package version from `pixi.lock`.

## 0.1.0

Initial collaborator-facing release candidate for the GRCh37/hg19 somatic
variant-calling pipeline.

### Scope

- Paired-end Illumina tumour-normal FASTQ input.
- Optional Trimmomatic read trimming.
- BWA-MEM alignment and duplicate marking.
- Somatic variant calling with Strelka2 and GATK Mutect2.
- Exact-allele consensus VCF generation from normalized Mutect2 and Strelka
  calls.
- Final consensus filtering by tumour allele fraction, tumour depth, and
  matched-normal depth.
- Mutational signature summaries for configured tissue signatures.
- Colibactin extended-context summary.
- MultiQC and config-only provenance output.
- Local execution and Snakemake 9 SLURM profile execution.

### Expected Inputs

- `config/samples.yaml` copied from `config/samples.example.yaml` and edited for
  local tumour-normal pairs.
- One R1/R2 FASTQ pair per sample. Multi-lane data should be concatenated per
  sample/read before running this version.
- GRCh37/hg19 reference FASTA at `resources/reference/genome.fa`.
- FASTA index, sequence dictionary, and BWA index for the configured reference.
- A bgzip-compressed, tabix-indexed regions BED at
  `resources/reference/regions.bed.gz`.
- GATK b37 gnomAD resource VCF and index at the configured paths.
- Singularity SIF images or Docker access for GATK and Strelka.

### Verification

Recommended pre-release checks:

```bash
pixi install --locked
bash tests/test_pipeline.sh
pixi run snakemake -n
```

The Snakemake dry run requires local reference and sample inputs. The smoke test
can run without full production data, but Docker-specific checks are skipped if
Docker is unavailable.

### Known Limitations

- The pipeline is configured for GRCh37/hg19/b37 resources.
- Sample configuration currently supports one paired FASTQ set per sample.
- Reference genome, target regions, test data, and container images are external
  assets and are not committed to the repository.
- Generated panel-of-normals output can be reused on reruns, but large cohorts
  should consider providing a curated external PON.
- The provenance table captures configured dependency constraints, not every
  resolved package version from `pixi.lock`.

### Intended Next Release Work

- Support multi-lane FASTQ inputs directly in `config/samples.yaml`, instead of
  requiring collaborators to concatenate lanes upstream.
- Add rule-level Snakemake logs for easier debugging, especially for SLURM
  execution where terminal output is harder to follow.
