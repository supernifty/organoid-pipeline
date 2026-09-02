# Low-pass organoid somatic SNV pipeline

Production-oriented Snakemake 9 workflow for mutations acquired between an ancestral early-passage monoclonal organoid baseline and later descendant organoids. GRCh38 WGS is the primary production profile; native GRCh37 WGS is supported for verified legacy alignments. Approximately 6× later samples, shared baselines, FASTQ and CRAM inputs, Mutect2, and Strelka2 are first-class.

The baseline is a biological time point, not an unrelated technical normal. A zero alternate count in a low-depth baseline is not proof that an allele was absent, so the workflow retains baseline depth/counts and makes filtering reversible.

## Workflow

```text
original FASTQ ─ FastQC ─ BWA-MEM ─ coordinate CRAM ─ mark duplicates ┐
validated CRAM ────────────────────────────────────────────────────────┤
                                                                     ├─ Mutect2 + contamination + orientation model ─┐
organoid + matched ancestral baseline + common WGS territory ────────┤                                                ├─ exact caller tiers
                                                                     └─ Strelka2 WGS ─────────────────────────────────┘
                                                                                                                        │
all organoid candidates ─ all cohort alignments ─ exact-allele recount ─ population + recurrence + evidence filters ────┤
                                                                                                                        └─ audit/retained/rejected + SBS96
```

Original FASTQs are aligned directly. There is no adapter trimming or quality trimming in the active DAG.

## Installation and dependency ownership

Pixi owns the outer environment, Snakemake, SLURM executor, and native bioinformatics tools. uv owns Python-only helpers and tests. GATK, Strelka2, and optional Somalier run in pinned containers through one shared runtime abstraction. PyYAML is deliberately present in both Python contexts: Snakemake parses the manifest inside Pixi, while unit tests import the same validator inside uv.

```bash
pixi install --locked
pixi run python-sync
pixi run test
```

Useful tasks are `python-sync`, `test`, `lint`, `format-check`, `dry-run`, `run-local`, and `run-slurm`. Keep both `pixi.lock` and `uv.lock` locked in production.

Pull configured images with:

```bash
./scripts/pull_images.sh
```

Docker is supported locally. Singularity/Apptainer is the intended cluster runtime. The default GATK and Strelka tags are inherited from the proven upstream pipeline; review and mirror them under site policy before a production release. Somalier and VEP images use immutable digests.

## Manifest

Copy `config/samples.example.yaml` to the ignored `config/samples.yaml`.

```yaml
samples:
  baseline_1:
    role: baseline
    donor: donor_1
    lineage: lineage_1
    cram: /data/baseline_1.cram
    crai: /data/baseline_1.cram.crai

  organoid_1:
    role: organoid
    donor: donor_1
    lineage: lineage_1
    condition: treated
    fastq_1: /data/organoid_1_R1.fastq.gz
    fastq_2: /data/organoid_1_R2.fastq.gz
    read_group:
      platform: ILLUMINA
      library: organoid_1
      unit: flowcell.lane

comparisons:
  organoid_1:
    baseline: baseline_1
```

Every organoid has exactly one baseline. One baseline may serve multiple descendants. Donor and lineage must match within each comparison. Identifiers must be filename-safe. A sample cannot mix FASTQ and CRAM inputs. FASTQs require both mates and complete read-group metadata.

Pre-aligned CRAMs must be coordinate sorted, indexed, duplicate marked, aligned against the identical FASTA, carry the complete ordered reference sequence dictionary, and use the configured `cram_version` (currently 3.0 for compatibility with pinned GATK 4.4/HTSJDK). Keep the CRAI at a standard samtools-discoverable path (`sample.cram.crai` or `sample.crai`); the manifest rejects arbitrary detached index paths. Preflight checks readable files, CRAM format, `SM`, sort order, exact dictionary equality, reference compatibility, and index usability. Use `bam_sample` only when a reviewed CRAM header name intentionally differs from the manifest key.

Strelka 2.9 explicitly accepts BAM or CRAM and documents `--normalBam`, `--tumorBam`, `--referenceFasta`, indexed bgzip `--callRegions`, and WGS defaults. Mutect2/GATK uses the reference-backed CRAM directly. The relevant upstream documentation is the [Strelka 2.9 user guide](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md) and [GATK Mutect2 workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2).

## GRCh37 and GRCh38 resources

The configured build must follow the alignment. Never relabel a BAM header, encode a GRCh37 BAM as CRAM with GRCh38, or use coordinate liftover as a substitute for read realignment. The preflight checks primary contig names and lengths, and downsampling additionally rejects mapped auxiliary contigs that the configured FASTA cannot represent.

### GRCh38

Edit `config/config.local.yaml` or supply a complete `--configfile`. All paths are configurable; nothing large is committed.

Required resources:

- GRCh38 FASTA, `.fai`, sequence dictionary, and BWA index including the ALT index when using the Broad assembly;
- WGS interval list and optional exclusions;
- AF-only gnomAD VCF for Mutect2 plus Tabix index;
- a population VCF with exact REF/ALT `AF` values plus index;
- biallelic common SNPs for contamination estimation plus index;
- optional independently constructed technical PoN;
- optional problematic, low-mappability, repeat, or segmental-duplication BED masks.

Every VCF and interval resource must use the same contig names, order, lengths, and GRCh38 assembly as the FASTA. The preflight compares dictionaries and the canonical GRCh38 primary-contig lengths. Record source version, URL, access date, checksum, and licence alongside locally managed resources. Broad GRCh38 resources are described by the [GATK resource bundle documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

The repository provides an explicit, plan-first provisioner for the pinned Broad GRCh38 v0 reference and GATK somatic-hg38 resources. Reference acquisition is intentionally separate from `pixi install`: the download is approximately 11.24 GiB, is shared across analyses, and belongs on a filesystem selected by the operator. First inspect the immutable object URLs, sizes, destination, and free-space calculation without changing the destination:

```bash
pixi run provision-grch38 --destination /absolute/shared/path/grch38
```

After reviewing the plan, start or resume the verified downloads explicitly:

```bash
pixi run provision-grch38 --destination /absolute/shared/path/grch38 --execute
```

The command downloads through `.partial` files, resumes when the server honors byte ranges, verifies the MD5 values published in Google Cloud Storage object metadata, atomically publishes complete files, constructs a primary-contig WGS interval list, validates resource dictionaries, and writes `resource-manifest.json`. It also writes `organoid-pipeline.reference.yaml`; copy its `reference:` mapping into the ignored `config/config.local.yaml`. The initial overlay deliberately uses the same AF-only gnomAD exact-allele VCF for Mutect2 and population filtering, avoiding a second much larger gnomAD download. A newer or more comprehensive exact-allele population resource can be configured later after compatibility review.

The configured paths should be canonical absolute host paths. Their containing directories are mounted into containers automatically; the resources do not need to be embedded in an image or copied into this repository.

### Native GRCh37 WGS

For verified legacy BAMs aligned to Broad-style GRCh37 with contigs `1`–`22`, `X`, and `Y`, use the separate plan-first provisioner. The official GRCh37 Mutect2 resource is distributed as a 13.1 GiB uncompressed VCF, so the complete pinned remote download is 21.08 GiB. From an empty destination, the provisioner reserves 5 GiB transformation headroom and requires approximately 31.3 GiB free after its 1.2× safety margin.

```bash
pixi run provision-grch37 --destination /absolute/shared/path/grch37
```

After reviewing the object generations, checksums, and capacity, execute with an appropriate number of compression threads:

```bash
pixi run provision-grch37 \
  --destination /absolute/shared/path/grch37 \
  --threads 8 \
  --execute
```

The command downloads the Broad `Homo_sapiens_assembly19` FASTA and BWA indexes plus the official GATK somatic-b37 resources. It verifies the published object checksums, injects the exact FASTA dictionary into the legacy VCF headers, BGZF-compresses and Tabix-indexes them, validates `INFO/AF` and reference compatibility, records derived SHA-256 checksums, and then removes the verified uncompressed VCF copies. Download and conversion progress is printed every 30 seconds. Interrupted downloads retain resumable `.partial` files.

Copy [config.grch37-wgs.example.yaml](config/config.grch37-wgs.example.yaml) to the ignored `config/config.local.yaml`, then replace its paths from `organoid-pipeline.grch37.reference.yaml`. Keep `mutational_signatures.reference_build`, `hotspots.reference_build`, `chromosomes`, and `analysis.wgs.contigs` on GRCh37. Use a fresh batch name: run management rejects reuse across builds.

The GRCh37 mode runs the same paired Mutect2, Strelka2, exact caller tiers, cohort recounting, reason-coded filtering, and SBS96 outputs as GRCh38. Its calls remain GRCh37-native. For GRCh38 production output, obtain original reads and realign them to GRCh38 rather than lifting the GRCh37 alignment.

Ancestral baselines are never automatically combined into a panel of normals. `reference.panel_of_normals` is optional and must point to an unrelated, reviewed technical PoN.

## Run

Validate local files, inspect the DAG, and then execute:

```bash
./scripts/check_sample_files.py config/samples.yaml
pixi run snakemake --dry-run --cores 1
pixi run run-local
```

For SLURM, install the locked environment on shared storage, copy `scripts/setup.sh.example` to the ignored `scripts/setup.sh`, and configure the local account, partition, ordered module prerequisites, container module, cache, and temporary paths. Some sites require a compiler module before the container-runtime module; follow the order reported by the site module system. Copy the configuration and manifest examples to their ignored local names and use absolute resource/input paths visible from compute nodes:

```bash
pixi install --locked
cp scripts/setup.sh.example scripts/setup.sh
cp config/config.local.example.yaml config/config.local.yaml
cp config/samples.example.yaml config/samples.yaml
```

Use Apptainer or Singularity on a compute node to pull the pinned images, then validate the manifest and DAG:

```bash
source scripts/setup.sh
./scripts/pull_images.sh apptainer  # or: singularity
pixi run python scripts/check_sample_files.py config/samples.yaml
pixi run snakemake --dry-run --cores 1 \
  --configfile config/config.local.yaml
```

Preview the persistent controller job without submitting it:

```bash
./scripts/run_driver_slurm.sh \
  --driver-dry-run --dry-run \
  --batch BATCH-preview \
  --samples config/samples.yaml
```

Submit a new batch, monitor it with the site's normal Slurm tools, and inspect the recorded state:

```bash
./scripts/run_driver_slurm.sh \
  --batch BATCH \
  --samples config/samples.yaml

squeue --me
pixi run python scripts/run_manager.py status --batch BATCH
```

For a small scheduler smoke run, use one organoid–baseline comparison and a reviewed small interval list through `reference.wgs_calling_regions`. Set one calling contig and one Mutect2 shard in the local overlay, and target that organoid's stringent catalog explicitly:

```yaml
container_runtime: apptainer
analysis:
  type: wgs
  wgs:
    contigs: [SMALL_CONTIG]
    scatter_count: 1
    max_concurrent_mutect2_shards: 1
reference:
  wgs_calling_regions: /absolute/shared/path/small.interval_list
germline:
  enabled: false
annotation:
  enabled: false
somalier:
  enabled: false
```

```bash
./scripts/run_driver_slurm.sh \
  --batch BATCH-mini \
  --samples config/samples.yaml \
  -- results/catalogs/ORGANOID_ID.stringent.vcf.gz
```

Resume failed or partial work without deleting outputs:

```bash
./scripts/run_driver_slurm.sh \
  --batch BATCH --resume \
  --samples config/samples.yaml
```

The profile uses Snakemake's native SLURM executor, scheduler resource limits, restart-safe outputs, per-rule logs, and benchmarks. Each batch records its effective configuration, samples, Git revision, launch history, and state under `runs/BATCH/`. Use `scripts/run_manager.py recover` only after the scheduler proves that a stale controller is no longer running. Guarded recovery refuses an active controller, clears the isolated batch's stale Snakemake lock with the pinned executable, and records the unlock outcome before permitting a later explicit `--resume`; do not run `snakemake --unlock` against an active batch. Long-running outputs publish from temporary paths after validation.

For planning, allow roughly the combined size of all input alignments plus 2–3× the largest active alignment for temporary CRAM, caller, and sorting work. Mutect2 uses approximately 20 Mb callable shards with a configurable concurrency cap and a 12-hour default walltime; Strelka requests 16 cores/32 GB and up to 24 hours for WGS; FASTQ alignment and optional HaplotypeCaller shards are capped at 24 hours. Actual wall time and scratch scale with callable territory, depth, filesystem, and scheduler. Check `df -h .` and site quotas before a real cohort or benchmark. Override a measured exception at submission with Snakemake syntax such as `--set-resources mutect2_chromosome:runtime=960` rather than restoring multi-day defaults globally.

## Outputs

Per organoid, `results/callers/` contains normalized indexed Mutect2 PASS, Strelka2 PASS, intersection, union, Mutect2-only, and Strelka2-only VCFs. `CALLER_SUPPORT` preserves provenance.

`results/cohort/candidates.union.vcf.gz` is the exact-allele candidate union. `allele_counts.tsv` contains depth, REF/ALT counts, VAF, strand counts, and mean observed base/mapping quality for every candidate/sample combination. Recounting uses one `samtools mpileup` pass per sample with mapping/base quality 20, default overlapping-mate removal, duplicate exclusion, and explicit normalized SNV/indel event matching. The [samtools mpileup manual](https://www.htslib.org/doc/samtools-mpileup.html) defines those semantics.

`results/catalogs/` provides:

- `audit`: every candidate and every annotation;
- `stringent`: excludes any configured population occurrence and other reason-coded failures;
- `sensitivity`: excludes population AF above 0.001 rather than every occurrence;
- `rejected`: candidates failing the stringent policy;
- `shared-lineage`: candidates called in multiple descendants of the lineage;
- `stage_counts`: counts through the main filters.

VCF and TSV outputs retain caller support, later and baseline evidence, exact population AF, total/same-lineage/unrelated recurrence, carriers, and `FILTER_REASONS`. The default later evidence is ALT ≥2 and VAF ≥0.20; baseline depth minimum is 6 and baseline ALT >1 is flagged. Exact recurrence across unrelated donors is excluded from stringent catalogs. Filters are configurable in `config/config.yaml`. Copy-number filtering is reserved for a future annotated stage.

`results/signatures/{sample}.sbs96.tsv` contains SBS96 counts for Mutect2, Strelka2, their intersection, and the stringent catalog. Signature fitting remains optional and separate; catalog completion is not evidence that either caller is biologically adequate at 6×.

QC includes original-read FastQC, samtools flagstat/stats, alignment and insert-size metrics, duplicate information, contamination, mosdepth mean autosomal depth, callable depth percentages, optional Somalier, and MultiQC. The approximately 6× expectation applies to organoids; deeper baselines are not warned against the same target.

## Tests and smoke dry run

```bash
pixi run test
pixi run lint
pixi run format-check
pixi run python tests/create_dryrun_fixture.py --build grch38 --output tmp/codex/dryrun-grch38
pixi run snakemake --dry-run --cores 1 --configfile tmp/codex/dryrun-grch38/config.yaml -- results/catalogs/O.stringent.vcf.gz
pixi run python tests/create_dryrun_fixture.py --build grch37 --output tmp/codex/dryrun-grch37
pixi run snakemake --dry-run --cores 1 --configfile tmp/codex/dryrun-grch37/config.yaml -- results/catalogs/O.stringent.vcf.gz
```

The fixtures are parse/DAG-only: placeholder CRAMs and resources must never be executed. Unit tests cover shared baselines, manifest rejection, build/reference mismatches, exact caller tiers, strand-aware SNV/indel recounting, SBS96 canonicalization, benchmark metrics, and absence of trimming.

## SEQC2

See `benchmarks/seqc2/README.md`. Ordinary tests never download SEQC2. The explicit workflow uses selected HCC1395/HCC1395BL alignments, deterministic template-preserving downsampling, measured achieved depth, tumour 6× versus baseline 6×/15×/30×, repeated seeds, a normal–normal negative control, v1.2 high-confidence regions/truth, exact-allele metrics, and the pipeline's SBS96 outputs.

The public truth and regions are hosted in the [NCBI SEQC2 somatic working-group release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/). Raw data are under SRA project `SRP162370`. Review remote sizes and available capacity before acquisition; do not download a complete release implicitly.

## Limitations

- Approximately 6× WGS has substantial allelic dropout. Sensitivity depends on clone fraction, local coverage, mapping, baseline depth, and caller model; retained counts must accompany interpretation.
- The strict cross-caller intersection sacrifices sensitivity and is one tier, not the only final result.
- A baseline is not a technical PoN. Shared biological alleles and low-depth uncertainty require lineage-aware interpretation.
- Copy-number-aware filters are not active until copy-number inputs are defined.
- Manta is not required for the initial SBS workflow; this may reduce Strelka indel sensitivity but does not block SNV output.
- SEQC2 validates technical performance on a cancer cell-line pair, not acquisition biology in monoclonal organoids.

The governing implementation contract and architecture assessment are in `specs/low-pass-organoid-somatic-snv.md`. This repository selectively derives from somatic pipeline commit `a533612`; large upstream data, results, caches, and unrelated changes were not copied.
