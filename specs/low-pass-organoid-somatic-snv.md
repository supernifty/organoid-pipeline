# Low-pass organoid somatic SNV pipeline

Status: accepted for implementation

Source specification: user request dated 2026-09-01

Upstream implementation base: `/Volumes/redback/cog/src/somatic_pipeline/pipeline` at commit `a533612`

## Problem and interpretation

The pipeline analyses low-pass whole-genome sequencing from monoclonal organoids. Each later organoid is compared with exactly one ancestral early-passage baseline. The desired variants are mutations acquired after that baseline, not the organoid lineage's complete lifetime catalog.

At approximately 6× depth, absence of an alternate read in the baseline is weak evidence. The pipeline must therefore preserve baseline depth and exact-allele counts, expose filtering reasons, and retain caller-specific and audit catalogs rather than presenting a caller intersection as complete truth.

## Goals

- Accept a cohort containing FASTQ and/or pre-aligned CRAM samples on GRCh38.
- Permit one baseline to be shared by multiple descendant organoids and support multiple donors or lineages.
- Align original paired FASTQs without adapter or quality trimming.
- Run paired Mutect2 and Strelka2 WGS calls for each organoid–baseline comparison.
- Produce normalized caller-specific, intersection, union, and caller-only exact-allele tiers.
- Recount every cohort candidate in relevant samples and preserve quantitative evidence.
- Apply configurable, reason-coded population, recurrence, baseline, later-sample, and region filters.
- Produce retained, rejected, shared-lineage, and complete audit catalogs.
- Produce caller-specific and stringent SBS96 catalogs and stage counts.
- Provide a small synthetic integration/benchmark test and opt-in SEQC2 benchmark support.
- Remain reproducible locally and on SLURM through Pixi, uv where useful, and pinned containers.

## Non-goals

- Downloading or committing reference genomes, FASTQs, CRAMs, BAMs, generated results, or the full SEQC2 dataset.
- Automatically constructing a technical panel of normals from ancestral baselines.
- Claiming biological validation from SEQC2 or from workflow completion.
- Adapter or quality trimming.
- Making signature fitting part of the required variant-calling path.
- Requiring copy-number inputs in the initial implementation.

## Architecture decision

Use this repository as a selective, source-tracked fork of the existing somatic pipeline rather than adding the workflow directly to the existing repository.

Reasons:

1. The upstream tracked source is small and reusable, while its 25 GB working tree contains local references, results, environments, archives, and unrelated untracked work that must not be copied.
2. The organoid workflow changes the domain model from generic tumour/normal pairs to explicit roles, lineage relationships, shared baselines, cohort recurrence, and cross-sample recounting.
3. Several upstream defaults are unsafe for this use case: trimming remains reachable, baselines can be used to create a PoN automatically, the primary somatic output is an intersection, and depth filters target deeper data.
4. Isolation prevents disruption of the existing WES/GRCh37 production workflow and its unrelated local changes.
5. Recording the upstream commit and keeping recognizable modules allows generally useful fixes to be proposed upstream later.

Only Git-tracked source will be imported. Existing generated data, reference resources, caches, environments, results, and untracked files will not be copied. The upstream repository will not be modified.

## Reusable upstream components

- Snakemake 9 modular workflow structure and native SLURM profile.
- FASTQ FastQC, BWA-MEM alignment, coordinate-sorted CRAM output, duplicate marking, and indexing.
- Container runtime, path binding, and image-selection helpers for Docker and Singularity/Apptainer.
- WGS territory preparation, interval scatter/gather, GRCh37/GRCh38 contig-length validation, resource-dictionary validation, and analysis manifest.
- Paired Mutect2 shards, F1R2 collection, stats gathering, contamination estimation, filtering, and normalization.
- Strelka2 somatic WGS configuration, execution, AF annotation, and normalization.
- Exact-allele intersection code and its tests.
- FastQC, Picard metrics, mosdepth coverage, Somalier, MultiQC, provenance, signature counting/fitting, test scaffolding, and run wrappers.
- Pixi lock, container pins, restart-safe publication patterns, and existing unit/shell test conventions.

## Gaps in the upstream implementation

- `Snakefile` explicitly rejects `reference.build: grch38` with `analysis.type: wgs`.
- GRCh38 explicitly rejects CRAM and final-VCF inputs.
- The manifest uses `tumours: tumour: normal` and lacks roles, donors, lineages, replicate metadata, and read-group fields.
- CRAM validation checks paths but not header sample identity, coordinate sorting, assembly compatibility, or index usability as a first-class preflight.
- FASTQ alignment has only a minimal read group and still conditionally depends on trimming.
- Normal samples may automatically become a Mutect2 PoN, which is prohibited for biological baselines.
- There is no explicit orientation-bias model in the inspected calling path.
- The existing final path centers on the intersection and does not produce a provenance-preserving union and both caller-only tiers.
- Current filtering drops records and uses tumour depth/AF defaults designed for deeper data.
- Exact population-allele annotation, cohort recurrence, cross-sample recounting, baseline evidence, and reason-coded retained/rejected catalogs are absent.
- Existing signatures center on the intersection and do not provide the four required SBS96 sources and stage counts.
- QC lacks a complete role-aware low-pass reporting contract.
- No SEQC2 acquisition/downsampling/benchmark workflow exists.

## Dependency ownership

- Pixi owns the outer reproducible environment, Snakemake 9, native tools (BWA, samtools, bcftools, htslib, bedtools, mosdepth, FastQC, MultiQC), and user-facing project tasks.
- uv owns Python-only runtime and test dependencies through `pyproject.toml` and `uv.lock`. Pixi invokes uv using its pinned `uv` executable. Python libraries are not independently declared in both systems unless a documented runtime constraint requires it.
- Containers own GATK, Strelka2, and Somalier. Images are pinned to immutable digests where upstream images support them, and runtime/bind logic remains centralized.

The imported project currently declares Python packages through Pixi's PyPI integration. Migration to uv must be atomic: add the Python project and lock, remove duplicated PyPI declarations from Pixi, and expose `uv sync --locked` plus `uv run` tasks through Pixi.

## Input contract

The YAML manifest has two top-level mappings: `samples` and `comparisons`.

Each sample requires:

- `role`: `baseline` or `organoid`;
- `donor` and `lineage`;
- exactly one input kind: paired `fastq_1`/`fastq_2`, or `cram` with a conventional samtools-discoverable CRAI path that may be stated explicitly;
- for FASTQ input, `read_group.platform`, `read_group.library`, and `read_group.unit`;
- optional `batch`, `condition`, and `replicate`.

Each comparison is keyed by organoid sample and defines exactly one `baseline`.

Validation fails before execution when:

- identifiers are duplicated or not filename-safe;
- roles are missing or inconsistent with comparisons;
- a comparison references an undefined sample;
- the organoid and baseline donor or lineage differ;
- inputs are mixed, paired FASTQs are incomplete/unreadable, or a required CRAI is unavailable;
- required read-group values are missing;
- an organoid lacks exactly one comparison.

Preflight validation of CRAMs uses samtools against the configured FASTA and fails on an unusable or detached index, non-coordinate sort, conflicting `SM` values, unexpected sample name without an explicit override, reference contig/length mismatch, or absent build evidence. The pipeline uses CRAM directly because the configured Mutect2, Strelka2, samtools, Picard, and mosdepth paths support reference-backed CRAM input.

## Reference contract

The primary profile is GRCh38 WGS with one consistent contig convention. Configuration provides FASTA, FAI, sequence dictionary, BWA indexes, AF-only gnomAD for Mutect2, exact-allele population VCF, contamination sites, optional technical PoN, optional masks, WGS territory, and exclusions.

Preflight validates path readability, sidecar indexes, contig names/order/lengths, VCF dictionaries, interval bounds/order, build identity, and resource compatibility. Provenance records configuration, tool/container identities, reference checksums, resource metadata, and commands.

Reference acquisition is a separate, explicit provisioning operation and is never an implicit side effect of `pixi install`, environment activation, ordinary tests, or pipeline execution. A tracked, site-neutral provisioning script must obtain one internally consistent Broad GRCh38 bundle containing the FASTA and its FAI, dictionary, BWA indexes, AF-only gnomAD Mutect2 resource, contamination-sites VCF, and indexes. It must also prepare a full-WGS interval list from the same dictionary and provide an exact-allele population VCF with an `AF` field. The script must:

- default to a non-mutating plan that reports the selected versions, URLs, expected sizes where discoverable, destination, and available capacity;
- require an explicit execution flag before downloading or generating anything;
- download into partial files, support safe restart, verify authoritative checksums when published, and atomically publish completed files;
- record a machine-readable manifest containing source URL, version, size, checksum, and access time for every acquired or generated artifact;
- reject incompatible or incomplete resources and run the repository reference preflight after preparation;
- accept an explicit destination so private cluster paths never appear in tracked source; and
- print a ready-to-copy local YAML overlay containing canonical absolute host paths, which the existing container binding abstraction makes visible inside Apptainer/Singularity or Docker.

The default first-run population resource may reuse the compatible AF-only gnomAD exact alleles when it exposes exact REF/ALT allele frequencies. A separate larger population resource can be substituted later through configuration. Estimated download and installed sizes must be documented so the operator can make the storage decision before execution.

### Native GRCh37 WGS pilot

Pre-aligned inputs must remain on their verified assembly. For a source organoid–baseline pair aligned to GRCh37, the pipeline supports a native GRCh37 WGS execution path rather than relabelling, coordinate-lifting, or encoding the alignments against GRCh38. The GRCh37 path must:

- require the canonical primary-contig lengths and the input convention `1`–`22`, `X`, and `Y`;
- use one pinned, internally consistent GRCh37 FASTA, FAI, dictionary, BWA indexes, Mutect2 germline resource, exact-allele population resource, contamination resource, and primary-contig WGS territory;
- validate each input BAM/CRAM against that exact FASTA before downsampling or calling;
- run the same paired Mutect2, Strelka2, caller-tier, filtering, recounting, and SBS96 vertical slice as GRCh38;
- record `grch37` in configuration, provenance, and signature-context metadata; and
- keep build-specific outputs isolated in a fresh batch so they cannot be mixed with GRCh38 results.

A separate plan-first GRCh37 provisioner follows the same restart, capacity, checksum, atomic-publication, manifest, and configuration-overlay requirements as the GRCh38 provisioner. The private first-pair preparation reads the build from its local configuration and must never assume GRCh38 in validation messages or provenance. Native GRCh37 output is suitable for this technical pilot and comparison with GRCh37 truth/catalogs; GRCh38 production analysis still requires original reads realigned to GRCh38 or a genuinely GRCh38-aligned input.

No technical PoN is generated by default. Automatic PoN construction from configured baselines is forbidden. An externally supplied technical PoN is optional.

## Workflow behavior

### FASTQ processing

Original FASTQs receive FastQC and are passed directly to BWA-MEM. Alignment read groups include `ID`, `SM`, `PL`, `LB`, and `PU`. Output is coordinate sorted, duplicate marked, validated, and indexed. Trimming is absent from the active DAG, configuration, targets, documentation, and new-workflow tests.

### Calling

For every comparison, Mutect2 receives the organoid as tumour and its early-passage baseline as normal. WGS calls are scattered over the common callable territory. F1R2 tarballs and shard stats are gathered; orientation bias, contamination, and `FilterMutectCalls` are applied. Raw, unfiltered, filtered, normalized, and indexed audit artifacts are retained as appropriate.

Strelka2 runs in somatic WGS mode for the same pair and territory. Its normalized indexed SNV and indel outputs retain SomaticEVS and quality fields. Manta is not required for the initial SBS vertical slice.

### Exact-allele tiers

Records are split and normalized before exact-key comparison using chromosome, position, REF, and ALT. Per organoid outputs are:

- Mutect2 PASS;
- Strelka2 PASS SNVs and indels;
- intersection;
- union with explicit caller support;
- Mutect2-only;
- Strelka2-only.

### Recounting and annotations

The union of exact candidate alleles is recounted in the matched baseline, its descendant organoids, and optionally all cohort samples. The method and caller/tool version are recorded. Output retains total, reference, alternate, forward/reverse counts where available, VAF, useful base/mapping-quality summaries, and caller support.

Filtering annotates rather than destroying the only record. Every audit candidate records semicolon-delimited reason codes and passes through explicit stages:

1. caller PASS/support;
2. later alternate reads, default ≥2;
3. later VAF, default ≥0.20;
4. baseline depth;
5. baseline alternate evidence;
6. exact population allele, with strict any-observation and sensitivity AF >0.001 catalogs;
7. normalized exact-allele recurrence;
8. optional masks;
9. optional copy-number logic when later configured.

Outputs include retained stringent, retained sensitivity, rejected, shared-lineage, and full audit VCF/TSV catalogs. Recurrence includes total, within-lineage, and unrelated-donor counts plus carrier names. Unrelated recurrence is excluded from signature-ready output by default; within-lineage shared variants are separated but preserved.

### QC and signatures

QC includes FastQC for FASTQs, flagstat, samtools stats, duplicate, insert-size and alignment metrics, depth summaries/distribution, callable-genome estimates, contamination, optional Somalier identity/relatedness, and cohort MultiQC. Expected-depth warnings are role-specific.

SBS96 counts are generated independently for Mutect2 PASS, Strelka2 PASS, intersection, and final stringent SNVs, with mutation counts before and after filter stages. Signature fitting remains optional and names its signature database/version with diagnostics.

## SEQC2 test and benchmark

Normal tests use tiny synthetic inputs and truth data. They do not download SEQC2.

An explicit opt-in benchmark command prepares selected HCC1395/HCC1395BL resources, verifies checksums, performs deterministic pair-preserving downsampling with configurable seeds, measures achieved depth, runs tumour at approximately 6× against baselines at approximately 6×, 15×, and 30×, includes a normal–normal negative control, and evaluates within high-confidence regions.

Reported metrics are precision, recall, F1, false positives per callable gigabase, metrics for truth VAF ≥0.25, SBS96 cosine similarity, false-positive spectrum, seed stability, and counts through every filtering stage. Documentation states that SEQC2 is a technical benchmark, not biological validation for organoids.

## Safe implementation sequence

1. Import tracked upstream source and record provenance.
2. Establish the GRCh38 WGS default profile and reference preflight.
3. Replace the sample model and validation; remove GRCh38 and CRAM restrictions.
4. Make untrimmed FASTQ and CRAM paths converge on validated alignments.
5. Adapt paired Mutect2/Strelka calls and prohibit automatic baseline PoN generation.
6. Add normalized caller tiers with caller-support annotations.
7. Add exact-population and recurrence annotation with reason-coded outputs.
8. Add cross-sample allele recounting.
9. Add role-aware QC and four SBS96 catalog sources.
10. Add synthetic end-to-end benchmark and explicit SEQC2 support.
11. Migrate Python dependency ownership to uv only if the resulting lock can be verified cleanly.
12. Update documentation and run verification in increasing cost order.

## Verification and acceptance

Required automated coverage includes manifest and baseline-sharing validation, reference and CRAM compatibility, FASTQ/CRAM resolution, exact-allele normalization, caller-support merge, population AF, recurrence, allele evidence filters, GRCh38 WGS configuration, benchmark metrics, absence of trimming from the active DAG, Snakemake dry run, and a miniature end-to-end test.

Verification order:

1. targeted Python tests;
2. configured lint/format checks;
3. full unit suite;
4. Snakemake dry run using generated tiny resources;
5. miniature integration test;
6. full SEQC2 benchmark only with explicit approval and adequate capacity.

Completion requires the functional and documentation acceptance criteria in the source specification. Any item not implemented or not executable in the local environment must be called out with the exact next command.

## Material open questions and chosen defaults

- Later-sample VAF default: use 0.20, configurable. This favors sensitivity at 6× while retaining a stringent caller-intersection tier.
- Baseline alternate evidence: annotate counts and use a configurable threshold/reason rather than treating zero alternate reads as proof of absence.
- CRAM conversion: prefer direct CRAM consumption when verified; create a temporary BAM only for a tool path that demonstrably requires it.
- Recounting implementation: choose a well-tested exact-allele pileup method after checking available tool semantics; document exclusions for overlapping pairs, base quality, and mapping quality.
- Repository strategy: selective fork is chosen as above. General-purpose fixes can later be proposed back to the upstream pipeline as isolated commits.

## Private site deployment runbook

The public repository must remain site-neutral. It documents only the generic path from a fresh SLURM checkout to a scheduler dry run and submission, without naming a private cluster, login host, account, filesystem, partition, allocation, or local data layout.

Detailed site deployment material belongs under the ignored `scripts/local/` tree and is copied manually rather than committed. The private runbook and overlays must:

- use the Git remote as the source of code and keep reference, sample, image, cache, run, and log data outside Git;
- capture the local login/transfer expectations, project allocation, filesystem placement, quota checks, modules, Pixi installation, locked environment creation, and Apptainer setup;
- provide an ignored `scripts/setup.sh`, `config/config.local.yaml`, and `config/samples.yaml` source material without embedding credentials;
- use absolute paths visible from compute nodes and avoid hardcoding unverified user, account, partition, or storage values;
- verify the checkout, environment, inputs, images, manifest, DAG, driver preview, submission, state, logs, restart, and retrieval in increasing cost order;
- provide a small one-comparison territory test and clearly distinguish it from the non-executable placeholder fixture;
- document safe update and recovery without deleting previous results.

Acceptance requires the public README and tracked examples to remain site-neutral, ignored private deployment files to be excluded by `git status`, shell syntax checks for changed scripts, the Python test suite, and a complete synthetic DAG dry run. The user executes remote commands interactively; authentication, allocation, data, or reference gaps must be reported rather than guessed.

The first private real-data scenario uses one existing high-depth organoid–baseline alignment pair. Preparation must validate and index the completed source BAMs, measure source depth over the configured full-WGS territory, estimate output storage, deterministically downsample both members to 6× while preserving read-pair sampling decisions, measure achieved depth, validate/index the resulting CRAMs, and generate a manifest using the actual alignment-header sample names. The analysis runs the normal default target so the deliverables include caller-specific and stringent variant calls plus the Mutect2, Strelka2, intersection, and stringent SBS96 context matrix. Source and generated BAM/CRAM data must be globally ignored by Git.

Reference provisioning acceptance requires unit coverage of plan/execute gating, capacity rejection, checksum rejection, restart behavior, configuration rendering, and manifest generation without network access. Documentation must keep reference provisioning separate from environment installation and show the operator the explicit plan and execute commands.

GRCh37 WGS acceptance additionally requires a synthetic native-build configuration test and a complete Snakemake dry run through the stringent catalog and SBS96 target. Tests must reject a GRCh37/GRCh38 alignment-reference mismatch before downsampling. The public documentation must explain that build selection follows the alignment, not operator preference, and must give explicit provisioning, configuration, dry-run, and execution commands for both supported WGS builds.
