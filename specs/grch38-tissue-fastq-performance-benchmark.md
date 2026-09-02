# GRCh38 tissue FASTQ performance benchmark

## Problem

The pipeline's principal use case is approximately 6× WGS, but its end-to-end
performance starting from genuine paired FASTQs has not yet been measured on the
target SLURM environment. Existing organoid pilot inputs are pre-aligned, so they
do not exercise FastQC, BWA-MEM, coordinate sorting, or duplicate marking.

A matched tissue WGS pair is available as one R1/R2 FASTQ pair per biological
sample. It is suitable for an operational GRCh38 FASTQ-path benchmark after
deterministic, read-pair-preserving downsampling to approximately 6×. It is not
an early-passage organoid comparison and has no somatic truth set.

## Goals

- Exercise the complete GRCh38 WGS pipeline from genuine paired FASTQs through
  FastQC, alignment, CRAM 3.0 creation, duplicate marking, Mutect2, Strelka2,
  caller tiers, cohort recounting, filtering, SBS96, mosdepth QC, MultiQC, and
  provenance.
- Downsample both members of the matched pair reproducibly without trimming or
  altering retained sequences and qualities.
- Calibrate the prepared inputs against measured mean autosomal depth rather
  than treating nominal sequenced bases as achieved coverage.
- Measure experienced turnaround separately from per-job compute use, queueing,
  memory, and I/O so later runs can be planned realistically.
- Preserve enough input, configuration, tool, container, and run provenance to
  reproduce and interpret the benchmark.

## Non-goals

- Do not use the tissue result to validate monoclonal organoid biology or the
  early-passage-baseline interpretation.
- Do not claim caller precision, recall, or biological accuracy without a truth
  set.
- Do not compare the tissue calls directly with the organoid CRAM-input calls.
- Do not benchmark full-depth tissue FASTQs in this phase.
- Do not derive FASTQs from the organoid BAMs or CRAMs in this phase.
- Do not enable optional VEP annotation, germline calling, Somalier, or known-
  signature fitting merely for this timing benchmark.

## Inputs and configuration

The benchmark uses:

- one case R1/R2 FASTQ pair;
- one matched-normal R1/R2 FASTQ pair;
- GRCh38 as the reference build;
- the pinned, validated GRCh38 resource bundle used by the production profile;
- explicit donor, condition, role, lineage, library, platform, and unit metadata;
- independent, recorded sampling seeds for case and normal;
- an initial target depth of 6× for both samples;
- Apptainer containers and the existing native Snakemake SLURM executor;
- `cram_version: "3.0"` for GATK 4.4/HTSJDK compatibility.

The public implementation and documentation must remain site-neutral. Exact
sample identifiers, cluster paths, accounts, partitions, and local module stacks
belong only in the ignored private scenario under `scripts/local/`.

Although the manifest must label the samples `organoid` and `baseline` to use
the paired somatic interface, the private scenario and provenance must clearly
identify them as a tissue case and matched tissue normal. The output must not be
described as mutations acquired after an ancestral organoid baseline.

## Deterministic paired-FASTQ downsampling

Add a site-neutral preparation command with a non-mutating plan mode and an
explicit execution flag. It must accept R1, R2, output paths, target depth,
calling territory, seed, threads, and report path.

Before execution it must:

1. Confirm that both inputs are readable gzip-compressed FASTQs.
2. Stream both inputs in lockstep and reject truncated records, malformed
   records, mismatched normalized read names, or unequal pair counts.
3. Count input read pairs and sequenced bases without loading the dataset into
   memory.
4. Calculate callable-territory bases from the configured, validated GRCh38 WGS
   territory.
5. Calculate the initial sampling fraction from target sequenced bases divided
   by observed input sequenced bases and reject a source that cannot supply the
   requested nominal depth.
6. Estimate output size and verify available capacity with a safety margin
   before creating task outputs.
7. Print the inputs, identities, seed, target, fraction, estimated outputs, and
   capacity without writing outputs unless execution was explicitly requested.

During execution, normalize `/1` and `/2` mate suffixes only for pair identity
comparison. Select or reject each pair from a stable cryptographic hash of the
canonical read name and seed. Both mates must always receive the same decision.
Retained FASTQ records must otherwise be byte-equivalent in identifier,
sequence, separator, and quality content. Adapter and quality trimming are
forbidden.

Write compressed R1 and R2 through task-specific temporary files, validate both
streams, and move them atomically to final paths only after success. Compute
output SHA-256 checksums while writing or validating. Interrupted or invalid
outputs must not appear complete.

The JSON report must contain at least:

- input absolute paths, sizes, and modification times;
- output absolute paths, sizes, and SHA-256 checksums;
- normalized sample label;
- seed and selection algorithm/version;
- target depth and callable-territory bases;
- input/output pair counts and sequenced bases;
- calculated sampling fraction and nominal output depth;
- tool version, command, start/end timestamps, and completion status.

Restart reuse is permitted only when both outputs validate and the report
matches the current input identities, seed, target, territory identity,
selection algorithm, and output checksums. Parameter or input changes must
cause atomic regeneration.

## Generated manifest and private scenario

Create an ignored private wrapper that invokes the generic downsampler for both
samples, checks aggregate capacity, and writes a generated sample manifest only
after both outputs validate. The manifest must contain one R1/R2 pair per
sample, complete read-group fields, consistent donor/lineage metadata, and one
case-to-normal comparison.

The wrapper must be restart-safe and must not modify or remove the source
FASTQs. Generated FASTQs, reports, manifests, and benchmark results must remain
ignored by Git. The private runbook must provide plan, execution, manifest
validation, DAG dry-run, calibration, full-run, monitoring, and reporting
commands.

## Coverage calibration

Nominal sequenced depth is not equivalent to mapped, duplicate-filtered
autosomal depth. Before launching callers:

1. Start an isolated run-manager batch targeting the two WGS mosdepth summary
   outputs only. This exercises FastQC, BWA-MEM, coordinate CRAM sorting,
   MarkDuplicates, preflight, and coverage QC without running callers.
2. Treat pipeline mosdepth mean autosomal depth as authoritative.
3. Accept each sample when achieved depth is from 5.4× through 6.6×,
   inclusive.
4. If either sample falls outside that range, calculate
   `new_fraction = old_fraction × 6 / observed_depth`, regenerate that sample's
   FASTQs atomically from the original source, and rerun calibration.
5. Reject an adjusted fraction greater than 1 rather than fabricating coverage.
6. Record every attempted fraction and achieved depth; do not overwrite the
   audit history.

Once both samples pass, resume the same batch with the default `rule all`
target. The validated alignment and QC products should be reused, while the
remaining callers, catalogs, signatures, aggregate QC, and provenance run
normally.

## Performance reporting

Add a site-neutral post-run summarizer that reads the run-manager record,
effective configuration, analysis manifest, and all Snakemake benchmark TSVs.
It must emit machine-readable TSV/JSON plus a compact human-readable report.

Report at least:

- Git revision and dirty-state provenance;
- reference build and identity, callable bases, achieved depths, FASTQ pair
  counts, read lengths, and compressed input sizes;
- container identities, tool versions, shard count, Mutect2 concurrency, and
  relevant rule resource requests;
- controller submission, start, completion, and queue-inclusive turnaround;
- per-rule and per-sample wall time, CPU time, maximum RSS, input/output I/O,
  completion state, and benchmark path;
- grouped costs for preparation, FastQC, BWA/sort, MarkDuplicates, alignment
  QC, Mutect2, Strelka2, recounting/filtering, SBS96, and aggregate QC;
- total CPU consumption and maximum observed memory;
- calibration attempts and the distinction between calibration time and the
  remaining full-DAG time;
- final output sizes and available storage where measurable.

Parallel job wall times must not be summed and presented as elapsed turnaround.
Present separately:

- experienced controller elapsed time;
- summed job wall time;
- summed CPU time;
- scheduler waiting where it can be established reliably.

If a metric cannot be derived reliably, report it as unavailable rather than
guessing. The report must state that results are cluster-, queue-, filesystem-,
sample-, and concurrency-dependent.

## Failure handling

- Malformed or unsynchronized FASTQs fail before output publication.
- Insufficient disk capacity fails before downsampling begins.
- Partial FASTQs and reports are never accepted as restart-complete.
- A generated alignment using a CRAM version other than configured 3.0 fails in
  production preflight before GATK/Picard work.
- A failed calibration or full run remains in its isolated run-manager batch
  with immutable configuration and launch history.
- Changes to FASTQ identity or sampling parameters invalidate reusable
  alignment outputs through existing run-manager fingerprints.
- Full-pipeline execution must not begin until both samples pass the achieved-
  depth gate.

## Verification

Verify in increasing cost order:

1. Unit tests for FASTQ parsing, canonical mate names, deterministic pair
   selection, different seeds, malformed/truncated input, unequal pairs,
   insufficient source depth, capacity rejection, report identity, checksums,
   atomic publication, and restart invalidation.
2. Unit tests for performance aggregation, including parallel wall-time versus
   elapsed-time handling and missing metrics.
3. Lint and formatting checks through existing Pixi tasks.
4. The complete Python test suite.
5. A synthetic GRCh38 FASTQ DAG dry run proving both samples traverse FastQC,
   BWA-MEM, sorting, MarkDuplicates, CRAM 3.0 preflight, mosdepth, both callers,
   filtering, SBS96, and MultiQC without trimming.
6. A tiny executable paired-FASTQ integration fixture, if supported by the
   existing integration-test conventions.
7. Private plan-mode downsampling and capacity review.
8. Private calibration DAG execution and achieved-depth review.
9. The full default tissue benchmark and post-run performance report.

## Acceptance criteria

The benchmark is complete when:

- the original tissue FASTQs remain unchanged;
- both generated FASTQ pairs are reproducible, pair-synchronized, checksummed,
  untrimmed, and provenance-recorded;
- both samples achieve 5.4–6.6× mean autosomal depth;
- the generated manifest passes validation and the full GRCh38 FASTQ DAG dry
  run includes every required alignment, calling, QC, filtering, and signature
  stage;
- the isolated default pipeline batch completes successfully;
- output includes Mutect2, Strelka2, union/intersection/caller-only tiers,
  reason-coded catalogs, SBS96 counts, mosdepth QC, MultiQC, and provenance;
- the performance report distinguishes elapsed, queued, wall, and CPU time and
  includes peak memory and storage observations;
- no site-specific identifiers or large data/results are committed; and
- documentation explicitly limits conclusions to operational performance on a
  matched GRCh38 tissue pair at approximately 6×.
