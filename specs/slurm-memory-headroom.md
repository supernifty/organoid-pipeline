# SLURM memory headroom for production rules

## Problem

The first end-to-end 6x paired-FASTQ WGS run demonstrated that
`samtools_alignment_qc` is under-provisioned. Spartan job `29966603` requested
2,048 MB, reached 1,909.94 MB RSS, and its `python3.13` job step was terminated
with `OUT_OF_MEMORY` after 3 minutes 30 seconds. The rule runs four-threaded
`samtools flagstat` and `samtools stats` over CRAM input, so its observed memory
at termination is only a lower bound on successful peak use.

The resource profile also contains rules whose configured Java maximum heap is
equal to the entire SLURM memory allocation. This leaves no space for JVM
non-heap memory, native libraries, the container runtime, or wrapper processes
and creates avoidable OOM risk.

The resumed run then exposed the same issue in the in-memory Python candidate
processing chain. Spartan job `30057077` ran `caller_tiers` with 2,048 MB,
reached 1,910.21 MB RSS, and was terminated with `OUT_OF_MEMORY` after 8
seconds. The script retains complete normalized caller VCF dictionaries plus
derived sets and output mappings. Its downstream cohort union, allele recount,
catalog filtering, and SBS96 steps also retain whole candidate collections or
captured subprocess output, so their previous 2–4 GiB requests are unsafe for
the same real WGS candidate volume.

The completed paired-FASTQ run subsequently measured successful peak RSS of
2,555.61 MB for `caller_tiers`, 2,994.10 MB for `cohort_allele_recount`, and
123.6 MB for `filter_organoid_catalog`. The production run is expected to have
20 samples. Cohort allele recount therefore retains its conservative 65,536 MB
allocation because both the candidate union and per-sample pileup work scale
with cohort size. Caller tiers is per tumour, while catalog filtering retains
enough headroom for the expected cohort after scaling the observed workload.

## Goals

- Give `samtools_alignment_qc` enough conservative headroom to resume and
  complete the observed low-depth WGS workload.
- Ensure FastQC and GATK jobs do not set a Java heap equal to their complete
  scheduler memory request.
- Keep in-rule resources and the native SLURM profile consistent where both
  specify memory.
- Provision the active in-memory WGS candidate-processing chain conservatively
  so the operational run does not discover each low limit serially.
- Add automated checks that prevent these unsafe allocations from returning.

## Scope

Increase `samtools_alignment_qc` from 2,048 MB to 8,192 MB because the failed
run cannot reveal how far memory would have grown past the enforced 2 GiB
limit. Increase FastQC from 2,048 MB to 4,096 MB; its four-thread invocation
configures a 2,048 MB Java heap. Increase `split_wgs_intervals`,
`contamination_sites`, and `mutect2_orientation_model` to 8,192 MB because
each configures a 4 GiB GATK heap. Apply the values in rule declarations and
the SLURM profile as applicable.

Other rules outside the candidate-processing chain are unchanged in this pass.
Heavy native tools already have larger allocations, and other fixed-heap GATK
rules have scheduler headroom. Defaults for lightweight rules remain 4,096 MB.

For the active candidate-processing chain, set `caller_tiers` to 16,384 MB and
retain `cohort_candidate_union` at 32,768 MB. Retain
`cohort_allele_recount` at 65,536 MB because it captures complete
`samtools mpileup` output and builds result rows in memory across the expected
20-sample cohort. Set `filter_organoid_catalog` to 8,192 MB; this remains more than 60
times the one-pair observed peak and leaves headroom for its cohort VCF
dictionaries, recount rows, annotations, and output classifications. Retain
`sbs96_catalogs` at 32,768 MB pending benchmark evidence from the corrected
run.

## Acceptance criteria

- `samtools_alignment_qc` requests 8,192 MB in both its rule and the SLURM
  profile.
- FastQC requests 4,096 MB in both its rule and the SLURM profile while its
  current four-thread invocation retains a 2,048 MB maximum heap.
- `split_wgs_intervals`, `contamination_sites`, and
  `mutect2_orientation_model` request 8,192 MB in both their rule declarations
  and the SLURM profile.
- `caller_tiers` requests 16,384 MB in both its rule declaration and the SLURM
  profile.
- `cohort_candidate_union` and `sbs96_catalogs` continue to request 32,768 MB
  in both their rule declarations and the SLURM profile.
- `cohort_allele_recount` continues to request 65,536 MB in both its rule
  declaration and the SLURM profile for the expected 20-sample cohort.
- `filter_organoid_catalog` requests 8,192 MB in both its rule declaration and
  the SLURM profile.
- Tests assert these production memory requests and the intended minimum
  scheduler-to-heap headroom.
- Existing outputs remain reusable when an inactive failed batch is explicitly
  resumed after updating to the new revision.

## Verification

1. Run focused resource-profile tests.
2. Run formatting and lint checks.
3. Run the complete Python test suite.
4. Run the repository's synthetic pipeline verification, including its
   Snakemake DAG checks.
5. Inspect the final diff and confirm unrelated working-tree files are not
   included in the commit.
