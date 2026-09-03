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

## Goals

- Give `samtools_alignment_qc` enough conservative headroom to resume and
  complete the observed low-depth WGS workload.
- Ensure FastQC and GATK jobs do not set a Java heap equal to their complete
  scheduler memory request.
- Keep in-rule resources and the native SLURM profile consistent where both
  specify memory.
- Add automated checks that prevent these unsafe allocations from returning.

## Scope

Increase `samtools_alignment_qc` from 2,048 MB to 8,192 MB because the failed
run cannot reveal how far memory would have grown past the enforced 2 GiB
limit. Increase FastQC from 2,048 MB to 4,096 MB; its four-thread invocation
configures a 2,048 MB Java heap. Increase `split_wgs_intervals`,
`contamination_sites`, and `mutect2_orientation_model` to 8,192 MB because
each configures a 4 GiB GATK heap. Apply the values in rule declarations and
the SLURM profile as applicable.

Other rules are unchanged in this pass. Heavy native tools already have larger
allocations, and other fixed-heap GATK rules have scheduler headroom. Defaults
for lightweight rules remain 4,096 MB. Successful benchmark and scheduler
measurements will be used for later right-sizing rather than speculative
increases.

## Acceptance criteria

- `samtools_alignment_qc` requests 8,192 MB in both its rule and the SLURM
  profile.
- FastQC requests 4,096 MB in both its rule and the SLURM profile while its
  current four-thread invocation retains a 2,048 MB maximum heap.
- `split_wgs_intervals`, `contamination_sites`, and
  `mutect2_orientation_model` request 8,192 MB in both their rule declarations
  and the SLURM profile.
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
