rule validate_alignment_input:
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        reference=config["reference"]["genome"],
        fai=lambda wildcards: f"{config['reference']['genome']}.fai"
    output: "results/qc/preflight/{sample}.alignment.json"
    log: "logs/preflight/{sample}.alignment.log"
    benchmark: "results/benchmarks/preflight/{sample}.alignment.tsv"
    params:
        sample=lambda wildcards: bam_sample_name(wildcards.sample),
        contigs=" ".join(f"--contig {contig}" for contig in ANALYSIS["contigs"]),
        temporary=tmp_path("preflight", "{sample}.alignment.json")
    threads: 1
    resources: mem_mb=1024, runtime=60, disk_mb=1024
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {params.temporary}) $(dirname {log})
        python3 {PIPELINE_DIR}/workflow/scripts/validate_alignment.py --cram {input.cram} --crai {input.crai} \
          --reference {input.reference} --fai {input.fai} --sample {params.sample} {params.contigs} \
          --output {params.temporary} > {log} 2>&1
        test -s {params.temporary}; mkdir -p results/qc/preflight; mv {params.temporary} {output}
        """
