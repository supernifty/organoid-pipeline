def read_group(wildcards):
    sample = wildcards.sample
    values = config["samples"]["samples"][sample]["read_group"]
    return (
        f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{values['platform']}"
        f"\\tLB:{values['library']}\\tPU:{values['unit']}"
    )


def get_paired_fastq_inputs(wildcards):
    sample = wildcards.sample
    return {
        "ref": config["reference"]["genome"],
        "bwa_indexes": bwa_index_inputs(),
        "r1": config["samples"]["samples"][sample]["fastq_1"],
        "r2": config["samples"]["samples"][sample]["fastq_2"],
    }


def merged_alignment_input(wildcards):
    return tmp_path(f"{wildcards.sample}.paired.cram")


rule bwa_mem_paired:
    input:
        unpack(get_paired_fastq_inputs)
    output:
        cram=temp(tmp_path("{sample}.paired.cram"))
    log: "logs/alignment/{sample}.bwa.log"
    benchmark: "results/benchmarks/alignment/{sample}.bwa.tsv"
    params:
        rg=read_group,
        cram_version=lambda wildcards: config.get("cram_version", "3.0"),
        temporary=tmp_path("{sample}.paired.publish.cram")
    threads: 32
    resources:
        mem_mb=49152,
        runtime=1920,
        disk_mb=lambda wildcards, input: max(8192, int(input.size_mb * 2.5))
    shell:
        r"""
        mkdir -p $(dirname {log})
        bwa mem -M -t 24 -R "{params.rg}" {input.ref} {input.r1} {input.r2} 2> {log} | \
          samtools sort -@ 8 -T {params.temporary}.sort -O CRAM --reference {input.ref} \
            --output-fmt-option version={params.cram_version} -o {params.temporary} - 2>> {log}
        samtools quickcheck {params.temporary}
        test -s {params.temporary}
        mv {params.temporary} {output.cram}
        """


rule mark_duplicates:
    input:
        cram=merged_alignment_input,
        ref=config["reference"]["genome"]
    output:
        cram="results/cram/{sample}.sorted.dups.cram",
        crai="results/cram/{sample}.sorted.dups.cram.crai",
        metrics="results/qc/metrics/{sample}.dups.metrics.txt"
    log: "logs/alignment/{sample}.mark_duplicates.log"
    benchmark: "results/benchmarks/alignment/{sample}.mark_duplicates.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("44g"),
        ref_container=container_path(config["reference"]["genome"]),
        cram_container=lambda wildcards: container_path(merged_alignment_input(wildcards)),
        temporary_cram=tmp_path("{sample}.sorted.dups.cram"),
        temporary_metrics=tmp_path("{sample}.dups.metrics.txt"),
        picard_tmp=lambda wildcards: scratch_path(f"picard_{wildcards.sample}")
    threads: 1
    resources:
        mem_mb=49152,
        runtime=480,
        disk_mb=lambda wildcards, input: max(16384, int(input.size_mb * 2.5))
    shell:
        """
        mkdir -p $(dirname {log})
        {params.gatk_cmd} gatk --java-options '{params.java_options}' MarkDuplicates \
          -I {params.cram_container} -O /data/{params.temporary_cram} \
          -M /data/{params.temporary_metrics} -R {params.ref_container} \
          --VALIDATION_STRINGENCY LENIENT --TMP_DIR {params.picard_tmp} > {log} 2>&1
        samtools index {params.temporary_cram}
        samtools quickcheck {params.temporary_cram}
        test -s {params.temporary_cram}; test -s {params.temporary_cram}.crai
        test -s {params.temporary_metrics}
        mkdir -p results/cram results/qc/metrics
        mv {params.temporary_cram} {output.cram}
        mv {params.temporary_cram}.crai {output.crai}
        mv {params.temporary_metrics} {output.metrics}
        """
