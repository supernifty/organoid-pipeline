rule samtools_alignment_qc:
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        reference=config["reference"]["genome"],
        validation="results/qc/preflight/{sample}.alignment.json"
    output:
        flagstat="results/qc/samtools/{sample}.flagstat.txt",
        stats="results/qc/samtools/{sample}.stats.txt"
    log: "logs/qc/{sample}.samtools.log"
    benchmark: "results/benchmarks/qc/{sample}.samtools.tsv"
    params:
        flagstat=tmp_path("qc", "{sample}.flagstat.txt"),
        stats=tmp_path("qc", "{sample}.stats.txt")
    threads: 4
    resources: mem_mb=8192, runtime=240, disk_mb=1024
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {params.flagstat}) results/qc/samtools $(dirname {log})
        samtools flagstat -@ {threads} {input.cram} > {params.flagstat} 2> {log}
        samtools stats -@ {threads} --reference {input.reference} {input.cram} > {params.stats} 2>> {log}
        test -s {params.flagstat}; test -s {params.stats}
        mv {params.flagstat} {output.flagstat}; mv {params.stats} {output.stats}
        """

rule alignment_summary:
    output:
        "results/qc/metrics/{sample}.alignment_metrics.txt"
    input:
        ref=config["reference"]["genome"],
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample)
    benchmark: "results/benchmarks/qc/{sample}.alignment_summary.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        ref_container=container_path(config["reference"]["genome"]),
        java_options=gatk_java_options("14g"),
        cram_container=lambda wildcards: container_path(aligned_cram(wildcards.sample)),
        tmp=tmp_path("{sample}.alignment_metrics.txt")
    threads: 2
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' CollectAlignmentSummaryMetrics \
            R={params.ref_container} \
            I={params.cram_container} \
            O=/data/{params.tmp}
        test -s {params.tmp}
        mkdir -p results/qc/metrics
        mv {params.tmp} {output}
        """

rule insert_size:
    output:
        "results/qc/metrics/{sample}.insert_size.txt"
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        ref=config["reference"]["genome"]
    benchmark: "results/benchmarks/qc/{sample}.insert_size.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        ref_container=container_path(config["reference"]["genome"]),
        cram_container=lambda wildcards: container_path(aligned_cram(wildcards.sample)),
        tmp=tmp_path("{sample}.insert_size.txt"),
        tmp_hist=tmp_path("{sample}.insert_size.histogram.pdf")
    threads: 2
    shell:
        """
        {params.gatk_cmd} \
            gatk CollectInsertSizeMetrics \
            I={params.cram_container} \
            R={params.ref_container} \
            O=/data/{params.tmp} \
            H=/data/{params.tmp_hist}
        test -s {params.tmp}
        test -s {params.tmp_hist}
        mkdir -p results/qc/metrics
        mv {params.tmp} {output}
        mv {params.tmp_hist} {output}.histogram.pdf
        """

rule final_variant_counts_qc:
    output:
        "results/qc/variants/final_variant_counts_mqc.tsv"
    input:
        vcfs=variant_outputs()
    benchmark: "results/benchmarks/aggregate/final_variant_counts.tsv"
    params:
        inputs=lambda wildcards: " ".join(
            f"--input {tumour}=results/catalogs/{tumour}.stringent.vcf.gz"
            for tumour in tumour_samples()
        ),
        tmp=tmp_path("final_variant_counts_mqc.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/final_variant_counts.py \
            {params.inputs} \
            --output {params.tmp}
        test -s {params.tmp}
        mkdir -p results/qc/variants
        mv {params.tmp} {output}
        """

rule wgs_per_contig_coverage:
    output:
        mqc="results/qc/coverage/{sample}.wgs_coverage_mqc.tsv",
        regions=temp(tmp_path("coverage", "wgs", "{sample}.regions.bed.gz")),
        thresholds=temp(tmp_path("coverage", "wgs", "{sample}.thresholds.bed.gz"))
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        reference=config["reference"]["genome"],
        territory=lambda wildcards: analysis_territory(),
        territory_tbi=lambda wildcards: analysis_territory_index(),
        manifest=rules.analysis_manifest.output
    benchmark: "results/benchmarks/qc/{sample}.wgs_coverage.tsv"
    params:
        prefix=tmp_path("coverage", "wgs", "{sample}"),
        tmp=tmp_path("{sample}.wgs_coverage_mqc.tsv"),
        role=lambda wildcards: config["samples"]["samples"][wildcards.sample]["role"],
        expected=lambda wildcards: config["coverage"].get(
            "expected_organoid_depth" if config["samples"]["samples"][wildcards.sample]["role"] == "organoid" else "expected_baseline_depth"
        ),
        expected_arg=lambda wildcards: (
            f"--expected-depth {config['coverage'].get('expected_organoid_depth' if config['samples']['samples'][wildcards.sample]['role'] == 'organoid' else 'expected_baseline_depth')}"
            if config["coverage"].get("expected_organoid_depth" if config["samples"]["samples"][wildcards.sample]["role"] == "organoid" else "expected_baseline_depth") is not None else ""
        )
    threads: 4
    resources:
        disk_mb=lambda wildcards, input: max(4096, int(input.size_mb * 0.1))
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        mosdepth --threads {threads} --no-per-base --by {input.territory} \
            --thresholds 1,3,5,10,20,50 --fasta {input.reference} {params.prefix} {input.cram}
        python3 {PIPELINE_DIR}/workflow/scripts/wgs_coverage_mqc.py \
            --sample {wildcards.sample} --role {params.role} {params.expected_arg} \
            --warning-fraction {config[coverage][depth_warning_fraction]} \
            --regions {params.prefix}.regions.bed.gz --thresholds {params.prefix}.thresholds.bed.gz \
            --output {params.tmp}
        test -s {params.tmp}
        test -s {output.regions}
        test -s {output.thresholds}
        mkdir -p results/qc/coverage
        mv {params.tmp} {output.mqc}
        """


rule prepare_exon_coverage_bed:
    output:
        temp(tmp_path("coverage", "exons.bed"))
    input:
        lambda wildcards: config["coverage"]["exon_bed"]["path"]
    benchmark: "results/benchmarks/qc/prepare_exon_coverage_bed.tsv"
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/coverage_qc.py prepare-bed --source {input} --output {output}
        """


rule exon_gene_coverage:
    output:
        exon="results/qc/coverage/{sample}.exon_coverage.tsv",
        gene="results/qc/coverage/{sample}.gene_coverage.tsv",
        regions=temp(tmp_path("coverage", "exon", "{sample}.regions.bed.gz")),
        thresholds=temp(tmp_path("coverage", "exon", "{sample}.thresholds.bed.gz"))
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        reference=config["reference"]["genome"],
        bed=rules.prepare_exon_coverage_bed.output
    benchmark: "results/benchmarks/qc/{sample}.exon_coverage.tsv"
    params:
        prefix=tmp_path("coverage", "exon", "{sample}"),
        exon_tmp=tmp_path("coverage", "exon", "{sample}.exon.tsv"),
        gene_tmp=tmp_path("coverage", "exon", "{sample}.gene.tsv"),
        warning_depth=lambda wildcards: config["coverage"][
            "tumour_complete_coverage_depth" if wildcards.sample in config["samples"].get("tumours", {}) else "normal_complete_coverage_depth"
        ]
    threads: 4
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        mosdepth --threads {threads} --no-per-base --by {input.bed} \
            --thresholds 10,20,50,100 --fasta {input.reference} {params.prefix} {input.cram}
        python3 {PIPELINE_DIR}/workflow/scripts/coverage_qc.py summarize \
            --sample {wildcards.sample} --regions {params.prefix}.regions.bed.gz \
            --thresholds {params.prefix}.thresholds.bed.gz --warning-depth {params.warning_depth} \
            --exon-output {params.exon_tmp} --gene-output {params.gene_tmp}
        mkdir -p results/qc/coverage
        mv {params.exon_tmp} {output.exon}
        mv {params.gene_tmp} {output.gene}
        mv {params.prefix}.regions.bed.gz {output.regions}
        mv {params.prefix}.thresholds.bed.gz {output.thresholds}
        """

rule multiqc_report:
    output:
        "results/aggregate/qc_summary.html"
    input:
        dups=expand("results/qc/metrics/{sample}.dups.metrics.txt", sample=pipeline_aligned_samples()),
        alignment=expand("results/qc/metrics/{sample}.alignment_metrics.txt", sample=config["samples"]["samples"]),
        insert_size=expand("results/qc/metrics/{sample}.insert_size.txt", sample=config["samples"]["samples"]),
        fastqc=fastqc_outputs(),
        samtools=expand("results/qc/samtools/{sample}.{kind}.txt", sample=config["samples"]["samples"], kind=["flagstat", "stats"]),
        final_variant_counts=rules.final_variant_counts_qc.output,
        wgs_coverage=(
            expand("results/qc/coverage/{sample}.wgs_coverage_mqc.tsv", sample=config["samples"]["samples"])
            if ANALYSIS_TYPE == "wgs" else []
        ),
        exon_coverage=(
            expand("results/qc/coverage/{sample}.{kind}_coverage.tsv", sample=config["samples"]["samples"], kind=["exon", "gene"])
            if config.get("coverage", {}).get("exon_enabled", False) else []
        ),
        somalier=somalier_outputs()
    benchmark: "results/benchmarks/aggregate/multiqc.tsv"
    params:
        outdir=tmp_path("multiqc"),
        html=tmp_path("multiqc", "qc_summary.html")
    threads: 1
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        multiqc --force --no-version-check -o {params.outdir} -n qc_summary.html results/qc/
        test -s {params.html}
        mkdir -p results/aggregate
        mv {params.html} {output}
        rm -rf {params.outdir}
        """
