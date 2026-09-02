"""Matched-normal germline calling with GATK HaplotypeCaller."""


GERMLINE_FILTER_DEFAULTS = {
    "snp": {
        "qd_min": 2.0,
        "qual_min": 30.0,
        "sor_max": 3.0,
        "fs_max": 60.0,
        "mq_min": 40.0,
        "mq_rank_sum_min": -12.5,
        "read_pos_rank_sum_min": -8.0,
    },
    "indel": {
        "qd_min": 2.0,
        "qual_min": 30.0,
        "fs_max": 200.0,
        "read_pos_rank_sum_min": -20.0,
    },
}

GERMLINE_CONFIG = config.get("germline", {})
if not isinstance(GERMLINE_CONFIG.get("enabled", True), bool):
    raise ValueError("germline.enabled must be true or false")
GERMLINE_MAX_CONCURRENT_SHARDS = GERMLINE_CONFIG.get(
    "max_concurrent_haplotypecaller_shards", 16
)
if (
    isinstance(GERMLINE_MAX_CONCURRENT_SHARDS, bool)
    or not isinstance(GERMLINE_MAX_CONCURRENT_SHARDS, int)
    or GERMLINE_MAX_CONCURRENT_SHARDS <= 0
):
    raise ValueError("germline.max_concurrent_haplotypecaller_shards must be a positive integer")


def germline_filter_values(variant_type):
    configured = config.get("germline", {}).get("hard_filters", {}).get(variant_type, {})
    if not isinstance(configured, dict):
        raise ValueError(f"germline.hard_filters.{variant_type} must be a mapping")
    unknown = sorted(set(configured) - set(GERMLINE_FILTER_DEFAULTS[variant_type]))
    if unknown:
        raise ValueError(
            f"Unknown germline.hard_filters.{variant_type} keys: {', '.join(unknown)}"
        )
    values = {**GERMLINE_FILTER_DEFAULTS[variant_type], **configured}
    if any(isinstance(value, bool) for value in values.values()):
        raise ValueError(f"germline.hard_filters.{variant_type} values must be numeric")
    try:
        return {key: float(value) for key, value in values.items()}
    except (TypeError, ValueError) as error:
        raise ValueError(f"germline.hard_filters.{variant_type} values must be numeric") from error


def germline_variant_filter_args(variant_type):
    values = germline_filter_values(variant_type)
    if variant_type == "snp":
        expressions = (
            ("SNP_QD", f"QD < {values['qd_min']}"),
            ("SNP_QUAL", f"QUAL < {values['qual_min']}"),
            ("SNP_SOR", f"SOR > {values['sor_max']}"),
            ("SNP_FS", f"FS > {values['fs_max']}"),
            ("SNP_MQ", f"MQ < {values['mq_min']}"),
            ("SNP_MQRankSum", f"MQRankSum < {values['mq_rank_sum_min']}"),
            ("SNP_ReadPosRankSum", f"ReadPosRankSum < {values['read_pos_rank_sum_min']}"),
        )
    else:
        expressions = (
            ("INDEL_QD", f"QD < {values['qd_min']}"),
            ("INDEL_QUAL", f"QUAL < {values['qual_min']}"),
            ("INDEL_FS", f"FS > {values['fs_max']}"),
            ("INDEL_ReadPosRankSum", f"ReadPosRankSum < {values['read_pos_rank_sum_min']}"),
        )
    return " ".join(
        f"--filter-name {name} --filter-expression '{expression}'"
        for name, expression in expressions
    )


rule germline_haplotypecaller_shard:
    output:
        vcf=temp(tmp_path("germline", "{normal}.{chromosome}.haplotypecaller.g.vcf.gz")),
        tbi=temp(tmp_path("germline", "{normal}.{chromosome}.haplotypecaller.g.vcf.gz.tbi"))
    input:
        cram=lambda wildcards: aligned_cram(wildcards.normal),
        crai=lambda wildcards: aligned_cram_index(wildcards.normal),
        reference=config["reference"]["genome"],
        shard_dependency=shard_dependency,
        manifest=rules.analysis_manifest.output
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("12g"),
        reference=lambda wildcards: container_path(config["reference"]["genome"]),
        cram=lambda wildcards: container_path(aligned_cram(wildcards.normal)),
        bam_sample=lambda wildcards: bam_sample_name(wildcards.normal),
        core=lambda wildcards: container_path(shard_core_interval(wildcards)),
        shard_contig_arg=shard_contig_arg,
        publish=lambda wildcards: tmp_path(
            "germline", f"{wildcards.normal}.{wildcards.chromosome}.haplotypecaller.publish.g.vcf.gz"
        ),
        publish_container=lambda wildcards: container_path(tmp_path(
            "germline", f"{wildcards.normal}.{wildcards.chromosome}.haplotypecaller.publish.g.vcf.gz"
        ))
    threads: 4
    resources:
        haplotypecaller_shards=1,
        mem_mb=16384,
        runtime=1440
    shell:
        """
        mkdir -p $(dirname {params.publish})
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' HaplotypeCaller \
            -R {params.reference} \
            -I {params.cram} \
            --sample-name {params.bam_sample} \
            -ERC GVCF \
            -L {params.core} \
            {params.shard_contig_arg} \
            --interval-padding 0 \
            -O {params.publish_container}

        gzip -t {params.publish}
        tabix -f -p vcf {params.publish}
        test -s {params.publish}.tbi
        mv {params.publish} {output.vcf}
        mv {params.publish}.tbi {output.tbi}
        """


rule germline_merge_gvcf:
    output:
        vcf="results/germline/{normal}.haplotypecaller.g.vcf.gz",
        tbi="results/germline/{normal}.haplotypecaller.g.vcf.gz.tbi"
    input:
        vcfs=lambda wildcards: [
            tmp_path("germline", f"{wildcards.normal}.{shard}.haplotypecaller.g.vcf.gz")
            for shard in SCATTER_IDS
        ],
        tbis=lambda wildcards: [
            tmp_path("germline", f"{wildcards.normal}.{shard}.haplotypecaller.g.vcf.gz.tbi")
            for shard in SCATTER_IDS
        ]
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        inputs=lambda wildcards: " ".join(
            f"-I {container_path(tmp_path('germline', f'{wildcards.normal}.{shard}.haplotypecaller.g.vcf.gz'))}"
            for shard in SCATTER_IDS
        ),
        publish=lambda wildcards: tmp_path("germline", f"{wildcards.normal}.haplotypecaller.publish.g.vcf.gz"),
        publish_container=lambda wildcards: container_path(
            tmp_path("germline", f"{wildcards.normal}.haplotypecaller.publish.g.vcf.gz")
        )
    threads: 2
    resources:
        mem_mb=8192,
        runtime=120
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' MergeVcfs \
            {params.inputs} \
            -O {params.publish_container}

        gzip -t {params.publish}
        tabix -f -p vcf {params.publish}
        test -s {params.publish}.tbi
        mkdir -p results/germline
        mv {params.publish} {output.vcf}
        mv {params.publish}.tbi {output.tbi}
        """


rule germline_genotype:
    output:
        vcf=temp(tmp_path("germline", "{normal}.haplotypecaller.genotyped.vcf.gz")),
        tbi=temp(tmp_path("germline", "{normal}.haplotypecaller.genotyped.vcf.gz.tbi"))
    input:
        gvcf="results/germline/{normal}.haplotypecaller.g.vcf.gz",
        gvcf_tbi="results/germline/{normal}.haplotypecaller.g.vcf.gz.tbi",
        reference=config["reference"]["genome"],
        manifest=rules.analysis_manifest.output
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("12g"),
        reference=lambda wildcards: container_path(config["reference"]["genome"]),
        gvcf=lambda wildcards: container_path(
            f"results/germline/{wildcards.normal}.haplotypecaller.g.vcf.gz"
        ),
        publish=lambda wildcards: tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.genotyped.publish.vcf.gz"
        ),
        publish_container=lambda wildcards: container_path(tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.genotyped.publish.vcf.gz"
        ))
    threads: 2
    resources:
        haplotypecaller_shards=1,
        mem_mb=16384,
        runtime=720
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' GenotypeGVCFs \
            -R {params.reference} \
            -V {params.gvcf} \
            -O {params.publish_container}

        gzip -t {params.publish}
        tabix -f -p vcf {params.publish}
        test -s {params.publish}.tbi
        mv {params.publish} {output.vcf}
        mv {params.publish}.tbi {output.tbi}
        """


rule germline_filter:
    output:
        vcf="results/germline/{normal}.haplotypecaller.filtered.vcf.gz",
        tbi="results/germline/{normal}.haplotypecaller.filtered.vcf.gz.tbi"
    input:
        vcf=tmp_path("germline", "{normal}.haplotypecaller.genotyped.vcf.gz"),
        tbi=tmp_path("germline", "{normal}.haplotypecaller.genotyped.vcf.gz.tbi"),
        reference=config["reference"]["genome"],
        manifest=rules.analysis_manifest.output
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        input_vcf=lambda wildcards: container_path(
            tmp_path("germline", f"{wildcards.normal}.haplotypecaller.genotyped.vcf.gz")
        ),
        snps=lambda wildcards: tmp_path("germline", f"{wildcards.normal}.haplotypecaller.snps.vcf.gz"),
        snps_container=lambda wildcards: container_path(
            tmp_path("germline", f"{wildcards.normal}.haplotypecaller.snps.vcf.gz")
        ),
        indels=lambda wildcards: tmp_path("germline", f"{wildcards.normal}.haplotypecaller.indels.vcf.gz"),
        indels_container=lambda wildcards: container_path(
            tmp_path("germline", f"{wildcards.normal}.haplotypecaller.indels.vcf.gz")
        ),
        filtered_snps=lambda wildcards: tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.snps.filtered.vcf.gz"
        ),
        filtered_snps_container=lambda wildcards: container_path(tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.snps.filtered.vcf.gz"
        )),
        filtered_indels=lambda wildcards: tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.indels.filtered.vcf.gz"
        ),
        filtered_indels_container=lambda wildcards: container_path(tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.indels.filtered.vcf.gz"
        )),
        snp_filters=germline_variant_filter_args("snp"),
        indel_filters=germline_variant_filter_args("indel"),
        publish=lambda wildcards: tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.filtered.publish.vcf.gz"
        ),
        publish_container=lambda wildcards: container_path(tmp_path(
            "germline", f"{wildcards.normal}.haplotypecaller.filtered.publish.vcf.gz"
        ))
    threads: 2
    resources:
        mem_mb=8192,
        runtime=240
    shell:
        """
        {params.gatk_cmd} gatk --java-options '{params.java_options}' SelectVariants \
            -V {params.input_vcf} \
            --exclude-non-variants true \
            --select-type-to-include SNP \
            -O {params.snps_container}

        {params.gatk_cmd} gatk --java-options '{params.java_options}' SelectVariants \
            -V {params.input_vcf} \
            --exclude-non-variants true \
            --select-type-to-include INDEL \
            --select-type-to-include MIXED \
            -O {params.indels_container}

        {params.gatk_cmd} gatk --java-options '{params.java_options}' VariantFiltration \
            -V {params.snps_container} \
            {params.snp_filters} \
            -O {params.filtered_snps_container}

        {params.gatk_cmd} gatk --java-options '{params.java_options}' VariantFiltration \
            -V {params.indels_container} \
            {params.indel_filters} \
            -O {params.filtered_indels_container}

        {params.gatk_cmd} gatk --java-options '{params.java_options}' MergeVcfs \
            -I {params.filtered_snps_container} \
            -I {params.filtered_indels_container} \
            -O {params.publish_container}

        gzip -t {params.publish}
        tabix -f -p vcf {params.publish}
        test -s {params.publish}.tbi
        mkdir -p results/germline
        mv {params.publish} {output.vcf}
        mv {params.publish}.tbi {output.tbi}
        rm -f \
            {params.snps} {params.snps}.tbi \
            {params.indels} {params.indels}.tbi \
            {params.filtered_snps} {params.filtered_snps}.tbi \
            {params.filtered_indels} {params.filtered_indels}.tbi
        """
