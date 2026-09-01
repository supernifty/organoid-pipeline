def annotation_resource(name):
    return config["annotation"]["resources"][name]["path"]


def annotation_index_inputs(wildcards=None):
    return [
        value["index"]["path"]
        for value in config.get("annotation", {}).get("resources", {}).values()
        if isinstance(value, dict) and isinstance(value.get("index"), dict) and value["index"].get("path")
    ] if ANNOTATION_ENABLED else []


def annotation_input_paths(wildcards=None):
    return [value["path"] for value in config.get("annotation", {}).get("resources", {}).values() if value.get("path")] if ANNOTATION_ENABLED else []


def vep_plugin_args(wildcards=None):
    if not ANNOTATION_ENABLED:
        return ""
    return " ".join((
        f"--custom {container_path(annotation_resource('clinvar'))},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN",
        f"--custom {container_path(annotation_resource('gnomad'))},gnomAD,vcf,exact,0,AF,AF_popmax",
        f"--plugin REVEL,file={container_path(annotation_resource('revel'))}",
        f"--plugin CADD,snv={container_path(annotation_resource('cadd_snv'))},indels={container_path(annotation_resource('cadd_indel'))}",
        f"--plugin AlphaMissense,file={container_path(annotation_resource('alphamissense'))}",
        f"--plugin SpliceAI,snv={container_path(annotation_resource('spliceai_snv'))},indel={container_path(annotation_resource('spliceai_indel'))}",
    ))


rule normalize_somatic_for_annotation:
    output:
        vcf=temp(tmp_path("annotation", "somatic", "{tumour}.normalized.vcf.gz")),
        tbi=temp(tmp_path("annotation", "somatic", "{tumour}.normalized.vcf.gz.tbi"))
    input:
        vcf=lambda wildcards: somatic_annotation_source(wildcards.tumour),
        reference=config["reference"]["genome"]
    params:
        tmp=tmp_path("annotation", "somatic", "{tumour}.publish.vcf.gz")
    threads: 4
    shell:
        """
        mkdir -p $(dirname {params.tmp})
        bcftools norm --threads {threads} -f {input.reference} -m -any -Oz -o {params.tmp} {input.vcf}
        tabix -f -p vcf {params.tmp}
        mv {params.tmp} {output.vcf}
        mv {params.tmp}.tbi {output.tbi}
        """


rule annotate_somatic_vep:
    output:
        vcf="results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz",
        tbi="results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz.tbi"
    input:
        vcf=rules.normalize_somatic_for_annotation.output.vcf,
        tbi=rules.normalize_somatic_for_annotation.output.tbi,
        reference=config["reference"]["genome"],
        resources=annotation_input_paths,
        indexes=annotation_index_inputs
    params:
        cmd=lambda wildcards: get_container_cmd(config["annotation"]),
        cache=lambda wildcards: container_path(annotation_resource("cache")),
        plugins=lambda wildcards: container_path(annotation_resource("plugins")),
        reference=lambda wildcards: container_path(config["reference"]["genome"]),
        input=lambda wildcards: container_path(tmp_path("annotation", "somatic", f"{wildcards.tumour}.normalized.vcf.gz")),
        tmp=tmp_path("annotation", "somatic", "{tumour}.annotated.publish.vcf.gz"),
        output=lambda wildcards: container_path(tmp_path("annotation", "somatic", f"{wildcards.tumour}.annotated.publish.vcf.gz")),
        pick_order=",".join(config.get("annotation", {}).get("pick_order", [])),
        plugin_args=vep_plugin_args
    threads: 8
    shell:
        """
        {params.cmd} env PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0 vep \
            --offline --cache --dir_cache {params.cache} --dir_plugins {params.plugins} \
            --assembly GRCh37 --fasta {params.reference} --format vcf --vcf \
            --compress_output bgzip --force_overwrite --fork {threads} \
            --everything --flag_pick_allele --pick_order {params.pick_order} \
            {params.plugin_args} --input_file {params.input} --output_file {params.output}
        tabix -f -p vcf {params.tmp}
        test -s {params.tmp}.tbi
        mkdir -p $(dirname {output.vcf})
        mv {params.tmp} {output.vcf}
        mv {params.tmp}.tbi {output.tbi}
        """


rule normalize_germline_for_annotation:
    output:
        vcf=temp(tmp_path("annotation", "germline", "{normal}.normalized.vcf.gz")),
        tbi=temp(tmp_path("annotation", "germline", "{normal}.normalized.vcf.gz.tbi"))
    input:
        vcf="results/germline/{normal}.haplotypecaller.filtered.vcf.gz",
        tbi="results/germline/{normal}.haplotypecaller.filtered.vcf.gz.tbi",
        reference=config["reference"]["genome"]
    params:
        tmp=tmp_path("annotation", "germline", "{normal}.publish.vcf.gz")
    threads: 4
    shell:
        """
        mkdir -p $(dirname {params.tmp})
        bcftools norm --threads {threads} -f {input.reference} -m -any -Oz -o {params.tmp} {input.vcf}
        tabix -f -p vcf {params.tmp}
        mv {params.tmp} {output.vcf}
        mv {params.tmp}.tbi {output.tbi}
        """


use rule annotate_somatic_vep as annotate_germline_vep with:
    input:
        vcf=rules.normalize_germline_for_annotation.output.vcf,
        tbi=rules.normalize_germline_for_annotation.output.tbi,
        reference=config["reference"]["genome"],
        resources=annotation_input_paths,
        indexes=annotation_index_inputs
    output:
        vcf="results/annotations/germline/{normal}.haplotypecaller.filtered.annotated.vcf.gz",
        tbi="results/annotations/germline/{normal}.haplotypecaller.filtered.annotated.vcf.gz.tbi"
    params:
        input=lambda wildcards: container_path(tmp_path("annotation", "germline", f"{wildcards.normal}.normalized.vcf.gz")),
        tmp=tmp_path("annotation", "germline", "{normal}.annotated.publish.vcf.gz"),
        output=lambda wildcards: container_path(tmp_path("annotation", "germline", f"{wildcards.normal}.annotated.publish.vcf.gz"))


rule aggregate_mutation_burden:
    output:
        "results/aggregate/mutation_burden.tsv"
    input:
        vcfs=expand("results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz", tumour=tumour_samples()),
        manifest=rules.analysis_manifest.output
    params:
        inputs=" ".join(f"--input {sample}=results/annotations/somatic/{sample}.intersect.annotated.vcf.gz" for sample in tumour_samples()),
        tmp=tmp_path("mutation_burden.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/aggregate_annotations.py burden {params.inputs} --manifest {input.manifest} --output {params.tmp}
        mkdir -p results/aggregate
        mv {params.tmp} {output}
        """


rule aggregate_recurrent_variants:
    output:
        "results/aggregate/recurrent_variants.tsv"
    input:
        vcfs=expand("results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz", tumour=tumour_samples())
    params:
        inputs=" ".join(f"--input {sample}=results/annotations/somatic/{sample}.intersect.annotated.vcf.gz" for sample in tumour_samples()),
        minimum=config.get("annotation", {}).get("recurrence_min_carriers", 2),
        tmp=tmp_path("recurrent_variants.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/aggregate_annotations.py recurrence {params.inputs} --minimum-carriers {params.minimum} --output {params.tmp}
        mkdir -p results/aggregate
        mv {params.tmp} {output}
        """


rule annotation_resource_report:
    output:
        "results/aggregate/annotation_resources.tsv"
    input:
        config=config.get("run_management", {}).get("config_file", os.path.join(PIPELINE_DIR, "config/config.yaml")),
        resources=annotation_input_paths,
        indexes=annotation_index_inputs
    params:
        tmp=tmp_path("annotation_resources.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/aggregate_annotations.py resources --config {input.config} --output {params.tmp}
        mkdir -p results/aggregate
        mv {params.tmp} {output}
        """
