CALLER_TIER_NAMES = ("mutect2.pass", "strelka.pass", "intersection", "union", "mutect2-only", "strelka-only")


def caller_tier_paths(sample):
    return {name.replace("-", "_").replace(".", "_"): f"results/callers/{sample}.{name}.vcf.gz" for name in CALLER_TIER_NAMES}


def cohort_alignment_args(wildcards):
    return " ".join(
        f"--sample {sample}={aligned_cram(sample)}" for sample in sorted(config["samples"]["samples"])
    )


def cohort_union_args(wildcards):
    return " ".join(
        f"--union {sample}=results/callers/{sample}.union.vcf.gz" for sample in tumour_samples()
    )


def optional_mask_inputs(wildcards):
    return [path for path in (
        config["reference"].get("problematic_regions"),
        config["reference"].get("low_mappability_regions"),
        config["reference"].get("repeat_regions"),
    ) if path]


def optional_mask_args(wildcards, input):
    labels = ("problematic", "low_mappability", "repeat")
    return " ".join(
        f"--mask {label}={config['reference'][key]}"
        for label, key in zip(labels, ("problematic_regions", "low_mappability_regions", "repeat_regions"))
        if config["reference"].get(key)
    )


rule caller_tiers:
    input:
        mutect2="results/variants/{tumour}.mutect2.somatic.vcf.gz",
        mutect2_tbi="results/variants/{tumour}.mutect2.somatic.vcf.gz.tbi",
        strelka_snvs=tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz"),
        strelka_snvs_tbi=tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz.tbi"),
        strelka_indels=tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz"),
        strelka_indels_tbi=tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz.tbi")
    output:
        mutect2_pass="results/callers/{tumour}.mutect2.pass.vcf.gz",
        mutect2_pass_tbi="results/callers/{tumour}.mutect2.pass.vcf.gz.tbi",
        strelka_pass="results/callers/{tumour}.strelka.pass.vcf.gz",
        strelka_pass_tbi="results/callers/{tumour}.strelka.pass.vcf.gz.tbi",
        intersection="results/callers/{tumour}.intersection.vcf.gz",
        intersection_tbi="results/callers/{tumour}.intersection.vcf.gz.tbi",
        union="results/callers/{tumour}.union.vcf.gz",
        union_tbi="results/callers/{tumour}.union.vcf.gz.tbi",
        mutect2_only="results/callers/{tumour}.mutect2-only.vcf.gz",
        mutect2_only_tbi="results/callers/{tumour}.mutect2-only.vcf.gz.tbi",
        strelka_only="results/callers/{tumour}.strelka-only.vcf.gz",
        strelka_only_tbi="results/callers/{tumour}.strelka-only.vcf.gz.tbi"
    log: "logs/caller_tiers/{tumour}.log"
    benchmark: "results/benchmarks/caller_tiers/{tumour}.tsv"
    params:
        prefix=tmp_path("caller_tiers", "{tumour}"),
        script=os.path.join(PIPELINE_DIR, "workflow/scripts/caller_tiers.py")
    threads: 1
    resources: mem_mb=2048, runtime=60, disk_mb=4096
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix}) results/callers $(dirname {log})
        python3 {params.script} --mutect2 {input.mutect2} --strelka {input.strelka_snvs} {input.strelka_indels} \
          --mutect2-pass {params.prefix}.mutect2.pass.vcf --strelka-pass {params.prefix}.strelka.pass.vcf \
          --intersection {params.prefix}.intersection.vcf --union {params.prefix}.union.vcf \
          --mutect2-only {params.prefix}.mutect2-only.vcf --strelka-only {params.prefix}.strelka-only.vcf > {log} 2>&1
        for tier in mutect2.pass strelka.pass intersection union mutect2-only strelka-only; do
          bgzip -f "{params.prefix}.${{tier}}.vcf"
          tabix -f -p vcf "{params.prefix}.${{tier}}.vcf.gz"
        done
        mv {params.prefix}.mutect2.pass.vcf.gz {output.mutect2_pass}; mv {params.prefix}.mutect2.pass.vcf.gz.tbi {output.mutect2_pass_tbi}
        mv {params.prefix}.strelka.pass.vcf.gz {output.strelka_pass}; mv {params.prefix}.strelka.pass.vcf.gz.tbi {output.strelka_pass_tbi}
        mv {params.prefix}.intersection.vcf.gz {output.intersection}; mv {params.prefix}.intersection.vcf.gz.tbi {output.intersection_tbi}
        mv {params.prefix}.union.vcf.gz {output.union}; mv {params.prefix}.union.vcf.gz.tbi {output.union_tbi}
        mv {params.prefix}.mutect2-only.vcf.gz {output.mutect2_only}; mv {params.prefix}.mutect2-only.vcf.gz.tbi {output.mutect2_only_tbi}
        mv {params.prefix}.strelka-only.vcf.gz {output.strelka_only}; mv {params.prefix}.strelka-only.vcf.gz.tbi {output.strelka_only_tbi}
        """


rule cohort_candidate_union:
    input:
        vcfs=expand("results/callers/{tumour}.union.vcf.gz", tumour=tumour_samples()),
        tbis=expand("results/callers/{tumour}.union.vcf.gz.tbi", tumour=tumour_samples())
    output:
        vcf="results/cohort/candidates.union.vcf.gz",
        tbi="results/cohort/candidates.union.vcf.gz.tbi"
    log: "logs/cohort/candidate_union.log"
    benchmark: "results/benchmarks/cohort_candidate_union.tsv"
    params:
        inputs=lambda wildcards: " ".join(f"--input {sample}=results/callers/{sample}.union.vcf.gz" for sample in tumour_samples()),
        temporary=tmp_path("cohort.candidates.union.vcf")
    threads: 1
    resources: mem_mb=2048, runtime=60, disk_mb=4096
    shell:
        """
        set -euo pipefail
        mkdir -p results/cohort $(dirname {log})
        python3 {PIPELINE_DIR}/workflow/scripts/cohort_union.py {params.inputs} --output {params.temporary} > {log} 2>&1
        bgzip -f {params.temporary}; tabix -f -p vcf {params.temporary}.gz
        mv {params.temporary}.gz {output.vcf}; mv {params.temporary}.gz.tbi {output.tbi}
        """


rule cohort_allele_recount:
    input:
        candidates=rules.cohort_candidate_union.output.vcf,
        candidates_tbi=rules.cohort_candidate_union.output.tbi,
        reference=config["reference"]["genome"],
        crams=lambda wildcards: [aligned_cram(sample) for sample in sorted(config["samples"]["samples"])],
        crais=lambda wildcards: [aligned_cram_index(sample) for sample in sorted(config["samples"]["samples"])]
    output: "results/cohort/allele_counts.tsv"
    log: "logs/cohort/allele_recount.log"
    benchmark: "results/benchmarks/cohort_allele_recount.tsv"
    params: sample_args=cohort_alignment_args, temporary=tmp_path("cohort.allele_counts.tsv")
    threads: 1
    resources: mem_mb=4096, runtime=1440, disk_mb=4096
    shell:
        """
        set -euo pipefail
        mkdir -p results/cohort $(dirname {log})
        python3 {PIPELINE_DIR}/workflow/scripts/allele_recount.py --candidates {input.candidates} \
          --reference {input.reference} {params.sample_args} --output {params.temporary} > {log} 2>&1
        test -s {params.temporary}; mv {params.temporary} {output}
        """


rule filter_organoid_catalog:
    input:
        union="results/callers/{tumour}.union.vcf.gz",
        cohort_unions=expand("results/callers/{tumour}.union.vcf.gz", tumour=tumour_samples()),
        counts=rules.cohort_allele_recount.output,
        population=config["reference"]["population_vcf"],
        population_tbi=lambda wildcards: f"{config['reference']['population_vcf']}.tbi",
        manifest=samples_path,
        masks=optional_mask_inputs
    output:
        audit_vcf="results/catalogs/{tumour}.audit.vcf.gz", audit_tbi="results/catalogs/{tumour}.audit.vcf.gz.tbi", audit_tsv="results/catalogs/{tumour}.audit.tsv",
        stringent_vcf="results/catalogs/{tumour}.stringent.vcf.gz", stringent_tbi="results/catalogs/{tumour}.stringent.vcf.gz.tbi", stringent_tsv="results/catalogs/{tumour}.stringent.tsv",
        sensitivity_vcf="results/catalogs/{tumour}.sensitivity.vcf.gz", sensitivity_tbi="results/catalogs/{tumour}.sensitivity.vcf.gz.tbi", sensitivity_tsv="results/catalogs/{tumour}.sensitivity.tsv",
        rejected_vcf="results/catalogs/{tumour}.rejected.vcf.gz", rejected_tbi="results/catalogs/{tumour}.rejected.vcf.gz.tbi", rejected_tsv="results/catalogs/{tumour}.rejected.tsv",
        shared_vcf="results/catalogs/{tumour}.shared-lineage.vcf.gz", shared_tbi="results/catalogs/{tumour}.shared-lineage.vcf.gz.tbi", shared_tsv="results/catalogs/{tumour}.shared-lineage.tsv",
        stage_counts="results/catalogs/{tumour}.stage_counts.tsv"
    log: "logs/catalogs/{tumour}.log"
    benchmark: "results/benchmarks/catalogs/{tumour}.tsv"
    params:
        union_args=cohort_union_args, mask_args=optional_mask_args,
        prefix=tmp_path("catalogs", "{tumour}"),
        af_field=lambda wildcards: config["reference"].get("population_af_field", "AF"),
        pop_af=config["filtering"]["population_af_threshold"],
        min_alt=config["filtering"]["min_later_alt_reads"], min_vaf=config["filtering"]["later_vaf_threshold"],
        min_base_dp=config["filtering"]["minimum_baseline_depth"], max_base_alt=config["filtering"]["maximum_baseline_alt_reads"]
    threads: 1
    resources: mem_mb=4096, runtime=240, disk_mb=4096
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/catalogs $(dirname {params.prefix}) $(dirname {log})
        python3 {PIPELINE_DIR}/workflow/scripts/catalog_filter.py --sample {wildcards.tumour} --manifest {input.manifest} \
          {params.union_args} --counts {input.counts} --population-vcf {input.population} --population-af-field {params.af_field} \
          --population-af-threshold {params.pop_af} --min-later-alt {params.min_alt} --min-later-vaf {params.min_vaf} \
          --min-baseline-depth {params.min_base_dp} --max-baseline-alt {params.max_base_alt} {params.mask_args} \
          --audit-vcf {params.prefix}.audit.vcf --audit-tsv {params.prefix}.audit.tsv \
          --stringent-vcf {params.prefix}.stringent.vcf --stringent-tsv {params.prefix}.stringent.tsv \
          --sensitivity-vcf {params.prefix}.sensitivity.vcf --sensitivity-tsv {params.prefix}.sensitivity.tsv \
          --rejected-vcf {params.prefix}.rejected.vcf --rejected-tsv {params.prefix}.rejected.tsv \
          --shared-lineage-vcf {params.prefix}.shared-lineage.vcf --shared-lineage-tsv {params.prefix}.shared-lineage.tsv \
          --stage-counts {params.prefix}.stage_counts.tsv > {log} 2>&1
        for tier in audit stringent sensitivity rejected shared-lineage; do bgzip -f "{params.prefix}.${{tier}}.vcf"; tabix -f -p vcf "{params.prefix}.${{tier}}.vcf.gz"; done
        mv {params.prefix}.audit.vcf.gz {output.audit_vcf}; mv {params.prefix}.audit.vcf.gz.tbi {output.audit_tbi}; mv {params.prefix}.audit.tsv {output.audit_tsv}
        mv {params.prefix}.stringent.vcf.gz {output.stringent_vcf}; mv {params.prefix}.stringent.vcf.gz.tbi {output.stringent_tbi}; mv {params.prefix}.stringent.tsv {output.stringent_tsv}
        mv {params.prefix}.sensitivity.vcf.gz {output.sensitivity_vcf}; mv {params.prefix}.sensitivity.vcf.gz.tbi {output.sensitivity_tbi}; mv {params.prefix}.sensitivity.tsv {output.sensitivity_tsv}
        mv {params.prefix}.rejected.vcf.gz {output.rejected_vcf}; mv {params.prefix}.rejected.vcf.gz.tbi {output.rejected_tbi}; mv {params.prefix}.rejected.tsv {output.rejected_tsv}
        mv {params.prefix}.shared-lineage.vcf.gz {output.shared_vcf}; mv {params.prefix}.shared-lineage.vcf.gz.tbi {output.shared_tbi}; mv {params.prefix}.shared-lineage.tsv {output.shared_tsv}
        mv {params.prefix}.stage_counts.tsv {output.stage_counts}
        """


rule sbs96_catalogs:
    input:
        reference=config["reference"]["genome"],
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai",
        mutect2="results/callers/{tumour}.mutect2.pass.vcf.gz",
        strelka="results/callers/{tumour}.strelka.pass.vcf.gz",
        intersection="results/callers/{tumour}.intersection.vcf.gz",
        stringent="results/catalogs/{tumour}.stringent.vcf.gz"
    output: "results/signatures/{tumour}.sbs96.tsv"
    log: "logs/signatures/{tumour}.sbs96.log"
    benchmark: "results/benchmarks/signatures/{tumour}.sbs96.tsv"
    params: temporary=tmp_path("signatures", "{tumour}.sbs96.tsv")
    threads: 1
    resources: mem_mb=2048, runtime=120, disk_mb=1024
    shell:
        """
        set -euo pipefail
        mkdir -p results/signatures $(dirname {params.temporary}) $(dirname {log})
        python3 {PIPELINE_DIR}/workflow/scripts/sbs96.py --reference {input.reference} \
          --catalog mutect2={input.mutect2} --catalog strelka2={input.strelka} \
          --catalog intersection={input.intersection} --catalog stringent={input.stringent} \
          --output {params.temporary} > {log} 2>&1
        test -s {params.temporary}; mv {params.temporary} {output}
        """
