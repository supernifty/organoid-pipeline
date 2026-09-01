MUTATIONAL_SIGNATURE_TOOL = os.path.join(PIPELINE_DIR, "tools/mutational_signature/mutational_signature")
MUTATIONAL_SIGNATURE_CONFIG = config["mutational_signatures"]


def tissue_signature_definition_output(signature_type):
    return f"results/signatures/definitions/tissue.{signature_type}.tsv"


def tissue_signatures(signature_type):
    return MUTATIONAL_SIGNATURE_CONFIG["tissue_signatures"][signature_type]


def mutational_signature_definition(signature_type):
    return MUTATIONAL_SIGNATURE_CONFIG[f"{signature_type}_definition"]


def tissue_signature_exposure_outputs(signature_type):
    return expand(
        f"results/signatures/exposures/{{tumour}}.intersect.tissue.{signature_type}.tsv",
        tumour=tumour_samples(),
    )


def strelka_tolerant_id_exposure_outputs():
    return expand(
        "results/signatures/exposures/{tumour}.strelka.tissue.id.tsv",
        tumour=tumour_samples(),
    )


def aggregate_tissue_signature_output(signature_type):
    return f"results/aggregate/signatures.tissue.{signature_type}.tsv"


def aggregate_strelka_tolerant_id_output():
    return "results/aggregate/signatures.strelka.tissue.id.tsv"


rule tissue_signature_definitions:
    output:
        sbs=tissue_signature_definition_output("sbs"),
        id=tissue_signature_definition_output("id"),
        dbs=tissue_signature_definition_output("dbs")
    input:
        sbs=mutational_signature_definition("sbs"),
        id=mutational_signature_definition("id"),
        dbs=mutational_signature_definition("dbs")
    params:
        sbs=" ".join(tissue_signatures("sbs")),
        id=" ".join(tissue_signatures("id")),
        dbs=" ".join(tissue_signatures("dbs")),
        tmp_sbs=tmp_path("tissue.sbs.tsv"),
        tmp_id=tmp_path("tissue.id.tsv"),
        tmp_dbs=tmp_path("tissue.dbs.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/filter_signature_definitions.py \
            --definition {input.sbs} \
            --output {params.tmp_sbs} \
            --signatures {params.sbs}
        test -s {params.tmp_sbs}
        mv {params.tmp_sbs} {output.sbs}

        python3 {PIPELINE_DIR}/workflow/scripts/filter_signature_definitions.py \
            --definition {input.id} \
            --output {params.tmp_id} \
            --signatures {params.id}
        test -s {params.tmp_id}
        mv {params.tmp_id} {output.id}

        python3 {PIPELINE_DIR}/workflow/scripts/filter_signature_definitions.py \
            --definition {input.dbs} \
            --output {params.tmp_dbs} \
            --signatures {params.dbs}
        test -s {params.tmp_dbs}
        mv {params.tmp_dbs} {output.dbs}
        """


rule mutational_signature_counts:
    output:
        counts="results/signatures/counts/{tumour}.intersect.counts.tsv"
    input:
        vcf="results/variants/{tumour}.intersect.vcf.gz",
        tbi="results/variants/{tumour}.intersect.vcf.gz.tbi",
        reference=config["reference"]["genome"]
    params:
        tmp_counts=tmp_path("{tumour}.intersect.counts.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/count.py \
            --indels \
            --doublets \
            --genome {input.reference} \
            --vcf {input.vcf} \
            > {params.tmp_counts}

        test -s {params.tmp_counts}
        mv {params.tmp_counts} {output.counts}
        """


rule mutational_signature_tissue_sbs:
    output:
        exposures="results/signatures/exposures/{tumour}.intersect.tissue.sbs.tsv"
    input:
        counts="results/signatures/counts/{tumour}.intersect.counts.tsv",
        definition=rules.tissue_signature_definitions.output.sbs
    params:
        tmp_exposures=tmp_path("{tumour}.intersect.tissue.sbs.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/decompose.py \
            --signatures {input.definition} \
            --counts {input.counts} \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule mutational_signature_tissue_id:
    output:
        exposures="results/signatures/exposures/{tumour}.intersect.tissue.id.tsv"
    input:
        counts="results/signatures/counts/{tumour}.intersect.counts.tsv",
        definition=rules.tissue_signature_definitions.output.id
    params:
        tmp_exposures=tmp_path("{tumour}.intersect.tissue.id.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/decompose.py \
            --signatures {input.definition} \
            --counts {input.counts} \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule strelka_tolerant_indel_signature_vcf:
    output:
        vcf=temp(tmp_path("{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz")),
        tbi=temp(tmp_path("{tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf.gz.tbi"))
    input:
        vcf=tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz")
    params:
        sample="TUMOR",
        af=MUTATIONAL_SIGNATURE_CONFIG["strelka_tolerant_id"]["af_threshold"],
        dp=MUTATIONAL_SIGNATURE_CONFIG["strelka_tolerant_id"]["dp_threshold"],
        annotated_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.strelka.somatic.indels.norm.annot.pass.af.vcf"),
        filtered_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.strelka.somatic.indels.norm.annot.pass.af.filter.vcf")
    threads: 1
    shell:
        """
        python3 scripts/annotate_strelka_af.py \
            --mode indel \
            --sample {params.sample} \
            --input {input.vcf} \
            --output {params.annotated_vcf}

        python3 scripts/filter_vcf.py \
            --input {params.annotated_vcf} \
            --output {params.filtered_vcf} \
            --sample {params.sample} \
            --af {params.af} \
            --dp {params.dp} \
            --info-af \
            --pass-only

        bgzip -f {params.filtered_vcf}
        gzip -t {params.filtered_vcf}.gz
        tabix -f -p vcf {params.filtered_vcf}.gz
        test -s {params.filtered_vcf}.gz.tbi
        mv {params.filtered_vcf}.gz {output.vcf}
        mv {params.filtered_vcf}.gz.tbi {output.tbi}
        rm -f {params.annotated_vcf}
        """


rule mutational_signature_strelka_tolerant_id_counts:
    output:
        counts="results/signatures/counts/{tumour}.strelka.id.counts.tsv"
    input:
        vcf=rules.strelka_tolerant_indel_signature_vcf.output.vcf,
        tbi=rules.strelka_tolerant_indel_signature_vcf.output.tbi,
        reference=config["reference"]["genome"]
    params:
        tmp_counts=tmp_path("{tumour}.strelka.id.counts.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/count.py \
            --indels \
            --just_indels \
            --genome {input.reference} \
            --vcf {input.vcf} \
            > {params.tmp_counts}

        test -s {params.tmp_counts}
        mv {params.tmp_counts} {output.counts}
        """


rule mutational_signature_strelka_tolerant_tissue_id:
    output:
        exposures="results/signatures/exposures/{tumour}.strelka.tissue.id.tsv"
    input:
        counts="results/signatures/counts/{tumour}.strelka.id.counts.tsv",
        definition=rules.tissue_signature_definitions.output.id
    params:
        tmp_exposures=tmp_path("{tumour}.strelka.tissue.id.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/decompose.py \
            --signatures {input.definition} \
            --counts {input.counts} \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule mutational_signature_tissue_dbs:
    output:
        exposures="results/signatures/exposures/{tumour}.intersect.tissue.dbs.tsv"
    input:
        counts="results/signatures/counts/{tumour}.intersect.counts.tsv",
        definition=rules.tissue_signature_definitions.output.dbs
    params:
        tmp_exposures=tmp_path("{tumour}.intersect.tissue.dbs.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/decompose.py \
            --signatures {input.definition} \
            --counts {input.counts} \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule aggregate_mutational_signature_tissue_sbs:
    output:
        exposures=aggregate_tissue_signature_output("sbs")
    input:
        exposures=tissue_signature_exposure_outputs("sbs"),
        definition=rules.tissue_signature_definitions.output.sbs
    params:
        tmp_exposures=tmp_path("signatures.tissue.sbs.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/combine_signatures.py \
            --signatures {input.definition} \
            --files {input.exposures} \
            | sed 's/\\.intersect\\.tissue\\.sbs\\.tsv//' \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule aggregate_mutational_signature_tissue_id:
    output:
        exposures=aggregate_tissue_signature_output("id")
    input:
        exposures=tissue_signature_exposure_outputs("id"),
        definition=rules.tissue_signature_definitions.output.id
    params:
        tmp_exposures=tmp_path("signatures.tissue.id.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/combine_signatures.py \
            --signatures {input.definition} \
            --files {input.exposures} \
            | sed 's/\\.intersect\\.tissue\\.id\\.tsv//' \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule aggregate_mutational_signature_strelka_tolerant_tissue_id:
    output:
        exposures=aggregate_strelka_tolerant_id_output()
    input:
        exposures=strelka_tolerant_id_exposure_outputs(),
        definition=rules.tissue_signature_definitions.output.id
    params:
        tmp_exposures=tmp_path("signatures.strelka.tissue.id.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/combine_signatures.py \
            --signatures {input.definition} \
            --files {input.exposures} \
            | sed 's/\\.strelka\\.tissue\\.id\\.tsv//' \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """


rule aggregate_mutational_signature_tissue_dbs:
    output:
        exposures=aggregate_tissue_signature_output("dbs")
    input:
        exposures=tissue_signature_exposure_outputs("dbs"),
        definition=rules.tissue_signature_definitions.output.dbs
    params:
        tmp_exposures=tmp_path("signatures.tissue.dbs.tsv")
    threads: 1
    shell:
        """
        python3 {MUTATIONAL_SIGNATURE_TOOL}/combine_signatures.py \
            --signatures {input.definition} \
            --files {input.exposures} \
            | sed 's/\\.intersect\\.tissue\\.dbs\\.tsv//' \
            > {params.tmp_exposures}

        test -s {params.tmp_exposures}
        mv {params.tmp_exposures} {output.exposures}
        """
