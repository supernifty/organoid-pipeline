def extended_contexts_output():
    return "results/aggregate/extended_contexts.tsv"


def extended_context_inputs():
    return [
        f"{tumour}=results/variants/{tumour}.intersect.vcf.gz"
        for tumour in tumour_samples()
    ]


rule aggregate_extended_contexts:
    output:
        tsv=extended_contexts_output()
    input:
        vcfs=variant_outputs(),
        tbis=expand("results/variants/{tumour}.intersect.vcf.gz.tbi", tumour=tumour_samples()),
        reference=config["reference"]["genome"],
        reference_index=f"{config['reference']['genome']}.fai"
    params:
        inputs=" ".join(extended_context_inputs()),
        tmp=tmp_path("extended_contexts.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/extended_contexts.py \
            --reference {input.reference} \
            --reference-index {input.reference_index} \
            --input {params.inputs} \
            --output {params.tmp}

        test -s {params.tmp}
        mv {params.tmp} {output.tsv}
        """
