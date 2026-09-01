def hotspots_output():
    return "results/aggregate/hotspots.tsv"


def hotspot_inputs():
    if VCF_ONLY_MODE:
        return [
            f"{tumour}=results/variants/{tumour}.intersect.vcf.gz"
            for tumour in tumour_samples()
        ]
    return [
        (
            f"{tumour}=results/variants/{tumour}.intersect.vcf.gz"
            f"=results/variants/{tumour}.mutect2.somatic.vcf.gz"
            f"={tmp_path(f'{tumour}.strelka.somatic.snvs.af.norm.vcf.gz')}"
        )
        for tumour in tumour_samples()
    ]


rule aggregate_hotspots:
    output:
        tsv=hotspots_output()
    input:
        hotspots=config.get("hotspots", {}).get("resource", "config/hotspots.tsv"),
        intersect=variant_outputs(),
        intersect_tbi=expand("results/variants/{tumour}.intersect.vcf.gz.tbi", tumour=tumour_samples()),
        mutect2=[] if VCF_ONLY_MODE else expand("results/variants/{tumour}.mutect2.somatic.vcf.gz", tumour=tumour_samples()),
        mutect2_tbi=[] if VCF_ONLY_MODE else expand("results/variants/{tumour}.mutect2.somatic.vcf.gz.tbi", tumour=tumour_samples()),
        strelka=[] if VCF_ONLY_MODE else expand(tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz"), tumour=tumour_samples()),
        strelka_tbi=[] if VCF_ONLY_MODE else expand(tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz.tbi"), tumour=tumour_samples())
    params:
        inputs=" ".join(hotspot_inputs()),
        tmp=tmp_path("hotspots.tsv")
    threads: 1
    shell:
        """
        mkdir -p results/aggregate

        python3 {PIPELINE_DIR}/workflow/scripts/hotspots.py \
            --hotspots {input.hotspots} \
            --input {params.inputs} \
            --output {params.tmp}

        test -s {params.tmp}
        mv {params.tmp} {output.tsv}
        """
