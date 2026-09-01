rule somalier_extract:
    output:
        "results/qc/somalier/extracted/{sample}.somalier"
    input:
        cram=lambda wildcards: aligned_cram(wildcards.sample),
        crai=lambda wildcards: aligned_cram_index(wildcards.sample),
        reference=config["reference"]["genome"],
        sites=lambda wildcards: config["somalier"]["sites_vcf"],
        sites_index=lambda wildcards: config["somalier"]["sites_vcf_index"]
    params:
        outdir=tmp_path("somalier", "extract", "{sample}"),
        cmd=get_container_cmd(config["somalier"])
    threads: 2
    shell:
        """
        mkdir -p {params.outdir}
        {params.cmd} env SOMALIER_SAMPLE_NAME={wildcards.sample} somalier extract \
            --sites {input.sites} --fasta {input.reference} --out-dir {params.outdir} {input.cram}
        test -s {params.outdir}/{wildcards.sample}.somalier
        mkdir -p results/qc/somalier/extracted
        mv {params.outdir}/{wildcards.sample}.somalier {output}
        """


rule somalier_expected_groups:
    output:
        "results/qc/somalier/expected_groups.txt"
    input:
        samples=samples_path
    params:
        tmp=tmp_path("somalier", "expected_groups.txt")
    threads: 1
    shell:
        """
        mkdir -p $(dirname {params.tmp})
        python3 {PIPELINE_DIR}/workflow/scripts/somalier_qc.py groups --samples {input.samples} --output {params.tmp}
        mkdir -p results/qc/somalier
        mv {params.tmp} {output}
        """


rule somalier_relate:
    output:
        pairs="results/qc/somalier/cohort.pairs.tsv",
        samples="results/qc/somalier/cohort.samples.tsv",
        groups="results/qc/somalier/cohort.groups.tsv",
        html="results/qc/somalier/cohort.html"
    input:
        extracted=expand("results/qc/somalier/extracted/{sample}.somalier", sample=config["samples"]["samples"]),
        groups=rules.somalier_expected_groups.output
    params:
        prefix=tmp_path("somalier", "cohort"),
        minimum_depth=config.get("somalier", {}).get("minimum_depth", 20),
        cmd=get_container_cmd(config["somalier"])
    threads: 1
    shell:
        """
        {params.cmd} env SOMALIER_REPORT_ALL_PAIRS=1 somalier relate --groups {input.groups} \
            --min-depth {params.minimum_depth} --output-prefix {params.prefix} {input.extracted}
        mkdir -p results/qc/somalier
        mv {params.prefix}.pairs.tsv {output.pairs}
        mv {params.prefix}.samples.tsv {output.samples}
        mv {params.prefix}.groups.tsv {output.groups}
        mv {params.prefix}.html {output.html}
        """


rule somalier_flags:
    output:
        "results/qc/somalier/somalier_flags.tsv"
    input:
        samples=samples_path,
        pairs=rules.somalier_relate.output.pairs
    params:
        expected=config.get("somalier", {}).get("expected_pair_concordance_min", 0.6),
        unexpected=config.get("somalier", {}).get("unexpected_pair_concordance_max", 0.4),
        tmp=tmp_path("somalier", "somalier_flags.tsv")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/somalier_qc.py flags --samples {input.samples} --pairs {input.pairs} \
            --expected-min {params.expected} --unexpected-max {params.unexpected} --output {params.tmp}
        mv {params.tmp} {output}
        """


rule somalier_ancestry:
    output:
        tsv="results/qc/somalier/ancestry.tsv",
        html="results/qc/somalier/ancestry.html"
    input:
        extracted=expand("results/qc/somalier/extracted/{sample}.somalier", sample=config["samples"]["samples"]),
        labels=lambda wildcards: config["somalier"]["ancestry"]["labels"]
    params:
        references=lambda wildcards: os.path.join(container_path(config["somalier"]["ancestry"]["reference_somalier_dir"]), "*.somalier"),
        prefix=tmp_path("somalier", "ancestry"),
        cmd=get_container_cmd(config["somalier"])
    threads: 1
    shell:
        """
        {params.cmd} somalier ancestry --labels {input.labels} --output-prefix {params.prefix} \
            "{params.references}" ++ {input.extracted}
        test -s {params.prefix}.somalier-ancestry.tsv
        test -s {params.prefix}.somalier-ancestry.html
        mkdir -p results/qc/somalier
        mv {params.prefix}.somalier-ancestry.tsv {output.tsv}
        mv {params.prefix}.somalier-ancestry.html {output.html}
        """
