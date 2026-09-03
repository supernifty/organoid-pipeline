rule fastqc:
    output:
        touch("results/qc/fastqc/{sample}/fastqc_done")
    input:
        fastqs=lambda wildcards: get_fastqs(wildcards.sample)
    benchmark: "results/benchmarks/fastqc/{sample}.tsv"
    params:
        outdir="results/qc/fastqc/{sample}"
    threads: 4
    resources: mem_mb=4096
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {threads} {input.fastqs}
        """
