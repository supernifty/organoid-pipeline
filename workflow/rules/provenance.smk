rule provenance:
    output:
        tsv="results/aggregate/provenance.tsv"
    input:
        config=config.get("run_management", {}).get("config_file", os.path.join(PIPELINE_DIR, "config/config.yaml")),
        local_config=[] if config.get("run_management") else (LOCAL_CONFIG_PATH if os.path.exists(LOCAL_CONFIG_PATH) else []),
        samples=samples_path,
        pixi=os.path.join(PIPELINE_DIR, "pixi.toml"),
        pixi_lock=os.path.join(PIPELINE_DIR, "pixi.lock"),
        pyproject=os.path.join(PIPELINE_DIR, "pyproject.toml"),
        uv_lock=os.path.join(PIPELINE_DIR, "uv.lock"),
        slurm=os.path.join(PIPELINE_DIR, "config/slurm/config.yaml")
    params:
        tmp=tmp_path("provenance.tsv"),
        local_config_arg=lambda wildcards, input: f"--config-overlay {input.local_config}" if input.local_config else ""
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/write_provenance.py \
            --config {input.config} \
            {params.local_config_arg} \
            --samples {input.samples} \
            --pixi {input.pixi} \
            --slurm-config {input.slurm} \
            --output {params.tmp}

        test -s {params.tmp}
        mv {params.tmp} {output.tsv}
        """
