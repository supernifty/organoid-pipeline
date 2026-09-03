def panel_of_normals_output():
    return config.get("mutect2", {}).get(
        "generated_panel_of_normals",
        "results/variants/mutect2.pon.vcf.gz",
    )

def panel_of_normals_rule_output():
    return panel_of_normals_output()

def configured_panel_of_normals():
    return config["reference"].get("panel_of_normals")

def should_build_panel_of_normals():
    mutect2_config = config.get("mutect2", {})
    return (
        mutect2_config.get("create_panel_of_normals", True)
        and configured_panel_of_normals() is None
        and len(normal_samples()) > 0
    )

def panel_of_normals_dependency(wildcards):
    if configured_panel_of_normals():
        return configured_panel_of_normals()
    if should_build_panel_of_normals():
        return panel_of_normals_output()
    return []

def panel_of_normals_index_dependency(wildcards):
    pon = panel_of_normals_dependency(wildcards)
    if isinstance(pon, str):
        return f"{pon}.tbi"
    return []

def mutect2_interval_padding():
    return config.get("mutect2", {}).get("interval_padding", 1000)

def final_region_padding():
    if ANALYSIS_TYPE == "wgs":
        return 0
    return config.get("filtering", {}).get("final_region_padding", 0)

def mutect2_pon_max_mnp_distance():
    return config.get("mutect2", {}).get("pon_max_mnp_distance", 0)

def mutect2_reference_args():
    if REFERENCE_BUILD == "grch38":
        return "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter"
    return ""

def pon_shard_slurm_extra():
    configured = os.environ.get("SLURM_EXTRA", "").strip()
    if re.search(r"(^|\s)--tmp(?:=|\s)", configured):
        return configured
    scratch_mb = int(config.get("storage", {}).get("pon_genomicsdb_scratch_mb", 20480))
    if scratch_mb <= 0:
        raise ValueError("storage.pon_genomicsdb_scratch_mb must be positive")
    return " ".join(value for value in (configured, f"--tmp={scratch_mb}M") if value)

def shard_core_interval(wildcards):
    if ANALYSIS_TYPE == "wgs":
        return tmp_path("analysis", "mutect2_shards", f"{wildcards.chromosome}-scattered.interval_list")
    return config["reference"]["regions"]

def shard_dependency(wildcards):
    if ANALYSIS_TYPE == "wgs":
        return rules.split_wgs_intervals.output.shards
    return config["reference"]["regions"]

def shard_contig_arg(wildcards):
    if ANALYSIS_TYPE == "wes":
        return f"-L {wildcards.chromosome} --interval-set-rule INTERSECTION"
    return ""

def wgs_exclusion_input(wildcards):
    return config["reference"].get("wgs_exclude_regions") or []

def reference_profile_vcfs():
    paths = [config["reference"]["gnomad"], config["reference"]["population_vcf"]]
    contamination = config["reference"].get("contamination_sites")
    pon = configured_panel_of_normals()
    for path in (contamination, pon):
        if path and path not in paths:
            paths.append(path)
    return paths

def reference_profile_vcf_indexes():
    return [f"{path}.tbi" for path in reference_profile_vcfs()]

rule prepare_wgs_territory:
    output:
        bed=tmp_path("analysis", "wgs.callable.bed"),
        intervals=tmp_path("analysis", "wgs.callable.interval_list"),
        metadata=tmp_path("analysis", "wgs.callable.json")
    input:
        source=lambda wildcards: config["reference"]["wgs_calling_regions"],
        fai=lambda wildcards: f"{config['reference']['genome']}.fai",
        exclude=wgs_exclusion_input
    params:
        contigs=" ".join(f"--contig {contig}" for contig in ANALYSIS["contigs"]),
        exclude_arg=lambda wildcards, input: f"--exclude {input.exclude}" if input.exclude else ""
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/analysis_mode.py prepare-territory \
            --source {input.source} \
            --fai {input.fai} \
            {params.contigs} \
            {params.exclude_arg} \
            --bed-output {output.bed} \
            --interval-output {output.intervals} \
            --metadata-output {output.metadata}
        """

rule index_wgs_territory:
    output:
        bed=tmp_path("analysis", "wgs.callable.bed.gz"),
        tbi=tmp_path("analysis", "wgs.callable.bed.gz.tbi")
    input:
        bed=rules.prepare_wgs_territory.output.bed
    threads: 1
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -f -p bed {output.bed}
        """

rule split_wgs_intervals:
    output:
        shards=directory(tmp_path("analysis", "mutect2_shards"))
    input:
        intervals=rules.prepare_wgs_territory.output.intervals,
        reference=config["reference"]["genome"]
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        ref_container=container_path(config["reference"]["genome"]),
        intervals_container=container_path(rules.prepare_wgs_territory.output.intervals),
        scatter_count=len(SCATTER_IDS),
        output_container=container_path(tmp_path("analysis", "mutect2_shards"))
    threads: 1
    resources: mem_mb=8192
    shell:
        """
        rm -rf {output.shards}
        mkdir -p {output.shards}
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' SplitIntervals \
            -R {params.ref_container} \
            -L {params.intervals_container} \
            --scatter-count {params.scatter_count} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --subdivision-mode INTERVAL_SUBDIVISION \
            -O {params.output_container}
        test "$(find {output.shards} -name '*-scattered.interval_list' | wc -l | tr -d ' ')" -eq {params.scatter_count}
        """

rule validate_reference_profile:
    output:
        tmp_path("analysis", "reference_profile.json")
    input:
        reference=config["reference"]["genome"],
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai",
        reference_dict=config["reference"]["genome_dict"],
        bwa_indexes=bwa_index_inputs(),
        territory=lambda wildcards: analysis_territory(),
        territory_tbi=lambda wildcards: analysis_territory_index(),
        vcfs=reference_profile_vcfs(),
        vcf_indexes=reference_profile_vcf_indexes(),
        capture_metadata=lambda wildcards: REFERENCE["capture_metadata"] or []
    params:
        build=REFERENCE_BUILD,
        chromosomes=" ".join(f"--chromosome {chrom}" for chrom in config["chromosomes"]),
        vcfs=" ".join(f"--vcf {path}" for path in reference_profile_vcfs()),
        capture_source=REFERENCE["capture_source"],
        capture_metadata_arg=lambda wildcards, input: (
            f"--capture-metadata {input.capture_metadata}" if input.capture_metadata else ""
        )
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/analysis_mode.py validate-reference-profile \
            --build {params.build} \
            --reference-fai {input.reference_fai} \
            --reference-dict {input.reference_dict} \
            --territory {input.territory} \
            --capture-source {params.capture_source} \
            {params.capture_metadata_arg} \
            {params.chromosomes} \
            {params.vcfs} \
            --output {output}
        """

rule analysis_manifest:
    output:
        "results/analysis_manifest.json"
    input:
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai",
        reference_dict=config["reference"]["genome_dict"],
        territory=lambda wildcards: analysis_territory(),
        territory_tbi=lambda wildcards: analysis_territory_index(),
        profile=rules.validate_reference_profile.output,
        capture_metadata=lambda wildcards: REFERENCE["capture_metadata"] or []
    benchmark: "results/benchmarks/aggregate/analysis_manifest.tsv"
    params:
        mode=ANALYSIS_TYPE,
        build=REFERENCE_BUILD,
        capture_source=REFERENCE["capture_source"],
        capture_metadata_arg=lambda wildcards, input: (
            f"--capture-metadata {input.capture_metadata}" if input.capture_metadata else ""
        ),
        tmp=tmp_path("analysis_manifest.json")
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/analysis_mode.py write-manifest \
            --mode {params.mode} \
            --build {params.build} \
            --capture-source {params.capture_source} \
            {params.capture_metadata_arg} \
            --reference-fai {input.reference_fai} \
            --reference-dict {input.reference_dict} \
            --territory {input.territory} \
            --output {params.tmp}
        mkdir -p results
        mv {params.tmp} {output}
        """

rule strelka_call_regions:
    output:
        bed=temp(tmp_path("strelka.callRegions.bed.gz")),
        tbi=temp(tmp_path("strelka.callRegions.bed.gz.tbi"))
    input:
        regions=lambda wildcards: analysis_territory(),
        regions_tbi=lambda wildcards: analysis_territory_index(),
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai"
    params:
        padding=mutect2_interval_padding() if ANALYSIS_TYPE == "wes" else 0,
        genome_file=tmp_path("strelka.callRegions.genome")
    threads: 1
    shell:
        """
        cut -f1,2 {input.reference_fai} > {params.genome_file}
        gzip -dc {input.regions} | \
            bedtools slop -i - -g {params.genome_file} -b {params.padding} | \
            bgzip -c > {output.bed}
        tabix -f -p bed {output.bed}
        rm -f {params.genome_file}
        """

rule final_filter_regions:
    output:
        bed=temp(tmp_path("final_filter_regions.bed.gz")),
        tbi=temp(tmp_path("final_filter_regions.bed.gz.tbi"))
    input:
        regions=lambda wildcards: analysis_territory(),
        regions_tbi=lambda wildcards: analysis_territory_index(),
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai"
    params:
        padding=final_region_padding(),
        genome_file=tmp_path("final_filter_regions.genome")
    threads: 1
    shell:
        """
        cut -f1,2 {input.reference_fai} > {params.genome_file}
        gzip -dc {input.regions} | \
            bedtools slop -i - -g {params.genome_file} -b {params.padding} | \
            bgzip -c > {output.bed}
        tabix -f -p bed {output.bed}
        rm -f {params.genome_file}
        """

rule mutect2_sample_pon_chromosome:
    output:
        vcf=temp(tmp_path("{normal}.{chromosome}.mutect2.pon.vcf.gz")),
        tbi=temp(tmp_path("{normal}.{chromosome}.mutect2.pon.vcf.gz.tbi"))
    input:
        cram=lambda wildcards: aligned_cram(wildcards.normal),
        crai=lambda wildcards: aligned_cram_index(wildcards.normal),
        reference=config["reference"]["genome"],
        shard_dependency=shard_dependency,
        manifest=rules.analysis_manifest.output
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        ref_container=container_path(config["reference"]["genome"]),
        core_container=lambda wildcards: container_path(shard_core_interval(wildcards)),
        shard_contig_arg=shard_contig_arg,
        cram_container=lambda wildcards: container_path(aligned_cram(wildcards.normal)),
        normal_bam_sample=lambda wildcards: bam_sample_name(wildcards.normal),
        raw_vcf=tmp_path("{normal}.{chromosome}.mutect2.pon.raw.vcf.gz"),
        tmp_vcf=tmp_path("{normal}.{chromosome}.mutect2.pon.publish.vcf.gz"),
        max_mnp_distance=mutect2_pon_max_mnp_distance(),
        interval_padding_arg=lambda wildcards: (
            f"--interval-padding {mutect2_interval_padding()}"
            if mutect2_interval_padding()
            else ""
        ),
        reference_args=mutect2_reference_args()
    threads: 4
    resources:
        mutect2_shards=1,
        mem_mb=16384,
        runtime=720
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' Mutect2 \
            -R {params.ref_container} \
            -I {params.cram_container} \
            --tumor-sample {params.normal_bam_sample} \
            {params.reference_args} \
            -O /data/{params.raw_vcf} \
            -L {params.core_container} \
            {params.shard_contig_arg} \
            --max-mnp-distance {params.max_mnp_distance} \
            {params.interval_padding_arg}

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' SelectVariants \
            -R {params.ref_container} \
            -V /data/{params.raw_vcf} \
            -L {params.core_container} \
            {params.shard_contig_arg} \
            --interval-padding 0 \
            -O /data/{params.tmp_vcf}

        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        rm -f {params.raw_vcf} {params.raw_vcf}.tbi {params.raw_vcf}.stats
        """

rule mutect2_pon_fingerprint:
    output:
        tmp_path("mutect2.pon.fingerprint.json")
    input:
        reference_fai=lambda wildcards: f"{config['reference']['genome']}.fai",
        reference_dict=config["reference"]["genome_dict"],
        territory=lambda wildcards: analysis_territory(),
        territory_tbi=lambda wildcards: analysis_territory_index(),
        samples=samples_path,
        manifest=rules.analysis_manifest.output
    params:
        mode=ANALYSIS_TYPE,
        build=REFERENCE_BUILD,
        normals=" ".join(f"--normal {normal}" for normal in normal_samples())
    threads: 1
    shell:
        """
        python3 {PIPELINE_DIR}/workflow/scripts/analysis_mode.py write-pon-fingerprint \
            --mode {params.mode} \
            --build {params.build} \
            --reference-fai {input.reference_fai} \
            --reference-dict {input.reference_dict} \
            --territory {input.territory} \
            --samples {input.samples} \
            {params.normals} \
            --output {output}
        """

rule mutect2_pon_shard:
    output:
        vcf=temp(tmp_path("pon", "{chromosome}.mutect2.pon.vcf.gz")),
        tbi=temp(tmp_path("pon", "{chromosome}.mutect2.pon.vcf.gz.tbi"))
    input:
        vcfs=lambda wildcards: [tmp_path(f"{normal}.{wildcards.chromosome}.mutect2.pon.vcf.gz") for normal in normal_samples()],
        tbis=lambda wildcards: [tmp_path(f"{normal}.{wildcards.chromosome}.mutect2.pon.vcf.gz.tbi") for normal in normal_samples()],
        shard_dependency=shard_dependency,
        reference=config["reference"]["genome"],
        fingerprint=rules.mutect2_pon_fingerprint.output
    params:
        inputs=lambda wildcards: " ".join(
            f"-V {container_path(tmp_path(f'{normal}.{wildcards.chromosome}.mutect2.pon.vcf.gz'))}"
            for normal in normal_samples()
        ),
        core_container=lambda wildcards: container_path(shard_core_interval(wildcards)),
        shard_contig_arg=shard_contig_arg,
        gatk_cmd=get_container_cmd(config["gatk"], bind_node_tmp=True),
        java_options=gatk_java_options("12g"),
        ref_container=container_path(config["reference"]["genome"]),
        workspace=lambda wildcards: scratch_path(f"pon_genomicsdb_{wildcards.chromosome}"),
        workspace_container=lambda wildcards: container_path(scratch_path(f"pon_genomicsdb_{wildcards.chromosome}")),
        node_scratch_template=lambda wildcards: f"/tmp/somatic-pipeline-${{SLURM_JOB_ID}}-pon-{wildcards.chromosome}.XXXXXX",
        raw_vcf=tmp_path("pon", "{chromosome}.mutect2.pon.raw.vcf.gz"),
        tmp_vcf=tmp_path("pon", "{chromosome}.mutect2.pon.publish.vcf.gz")
    threads: 4
    resources:
        mutect2_shards=1,
        mem_mb=16384,
        runtime=720,
        slurm_extra=pon_shard_slurm_extra()
    shell:
        """
        set -euo pipefail
        if [ -n "${{SLURM_JOB_ID:-}}" ]; then
            scratch_root=$(mktemp -d "{params.node_scratch_template}")
            workspace="$scratch_root/genomicsdb"
            workspace_container="/local_scratch${{scratch_root#/tmp}}/genomicsdb"
            cleanup_target="$scratch_root"
        else
            workspace="{params.workspace}"
            workspace_container="{params.workspace_container}"
            cleanup_target="$workspace"
        fi
        cleanup_pon_workspace() {{
            status=$?
            rm -rf -- "$cleanup_target"
            exit "$status"
        }}
        trap cleanup_pon_workspace EXIT
        rm -rf -- "$workspace"
        mkdir -p "$(dirname "$workspace")" "$(dirname {params.tmp_vcf})"
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' GenomicsDBImport \
            {params.inputs} \
            --genomicsdb-workspace-path "$workspace_container" \
            --merge-input-intervals true \
            -L {params.core_container} \
            {params.shard_contig_arg}
        test -s "$workspace/callset.json" || {{ echo "ERROR: GenomicsDBImport did not create callset.json in $workspace" >&2; exit 1; }}
        test -s "$workspace/vidmap.json" || {{ echo "ERROR: GenomicsDBImport did not create vidmap.json in $workspace" >&2; exit 1; }}
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' CreateSomaticPanelOfNormals \
            -R {params.ref_container} \
            -V "gendb://$workspace_container" \
            -O /data/{params.raw_vcf}
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' SelectVariants \
            -R {params.ref_container} \
            -V /data/{params.raw_vcf} \
            -L {params.core_container} \
            {params.shard_contig_arg} \
            --interval-padding 0 \
            -O /data/{params.tmp_vcf}
        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        rm -f {params.raw_vcf} {params.raw_vcf}.tbi
        """

rule mutect2_panel_of_normals:
    output:
        vcf=panel_of_normals_rule_output(),
        tbi=f"{panel_of_normals_rule_output()}.tbi",
        manifest=f"{panel_of_normals_rule_output()}.manifest.json"
    input:
        vcfs=lambda wildcards: [tmp_path("pon", f"{shard}.mutect2.pon.vcf.gz") for shard in SCATTER_IDS],
        tbis=lambda wildcards: [tmp_path("pon", f"{shard}.mutect2.pon.vcf.gz.tbi") for shard in SCATTER_IDS],
        reference=config["reference"]["genome"],
        fingerprint=rules.mutect2_pon_fingerprint.output
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options(),
        ref_container=container_path(config["reference"]["genome"]),
        inputs=" ".join(
            f"-I {container_path(tmp_path('pon', f'{shard}.mutect2.pon.vcf.gz'))}"
            for shard in SCATTER_IDS
        ),
        tmp_vcf=tmp_path("mutect2.pon.publish.vcf.gz")
    threads: 2
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' MergeVcfs \
            {params.inputs} \
            -O /data/{params.tmp_vcf}

        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        mkdir -p results/variants
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        cp {input.fingerprint} {output.manifest}
        """

rule strelka_somatic:
    output:
        snvs="results/variants/{tumour}.strelka.somatic.snvs.vcf.gz",
        snvs_tbi="results/variants/{tumour}.strelka.somatic.snvs.vcf.gz.tbi",
        indels="results/variants/{tumour}.strelka.somatic.indels.vcf.gz",
        indels_tbi="results/variants/{tumour}.strelka.somatic.indels.vcf.gz.tbi"
    input:
        unpack(tumour_normal_crams),
        reference=config["reference"]["genome"],
        regions=rules.strelka_call_regions.output.bed,
        regions_tbi=rules.strelka_call_regions.output.tbi,
        manifest=rules.analysis_manifest.output
    log: "logs/strelka/{tumour}.log"
    benchmark: "results/benchmarks/strelka/{tumour}.tsv"
    params:
        ref=config["reference"]["genome"],
        strelka_params=strelka_mode_params(),
        strelka_image=get_container_image(config["strelka"]),
        run_dir=scratch_path("strelka_{tumour}"),
        run_dir_container=lambda wildcards: container_path(scratch_path(f"strelka_{wildcards.tumour}")),
        strelka_cmd=get_container_cmd(config["strelka"]),
        ref_container=container_path(config["reference"]["genome"]),
        regions_container=container_path(rules.strelka_call_regions.output.bed),
        tumour_container=lambda wildcards: container_path(aligned_cram(wildcards.tumour)),
        normal_container=lambda wildcards: container_path(aligned_cram(config["samples"]["tumours"][wildcards.tumour]))
    threads: 16 if ANALYSIS_TYPE == "wgs" else 8
    resources:
        disk_mb=lambda wildcards, input: max(16384, int(input.size_mb * 0.75)),
        mem_mb=32768,
        runtime=1440 if ANALYSIS_TYPE == "wgs" else 360
    shell:
        """
        mkdir -p {params.run_dir} $(dirname {log})

        {params.strelka_cmd} \
            configureStrelkaSomaticWorkflow.py \
            --referenceFasta {params.ref_container} \
            --tumorBam {params.tumour_container} \
            --normalBam {params.normal_container} \
            --callRegions {params.regions_container} \
            --runDir {params.run_dir_container} \
            {params.strelka_params} > {log} 2>&1

        {params.strelka_cmd} \
            python {params.run_dir_container}/runWorkflow.py -m local -j {threads} >> {log} 2>&1

        gzip -t {params.run_dir}/results/variants/somatic.snvs.vcf.gz
        gzip -t {params.run_dir}/results/variants/somatic.indels.vcf.gz
        tabix -f -p vcf {params.run_dir}/results/variants/somatic.snvs.vcf.gz
        tabix -f -p vcf {params.run_dir}/results/variants/somatic.indels.vcf.gz
        test -s {params.run_dir}/results/variants/somatic.snvs.vcf.gz.tbi
        test -s {params.run_dir}/results/variants/somatic.indels.vcf.gz.tbi
        mkdir -p results/variants
        mv {params.run_dir}/results/variants/somatic.snvs.vcf.gz {output.snvs}
        mv {params.run_dir}/results/variants/somatic.indels.vcf.gz {output.indels}
        mv {params.run_dir}/results/variants/somatic.snvs.vcf.gz.tbi {output.snvs_tbi}
        mv {params.run_dir}/results/variants/somatic.indels.vcf.gz.tbi {output.indels_tbi}

        rm -rf {params.run_dir}
        """

rule strelka_annotate_af_snvs:
    output:
        vcf=temp(tmp_path("{tumour}.strelka.somatic.snvs.af.vcf.gz")),
        tbi=temp(tmp_path("{tumour}.strelka.somatic.snvs.af.vcf.gz.tbi"))
    input:
        vcf="results/variants/{tumour}.strelka.somatic.snvs.vcf.gz"
    benchmark: "results/benchmarks/strelka/{tumour}.annotate_af.tsv"
    params:
        af_vcf=tmp_path("{tumour}.strelka.somatic.snvs.af.vcf")
    threads: 2
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/annotate_strelka_af.py \
            --mode snv \
            --sample TUMOR \
            --input {input.vcf} \
            --output {params.af_vcf}

        bgzip -f {params.af_vcf}
        gzip -t {output.vcf}
        tabix -f -p vcf {output.vcf}
        test -s {output.tbi}
        """

rule strelka_normalize:
    output:
        snvs=temp(tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz")),
        snvs_tbi=temp(tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz.tbi")),
        indels=temp(tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz")),
        indels_tbi=temp(tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz.tbi"))
    input:
        reference=config["reference"]["genome"],
        snvs=tmp_path("{tumour}.strelka.somatic.snvs.af.vcf.gz"),
        indels="results/variants/{tumour}.strelka.somatic.indels.vcf.gz"
    benchmark: "results/benchmarks/strelka/{tumour}.normalize.tsv"
    params:
        vt_decompose=config.get("vt", {}).get("decompose_params", "-s")
    threads: 2
    shell:
        """
        vt decompose {params.vt_decompose} {input.snvs} | vt normalize -n -r {input.reference} - -o {output.snvs}
        vt decompose {params.vt_decompose} {input.indels} | vt normalize -n -r {input.reference} - -o {output.indels}
        gzip -t {output.snvs}
        gzip -t {output.indels}
        tabix -f -p vcf {output.snvs}
        tabix -f -p vcf {output.indels}
        test -s {output.snvs_tbi}
        test -s {output.indels_tbi}
        """

rule mutect2_chromosome:
    output:
        vcf=temp(tmp_path("{tumour}.{chromosome}.mutect2.vcf.gz")),
        tbi=temp(tmp_path("{tumour}.{chromosome}.mutect2.vcf.gz.tbi")),
        stats=temp(tmp_path("{tumour}.{chromosome}.mutect2.vcf.gz.stats")),
        f1r2=temp(tmp_path("{tumour}.{chromosome}.mutect2.f1r2.tar.gz"))
    input:
        unpack(tumour_normal_crams),
        reference=config["reference"]["genome"],
        pon=panel_of_normals_dependency,
        pon_tbi=panel_of_normals_index_dependency,
        gnomad=config["reference"]["gnomad"],
        gnomad_tbi=lambda wildcards: f"{config['reference']['gnomad']}.tbi",
        shard_dependency=shard_dependency,
        manifest=rules.analysis_manifest.output
    log: "logs/mutect2/{tumour}.{chromosome}.log"
    benchmark: "results/benchmarks/mutect2/{tumour}.{chromosome}.tsv"
    params:
        normal=lambda wildcards: config["samples"]["tumours"][wildcards.tumour],
        tumour_bam_sample=lambda wildcards: bam_sample_name(wildcards.tumour),
        normal_bam_sample=lambda wildcards: bam_sample_name(config["samples"]["tumours"][wildcards.tumour]),
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        ref_container=container_path(config["reference"]["genome"]),
        core_container=lambda wildcards: container_path(shard_core_interval(wildcards)),
        shard_contig_arg=shard_contig_arg,
        gnomad_container=container_path(config["reference"]["gnomad"]),
        tumour_container=lambda wildcards: container_path(aligned_cram(wildcards.tumour)),
        normal_container=lambda wildcards: container_path(aligned_cram(config["samples"]["tumours"][wildcards.tumour])),
        raw_vcf=tmp_path("{tumour}.{chromosome}.mutect2.raw.publish.vcf.gz"),
        tmp_vcf=tmp_path("{tumour}.{chromosome}.mutect2.publish.vcf.gz"),
        interval_padding_arg=lambda wildcards: (
            f"--interval-padding {mutect2_interval_padding()}"
            if mutect2_interval_padding()
            else ""
        ),
        pon_arg=lambda wildcards: (
            f"--panel-of-normals {container_path(configured_panel_of_normals())}"
            if configured_panel_of_normals()
            else (
                f"--panel-of-normals /data/{panel_of_normals_output()}"
                if should_build_panel_of_normals()
                else ""
            )
        ),
        af_not_in_resource_arg=lambda wildcards: (
            f"--af-of-alleles-not-in-resource {config.get('mutect2', {}).get('af_of_alleles_not_in_resource')}"
            if config.get("mutect2", {}).get("af_of_alleles_not_in_resource") is not None
            else ""
        ),
        reference_args=mutect2_reference_args()
    threads: 4
    resources:
        mutect2_shards=1,
        mem_mb=16384,
        runtime=720
    shell:
        """
        mkdir -p $(dirname {log})
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' Mutect2 \
            -R {params.ref_container} \
            -I {params.tumour_container} \
            -I {params.normal_container} \
            --tumor-sample {params.tumour_bam_sample} \
            --normal-sample {params.normal_bam_sample} \
            {params.reference_args} \
            --germline-resource {params.gnomad_container} \
            {params.af_not_in_resource_arg} \
            {params.pon_arg} \
            -O /data/{params.raw_vcf} \
            -L {params.core_container} \
            {params.shard_contig_arg} \
            {params.interval_padding_arg} \
            --genotype-pon-sites false \
            --f1r2-tar-gz /data/{output.f1r2} > {log} 2>&1

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' SelectVariants \
            -R {params.ref_container} \
            -V /data/{params.raw_vcf} \
            -L {params.core_container} \
            {params.shard_contig_arg} \
            --interval-padding 0 \
            -O /data/{params.tmp_vcf} >> {log} 2>&1

        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        test -s {params.raw_vcf}.stats
        test -s {output.f1r2}
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        mv {params.raw_vcf}.stats {output.stats}
        rm -f {params.raw_vcf} {params.raw_vcf}.tbi
        """

rule mutect2_merge:
    output:
        vcf=protected("results/variants/{tumour}.mutect2.unfiltered.vcf.gz"),
        tbi=protected("results/variants/{tumour}.mutect2.unfiltered.vcf.gz.tbi"),
        stats=protected("results/variants/{tumour}.mutect2.unfiltered.vcf.gz.stats")
    input:
        vcfs=lambda wildcards: [tmp_path(f"{wildcards.tumour}.{chrom}.mutect2.vcf.gz") for chrom in SCATTER_IDS],
        tbis=lambda wildcards: [tmp_path(f"{wildcards.tumour}.{chrom}.mutect2.vcf.gz.tbi") for chrom in SCATTER_IDS],
        stats=lambda wildcards: [tmp_path(f"{wildcards.tumour}.{chrom}.mutect2.vcf.gz.stats") for chrom in SCATTER_IDS]
    log: "logs/mutect2/{tumour}.merge.log"
    benchmark: "results/benchmarks/mutect2/{tumour}.merge.tsv"
    params:
        inputs=lambda wildcards: ' '.join([f"-I /data/{tmp_path(f'{wildcards.tumour}.{chrom}.mutect2.vcf.gz')}" for chrom in SCATTER_IDS]),
        stats_inputs=lambda wildcards: ' '.join([f"--stats /data/{tmp_path(f'{wildcards.tumour}.{chrom}.mutect2.vcf.gz.stats')}" for chrom in SCATTER_IDS]),
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options(),
        tmp_vcf=tmp_path("{tumour}.mutect2.raw.publish.vcf.gz"),
        tmp_stats=tmp_path("{tumour}.mutect2.raw.publish.vcf.gz.stats")
    threads: 2
    shell:
        """
        mkdir -p $(dirname {log})
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' MergeVcfs \
            {params.inputs} \
            -O /data/{params.tmp_vcf} > {log} 2>&1

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' MergeMutectStats \
            {params.stats_inputs} \
            -O /data/{params.tmp_stats} >> {log} 2>&1

        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        test -s {params.tmp_stats}
        mkdir -p results/variants
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        mv {params.tmp_stats} {output.stats}
        """

rule mutect2_orientation_model:
    input:
        f1r2=lambda wildcards: [tmp_path(f"{wildcards.tumour}.{chrom}.mutect2.f1r2.tar.gz") for chrom in SCATTER_IDS]
    output:
        protected("results/variants/{tumour}.mutect2.orientation-model.tar.gz")
    log: "logs/mutect2/{tumour}.orientation-model.log"
    benchmark: "results/benchmarks/mutect2/{tumour}.orientation-model.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        inputs=lambda wildcards: " ".join(f"-I /data/{tmp_path(f'{wildcards.tumour}.{chrom}.mutect2.f1r2.tar.gz')}" for chrom in SCATTER_IDS),
        temporary=tmp_path("{tumour}.mutect2.orientation-model.tar.gz")
    threads: 1
    resources: mem_mb=8192, runtime=120, disk_mb=4096
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log}) results/variants
        {params.gatk_cmd} gatk --java-options '{params.java_options}' LearnReadOrientationModel \
          {params.inputs} -O /data/{params.temporary} > {log} 2>&1
        test -s {params.temporary}; mv {params.temporary} {output}
        """

rule contamination_sites:
    output:
        vcf=tmp_path("analysis", "contamination.sites.vcf.gz"),
        tbi=tmp_path("analysis", "contamination.sites.vcf.gz.tbi")
    input:
        source=lambda wildcards: config["reference"].get("contamination_sites") or config["reference"]["gnomad"],
        source_tbi=lambda wildcards: f"{config['reference'].get('contamination_sites') or config['reference']['gnomad']}.tbi",
        territory=lambda wildcards: analysis_territory(),
        territory_tbi=lambda wildcards: analysis_territory_index(),
        reference=config["reference"]["genome"]
    benchmark: "results/benchmarks/mutect2/contamination_sites.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("4g"),
        ref_container=container_path(config["reference"]["genome"]),
        source_container=lambda wildcards, input: container_path(input.source),
        territory_container=lambda wildcards: container_path(analysis_territory()),
        tmp_vcf=tmp_path("analysis", "contamination.sites.publish.vcf.gz")
    threads: 1
    resources: mem_mb=8192
    shell:
        """
        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' SelectVariants \
            -R {params.ref_container} \
            -V {params.source_container} \
            -L {params.territory_container} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            --sites-only-vcf-output true \
            -O /data/{params.tmp_vcf}
        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        """

rule mutect2_filter:
    output:
        vcf=temp(tmp_path("{tumour}.mutect2.filtered.vcf.gz")),
        tbi=temp(tmp_path("{tumour}.mutect2.filtered.vcf.gz.tbi")),
        contamination="results/qc/contamination/{tumour}.contamination.table",
        tumour_pileup="results/qc/contamination/{tumour}.organoid.pileup.table",
        baseline_pileup="results/qc/contamination/{tumour}.baseline.pileup.table"
    input:
        vcf="results/variants/{tumour}.mutect2.unfiltered.vcf.gz",
        tbi="results/variants/{tumour}.mutect2.unfiltered.vcf.gz.tbi",
        stats="results/variants/{tumour}.mutect2.unfiltered.vcf.gz.stats",
        orientation=rules.mutect2_orientation_model.output,
        cram=lambda wildcards: aligned_cram(wildcards.tumour),
        crai=lambda wildcards: aligned_cram_index(wildcards.tumour),
        normal_cram=lambda wildcards: aligned_cram(config["samples"]["tumours"][wildcards.tumour]),
        normal_crai=lambda wildcards: aligned_cram_index(config["samples"]["tumours"][wildcards.tumour]),
        reference=config["reference"]["genome"],
        contamination_sites=rules.contamination_sites.output.vcf,
        contamination_sites_tbi=rules.contamination_sites.output.tbi
    log: "logs/mutect2/{tumour}.filter.log"
    benchmark: "results/benchmarks/mutect2/{tumour}.filter.tsv"
    params:
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options("12g"),
        sites_container=container_path(rules.contamination_sites.output.vcf),
        ref_container=container_path(config["reference"]["genome"]),
        regions_container=container_path(config["reference"]["regions"]),
        cram_container=lambda wildcards: container_path(aligned_cram(wildcards.tumour)),
        normal_cram_container=lambda wildcards: container_path(aligned_cram(config["samples"]["tumours"][wildcards.tumour])),
        run_dir=tmp_path("{tumour}.mutect2_filter"),
        pileup_table=tmp_path("{tumour}.mutect2_filter", "{tumour}.pileup.table"),
        normal_pileup_table=tmp_path("{tumour}.mutect2_filter", "{tumour}.normal.pileup.table"),
        contamination_table=tmp_path("{tumour}.mutect2_filter", "{tumour}.contamination.table"),
        filtered_vcf=tmp_path("{tumour}.mutect2_filter", "{tumour}.filtered.vcf.gz")
    threads: 2
    shell:
        """
        rm -rf {params.run_dir}
        mkdir -p {params.run_dir} $(dirname {log})

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' GetPileupSummaries \
            -R {params.ref_container} \
            -I {params.cram_container} \
            -V {params.sites_container} \
            -L {params.sites_container} \
            -O /data/{params.pileup_table} > {log} 2>&1

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' GetPileupSummaries \
            -R {params.ref_container} \
            -I {params.normal_cram_container} \
            -V {params.sites_container} \
            -L {params.sites_container} \
            -O /data/{params.normal_pileup_table} >> {log} 2>&1

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' CalculateContamination \
            -I /data/{params.pileup_table} \
            --matched-normal /data/{params.normal_pileup_table} \
            -O /data/{params.contamination_table} >> {log} 2>&1

        python3 {PIPELINE_DIR}/workflow/scripts/analysis_mode.py validate-contamination \
            --pileup {params.pileup_table} \
            --contamination {params.contamination_table}

        {params.gatk_cmd} \
            gatk --java-options '{params.java_options}' FilterMutectCalls \
            -R {params.ref_container} \
            -V /data/{input.vcf} \
            --stats /data/{input.stats} \
            --contamination-table /data/{params.contamination_table} \
            --ob-priors /data/{input.orientation} \
            -O /data/{params.filtered_vcf} >> {log} 2>&1

        gzip -t {params.filtered_vcf}
        tabix -f -p vcf {params.filtered_vcf}
        test -s {params.filtered_vcf}.tbi
        test -s {params.contamination_table}
        mkdir -p results/qc/contamination
        mv {params.filtered_vcf} {output.vcf}
        mv {params.filtered_vcf}.tbi {output.tbi}
        mv {params.contamination_table} {output.contamination}
        mv {params.pileup_table} {output.tumour_pileup}
        mv {params.normal_pileup_table} {output.baseline_pileup}
        rm -rf {params.run_dir}
        """

rule mutect2_normalize:
    output:
        vcf="results/variants/{tumour}.mutect2.somatic.vcf.gz",
        tbi="results/variants/{tumour}.mutect2.somatic.vcf.gz.tbi"
    input:
        reference=config["reference"]["genome"],
        vcf=tmp_path("{tumour}.mutect2.filtered.vcf.gz"),
        tbi=tmp_path("{tumour}.mutect2.filtered.vcf.gz.tbi")
    benchmark: "results/benchmarks/mutect2/{tumour}.normalize.tsv"
    params:
        vt_decompose=config.get("vt", {}).get("decompose_params", "-s"),
        tmp_vcf=tmp_path("{tumour}.mutect2.somatic.vcf.gz")
    threads: 2
    shell:
        """
        vt decompose {params.vt_decompose} {input.vcf} | vt normalize -n -r {input.reference} - -o {params.tmp_vcf}
        gzip -t {params.tmp_vcf}
        tabix -f -p vcf {params.tmp_vcf}
        test -s {params.tmp_vcf}.tbi
        mkdir -p results/variants
        mv {params.tmp_vcf} {output.vcf}
        mv {params.tmp_vcf}.tbi {output.tbi}
        """

rule intersect_somatic_callers_raw:
    output:
        vcf=temp(tmp_path("{tumour}.intersect.raw.vcf.gz")),
        tbi=temp(tmp_path("{tumour}.intersect.raw.vcf.gz.tbi"))
    input:
        mutect2="results/variants/{tumour}.mutect2.somatic.vcf.gz",
        mutect2_tbi="results/variants/{tumour}.mutect2.somatic.vcf.gz.tbi",
        strelka_snvs=tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz"),
        strelka_snvs_tbi=tmp_path("{tumour}.strelka.somatic.snvs.af.norm.vcf.gz.tbi"),
        strelka_indels=tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz"),
        strelka_indels_tbi=tmp_path("{tumour}.strelka.somatic.indels.norm.vcf.gz.tbi")
    benchmark: "results/benchmarks/caller_tiers/{tumour}.raw_intersection.tsv"
    params:
        raw_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.intersect.raw.vcf")
    threads: 2
    shell:
        """
        python3 {PIPELINE_DIR}/scripts/vcf_intersect.py \
            --mutect2-vcf {input.mutect2} \
            --strelka-vcf {input.strelka_snvs} {input.strelka_indels} \
            --allowed-filters str_contraction LowDepth \
            --output-vcf {params.raw_vcf}

        bgzip -f {params.raw_vcf}
        gzip -t {output.vcf}
        tabix -p vcf {output.vcf}
        test -s {output.tbi}
        """

rule intersect_somatic_callers:
    output:
        vcf="results/variants/{tumour}.intersect.vcf.gz",
        tbi="results/variants/{tumour}.intersect.vcf.gz.tbi"
    input:
        vcf=lambda wildcards: final_vcf(wildcards.tumour) if sample_has_final_vcf(wildcards.tumour) else tmp_path(f"{wildcards.tumour}.intersect.raw.vcf.gz"),
        tbi=lambda wildcards: [] if sample_has_final_vcf(wildcards.tumour) else tmp_path(f"{wildcards.tumour}.intersect.raw.vcf.gz.tbi"),
        cram=lambda wildcards: [] if sample_has_final_vcf(wildcards.tumour) else aligned_cram(wildcards.tumour),
        crai=lambda wildcards: [] if sample_has_final_vcf(wildcards.tumour) else aligned_cram_index(wildcards.tumour),
        reference=config["reference"]["genome"],
        final_regions=lambda wildcards: [] if sample_has_final_vcf(wildcards.tumour) else rules.final_filter_regions.output.bed,
        final_regions_tbi=lambda wildcards: [] if sample_has_final_vcf(wildcards.tumour) else rules.final_filter_regions.output.tbi
    benchmark: "results/benchmarks/caller_tiers/{tumour}.filtered_intersection.tsv"
    params:
        final_vcf_mode=lambda wildcards: "true" if sample_has_final_vcf(wildcards.tumour) else "false",
        gatk_cmd=get_container_cmd(config["gatk"]),
        java_options=gatk_java_options(),
        ref_container=container_path(config["reference"]["genome"]),
        cram_container=lambda wildcards: "" if sample_has_final_vcf(wildcards.tumour) else container_path(aligned_cram(wildcards.tumour)),
        af=config["filtering"]["af_threshold"],
        tumour_dp=config["filtering"]["tumour_dp_threshold"],
        normal_dp=config["filtering"]["normal_dp_threshold"],
        normal=lambda wildcards: "" if sample_has_final_vcf(wildcards.tumour) else config["samples"]["tumours"][wildcards.tumour],
        tumour_bam_sample=lambda wildcards: "" if sample_has_final_vcf(wildcards.tumour) else bam_sample_name(wildcards.tumour),
        normal_bam_sample=lambda wildcards: "" if sample_has_final_vcf(wildcards.tumour) else bam_sample_name(config["samples"]["tumours"][wildcards.tumour]),
        vt_decompose=config.get("vt", {}).get("decompose_params", "-s"),
        depth_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.intersect.raw.depth.vcf.gz"),
        filtered_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.intersect.filtered.vcf"),
        region_filtered_vcf=lambda wildcards: tmp_path(f"{wildcards.tumour}.intersect.filtered.regions.vcf"),
        region_filtered_vcfgz=lambda wildcards: tmp_path(f"{wildcards.tumour}.intersect.filtered.regions.vcf.gz")
    threads: 2
    shell:
        """
        if [ "{params.final_vcf_mode}" = "true" ]; then
            vt decompose {params.vt_decompose} {input.vcf} | vt normalize -n -r {input.reference} - -o {params.region_filtered_vcfgz}
        else
            {params.gatk_cmd} \
                gatk --java-options '{params.java_options}' AnnotateVcfWithBamDepth \
                --lenient \
                -R {params.ref_container} \
                -I {params.cram_container} \
                -V /data/{input.vcf} \
                -O /data/{params.depth_vcf}

            python3 {PIPELINE_DIR}/scripts/filter_vcf.py \
                --input {params.depth_vcf} \
                --output {params.filtered_vcf} \
                --tumour-sample {params.tumour_bam_sample} \
                --normal-sample {params.normal_bam_sample} \
                --af {params.af} \
                --tumour-dp {params.tumour_dp} \
                --normal-dp {params.normal_dp}

            bedtools intersect \
                -header \
                -wa \
                -a {params.filtered_vcf} \
                -b {input.final_regions} \
                > {params.region_filtered_vcf}

            bgzip -f {params.region_filtered_vcf}
        fi

        gzip -t {params.region_filtered_vcfgz}
        tabix -f -p vcf {params.region_filtered_vcfgz}
        test -s {params.region_filtered_vcfgz}.tbi
        mkdir -p results/variants
        mv {params.region_filtered_vcfgz} {output.vcf}
        mv {params.region_filtered_vcfgz}.tbi {output.tbi}
        if [ "{params.final_vcf_mode}" != "true" ]; then
            rm -f {params.depth_vcf} {params.depth_vcf}.tbi
            rm -f {params.filtered_vcf}
        fi
        """
