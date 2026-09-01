import yaml
import os
import re
import sys
from pathlib import Path

PIPELINE_DIR = os.path.realpath(os.path.abspath(workflow.basedir))
configfile: os.path.join(PIPELINE_DIR, "config/config.yaml")

LOCAL_CONFIG_PATH = os.path.join(PIPELINE_DIR, "config/config.local.yaml")
if os.path.exists(LOCAL_CONFIG_PATH):
    configfile: LOCAL_CONFIG_PATH

sys.path.insert(0, os.path.join(PIPELINE_DIR, "workflow", "scripts"))
from analysis_mode import (
    analysis_settings,
    merge_adjacent,
    read_bed,
    read_fai,
    read_interval_list,
    reference_settings,
    shard_count,
    subtract_intervals,
    validate_existing_manifest,
    validate_intervals,
)
from workflow.scripts.sample_inputs import (
    aligned_cram as resolve_aligned_cram,
    aligned_cram_index as resolve_aligned_cram_index,
    sample_config as resolve_sample_config,
    sample_has_cram as resolve_sample_has_cram,
    sample_has_final_vcf as resolve_sample_has_final_vcf,
    sample_has_fastqs as resolve_sample_has_fastqs,
    bam_sample_name as resolve_bam_sample_name,
    final_vcf as resolve_final_vcf,
    final_vcf_samples as resolve_final_vcf_samples,
    matched_normal_samples as resolve_matched_normal_samples,
    comparison_map as resolve_comparison_map,
    vcf_only_mode as resolve_vcf_only_mode,
    validate_samples,
)
from annotation import validate_annotation_config
from coverage_qc import validate_exon_config
from somalier_qc import validate_somalier_config

storage_config = config.setdefault("storage", {})
LOCAL_TMPDIR = storage_config.setdefault("tmp_dir", "tmp")
LOCAL_SCRATCH = storage_config.setdefault("local_scratch", os.path.join(LOCAL_TMPDIR, "local_scratch"))

ANALYSIS = analysis_settings(config)
ANALYSIS_TYPE = ANALYSIS["type"]
REFERENCE = reference_settings(config)
REFERENCE_BUILD = REFERENCE["build"]
validate_existing_manifest(Path("results"), ANALYSIS_TYPE, REFERENCE_BUILD)

samples_path = config.get("run_management", {}).get("samples_file")
if not samples_path:
    samples_path = os.path.join(PIPELINE_DIR, "config/samples.yaml")
    if not os.path.exists(samples_path):
        samples_path = os.path.join(PIPELINE_DIR, "config/samples.example.yaml")

with open(samples_path) as f:
    samples = yaml.safe_load(f) or {}

config["samples"] = samples
validate_samples(config["samples"], samples_path)
# Internal compatibility mapping for the proven paired-caller rules. The public
# manifest remains role/comparison based and is the only accepted input schema.
config["samples"]["tumours"] = resolve_comparison_map(config["samples"], samples_path)

os.makedirs(LOCAL_TMPDIR, exist_ok=True)
os.makedirs(LOCAL_SCRATCH, exist_ok=True)


def tmp_path(*parts):
    if not parts:
        return LOCAL_TMPDIR
    return os.path.join(LOCAL_TMPDIR, *parts)

def scratch_path(*parts):
    return os.path.join(LOCAL_SCRATCH, *parts)

def bwa_index_inputs():
    reference = config["reference"]["genome"]
    suffix = config["reference"].get("bwa_index_suffix", "")
    extensions = ["amb", "ann", "bwt", "pac", "sa"]
    indexes = [f"{reference}{suffix}.{extension}" for extension in extensions]
    if REFERENCE_BUILD == "grch38":
        indexes.append(f"{reference}{suffix}.alt")
    return indexes

shell.executable("/bin/bash")
shell.prefix(
    f"set -euo pipefail; mkdir -p {LOCAL_TMPDIR} {LOCAL_SCRATCH} && "
    f"export TMPDIR=$(pwd)/{LOCAL_TMPDIR} && "
    f"export TMP=$(pwd)/{LOCAL_TMPDIR} && "
    f"export TEMP=$(pwd)/{LOCAL_TMPDIR} && "
)

# Load cluster config if it exists (for SLURM submission)
cluster_config_path = os.path.join(PIPELINE_DIR, "config/cluster.yaml")
if os.path.exists(cluster_config_path):
    with open(cluster_config_path) as f:
        cluster_config = yaml.safe_load(f)
else:
    cluster_config = {}

def configured_scatter_ids():
    if ANALYSIS_TYPE == "wes":
        return [str(value) for value in config["chromosomes"]]
    source = Path(config["reference"]["wgs_calling_regions"])
    fai = Path(f"{config['reference']['genome']}.fai")
    if not source.exists() or not fai.exists():
        missing = [str(path) for path in (source, fai) if not path.exists()]
        raise ValueError("WGS mode requires readable reference resources: " + ", ".join(missing))
    order, lengths = read_fai(fai)
    dictionary, intervals = read_interval_list(source)
    for contig in ANALYSIS["contigs"]:
        if dictionary.get(contig) != lengths.get(contig):
            raise ValueError(f"WGS interval-list/reference dictionary mismatch for {contig}")
    territory = validate_intervals(
        [value for value in intervals if value[0] in set(ANALYSIS["contigs"])], order, lengths
    )
    exclusion_path = config["reference"].get("wgs_exclude_regions")
    if exclusion_path:
        exclusions = merge_adjacent(sorted(
            read_bed(Path(exclusion_path)), key=lambda value: (order.index(value[0]), value[1], value[2])
        ))
        validate_intervals(exclusions, order, lengths)
        territory = validate_intervals(subtract_intervals(territory, exclusions), order, lengths)
    bases = sum(end - start for _, start, end in territory)
    count = shard_count(bases, ANALYSIS["target_bases_per_shard"], ANALYSIS["scatter_count"])
    return [f"{index:04d}" for index in range(count)]

SCATTER_IDS = configured_scatter_ids()

def analysis_territory():
    if ANALYSIS_TYPE == "wgs":
        return tmp_path("analysis", "wgs.callable.bed.gz")
    return config["reference"]["regions"]

def analysis_territory_index():
    return f"{analysis_territory()}.tbi"

def strelka_mode_params():
    configured = str(config.get("strelka", {}).get("params", "")).strip()
    pieces = [configured] if configured else []
    if ANALYSIS_TYPE == "wes" and "--exome" not in configured.split():
        pieces.append("--exome")
    if ANALYSIS_TYPE == "wgs" and "--exome" in configured.split():
        raise ValueError("strelka.params must not contain --exome in WGS mode")
    return " ".join(pieces)

def sample_config(sample):
    return resolve_sample_config(config["samples"], sample, samples_path)

def sample_has_fastqs(sample):
    return resolve_sample_has_fastqs(config["samples"], sample, samples_path)

def sample_has_cram(sample):
    return resolve_sample_has_cram(config["samples"], sample, samples_path)

def sample_has_final_vcf(sample):
    return resolve_sample_has_final_vcf(config["samples"], sample, samples_path)

VCF_ONLY_MODE = resolve_vcf_only_mode(config["samples"], samples_path)
ANNOTATION_RESOURCES = validate_annotation_config(config, REFERENCE_BUILD)
ANNOTATION_ENABLED = bool(config.get("annotation", {}).get("enabled", False))
validate_exon_config(config, REFERENCE_BUILD)
validate_somalier_config(config, REFERENCE_BUILD)
signature_build = str(config.get("mutational_signatures", {}).get(
    "reference_build", "grch37"
)).lower()
if signature_build != REFERENCE_BUILD:
    raise ValueError(
        f"mutational_signatures.reference_build is {signature_build!r}, expected {REFERENCE_BUILD!r}"
    )
hotspot_build = str(config.get("hotspots", {}).get("reference_build", "grch37")).lower()
if hotspot_build != REFERENCE_BUILD:
    raise ValueError(
        f"hotspots.reference_build is {hotspot_build!r}, expected {REFERENCE_BUILD!r}"
    )

def external_cram_bindings():
    bind_dirs = []
    seen = set()
    for sample_data in config["samples"]["samples"].values():
        for key in ("cram", "crai"):
            path = sample_data.get(key)
            if path and os.path.isabs(path):
                bind_dir = os.path.dirname(path)
                if bind_dir and bind_dir not in seen:
                    seen.add(bind_dir)
                    bind_dirs.append(bind_dir)
    return [
        (bind_dir, f"/external_inputs/{index:03d}")
        for index, bind_dir in enumerate(bind_dirs)
    ]

def external_resource_bindings():
    paths = []
    for value in config.get("reference", {}).values():
        if isinstance(value, str):
            paths.append(value)
    for resource in config.get("annotation", {}).get("resources", {}).values():
        if isinstance(resource, dict):
            paths.append(resource.get("path"))
            if isinstance(resource.get("index"), dict):
                paths.append(resource["index"].get("path"))
    coverage_resource = config.get("coverage", {}).get("exon_bed", {})
    paths.append(coverage_resource.get("path") if isinstance(coverage_resource, dict) else None)
    somalier = config.get("somalier", {})
    paths.extend((somalier.get("sites_vcf"), somalier.get("sites_vcf_index")))
    ancestry = somalier.get("ancestry", {})
    paths.extend((ancestry.get("labels"), ancestry.get("reference_somalier_dir")))
    directories = []
    for value in paths:
        if not value or not os.path.isabs(value):
            continue
        path = Path(value)
        directory = str(path if path.is_dir() else path.parent)
        if directory not in directories:
            directories.append(directory)
    return [(directory, f"/external_resources/{index:03d}") for index, directory in enumerate(directories)]

def extra_container_bindings():
    bindings = external_cram_bindings() + external_resource_bindings()
    bound_hosts = {host for host, _ in bindings}
    if PIPELINE_DIR not in bound_hosts:
        bindings.append((PIPELINE_DIR, "/pipeline"))
        bound_hosts.add(PIPELINE_DIR)
    if os.path.isabs(LOCAL_SCRATCH) and LOCAL_SCRATCH not in bound_hosts:
        bindings.append((LOCAL_SCRATCH, "/local_scratch"))
    return bindings

def get_container_image(tool_config):
    """Return the configured image for the selected runtime."""
    runtime = config.get("container_runtime", "singularity")
    if runtime in {"apptainer", "singularity"}:
        return tool_config.get("singularity_image", tool_config.get("image"))
    if runtime == "docker":
        image = tool_config.get("docker_image")
        if image:
            return image
        legacy_image = tool_config.get("image")
        if legacy_image and not legacy_image.endswith(".sif"):
            return legacy_image
        raise ValueError("Docker runtime requires a docker_image value in tool config")
    raise ValueError(f"Unknown container_runtime: {runtime}")

def get_container_cmd(tool_config, bind_node_tmp=False):
    """Generate container command based on configured runtime."""
    runtime = config.get("container_runtime", "singularity")
    image_path = get_container_image(tool_config)
    if not image_path:
        raise ValueError(f"No container image configured for runtime: {runtime}")
    tmp_env = f"env TMPDIR={CONTAINER_TMPDIR} TMP={CONTAINER_TMPDIR} TEMP={CONTAINER_TMPDIR}"
    external_bindings = extra_container_bindings()
    if bind_node_tmp:
        external_bindings = [binding for binding in external_bindings if binding[1] != "/local_scratch"]
        external_bindings.append(("/tmp", "/local_scratch"))

    if runtime in {"apptainer", "singularity"}:
        # For relative SIF paths, prefix with $(pwd) so they resolve at execution time.
        if not os.path.isabs(image_path):
            image_path = f"$(pwd)/{image_path}"
        bind_arg = ",".join(
            ["$(pwd):/data"]
            + [f"{host}:{container}" for host, container in external_bindings]
        )
        return f"{runtime} exec --bind {bind_arg} --pwd /data {image_path} {tmp_env}"
    elif runtime == "docker":
        bind_args = " ".join(
            ["-v $(pwd):/data"]
            + [f"-v {host}:{container}" for host, container in external_bindings]
        )
        return f"docker run --rm --platform linux/amd64 {bind_args} -w /data {image_path} {tmp_env}"
    else:
        raise ValueError(f"Unknown container_runtime: {runtime}")

def gatk_java_options(heap=None):
    options = [f"-Djava.io.tmpdir={CONTAINER_TMPDIR}"]
    if heap:
        options.insert(0, f"-Xmx{heap}")
    return " ".join(options)

def container_path(path):
    """Convert host path to container path."""
    if os.path.isabs(path):
        absolute_path = os.path.abspath(path)
        pipeline_dir = os.path.abspath(os.getcwd())
        try:
            if os.path.commonpath([absolute_path, pipeline_dir]) == pipeline_dir:
                return f"/data/{os.path.relpath(absolute_path, pipeline_dir)}"
        except ValueError:
            pass

        for host_dir, container_dir in extra_container_bindings():
            absolute_host_dir = os.path.abspath(host_dir)
            try:
                if os.path.commonpath([absolute_path, absolute_host_dir]) == absolute_host_dir:
                    relative_path = os.path.relpath(absolute_path, absolute_host_dir)
                    if relative_path == ".":
                        return container_dir
                    return f"{container_dir}/{relative_path}"
            except ValueError:
                continue

        # Other absolute paths must already be visible inside the container.
        return path
    else:
        # Relative path - convert to /data-relative
        return f"/data/{path}"


CONTAINER_TMPDIR = container_path(LOCAL_TMPDIR)

def aligned_cram(sample):
    return resolve_aligned_cram(config["samples"], sample, samples_path)

def aligned_cram_index(sample):
    return resolve_aligned_cram_index(config["samples"], sample, samples_path)

def bam_sample_name(sample):
    return resolve_bam_sample_name(config["samples"], sample, samples_path)

def final_vcf(sample):
    return resolve_final_vcf(config["samples"], sample, samples_path)

def tumour_normal_crams(wildcards):
    tumour = wildcards.tumour
    normal = config["samples"]["tumours"][tumour]
    return {
        "tumour": aligned_cram(tumour),
        "normal": aligned_cram(normal),
        "tumour_index": aligned_cram_index(tumour),
        "normal_index": aligned_cram_index(normal),
        "tumour_validation": f"results/qc/preflight/{tumour}.alignment.json",
        "normal_validation": f"results/qc/preflight/{normal}.alignment.json",
    }

def get_fastqs(sample):
    s = sample_config(sample)
    return [s["fastq_1"], s["fastq_2"]]

def normal_samples():
    if VCF_ONLY_MODE:
        return []
    all_samples = set(config["samples"]["samples"])
    tumours = set(config["samples"]["tumours"])
    return sorted(all_samples - tumours)

def matched_normal_samples():
    if not config.get("germline", {}).get("enabled", True):
        return []
    return resolve_matched_normal_samples(config["samples"], samples_path)

def tumour_samples():
    if VCF_ONLY_MODE:
        return resolve_final_vcf_samples(config["samples"], samples_path)
    return sorted(config["samples"]["tumours"])

def wildcard_regex(values):
    values = sorted(str(value) for value in values)
    if not values:
        return "(?!)"
    return "|".join(re.escape(value) for value in values)

wildcard_constraints:
    sample=wildcard_regex(config["samples"]["samples"]),
    tumour=wildcard_regex(tumour_samples()),
    normal=wildcard_regex(set(normal_samples()) | set(matched_normal_samples())),
    chromosome=wildcard_regex(SCATTER_IDS)

def fastq_samples():
    if VCF_ONLY_MODE:
        return []
    return sorted(sample for sample in config["samples"]["samples"] if sample_has_fastqs(sample))

def pipeline_aligned_samples():
    if VCF_ONLY_MODE:
        return []
    return sorted(sample for sample in config["samples"]["samples"] if not sample_has_cram(sample))

def tissue_signature_outputs(signature_type):
    return expand(
        f"results/signatures/exposures/{{tumour}}.intersect.tissue.{signature_type}.tsv",
        tumour=tumour_samples(),
    )


def aggregate_tissue_signature_outputs():
    outputs = [
        "results/aggregate/signatures.tissue.sbs.tsv",
        "results/aggregate/signatures.tissue.id.tsv",
        "results/aggregate/signatures.tissue.dbs.tsv",
    ]
    if not VCF_ONLY_MODE:
        outputs.insert(2, "results/aggregate/signatures.strelka.tissue.id.tsv")
    return outputs


def extended_contexts_output():
    return "results/aggregate/extended_contexts.tsv"


def hotspots_output():
    return "results/aggregate/hotspots.tsv"


def fastqc_outputs():
    return expand(
        "results/qc/fastqc/{sample}/fastqc_done",
        sample=fastq_samples(),
    )


def variant_outputs():
    return expand(
        "results/catalogs/{tumour}.stringent.vcf.gz",
        tumour=tumour_samples(),
    )

def somatic_annotation_source(sample):
    return final_vcf(sample) if VCF_ONLY_MODE else f"results/variants/{sample}.intersect.vcf.gz"

def somatic_annotation_outputs():
    if not ANNOTATION_ENABLED:
        return []
    patterns = (
        "results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz",
        "results/annotations/somatic/{tumour}.intersect.annotated.vcf.gz.tbi",
    )
    return [path for pattern in patterns for path in expand(pattern, tumour=tumour_samples())]

def germline_annotation_outputs():
    if not ANNOTATION_ENABLED or VCF_ONLY_MODE:
        return []
    patterns = (
        "results/annotations/germline/{normal}.haplotypecaller.filtered.annotated.vcf.gz",
        "results/annotations/germline/{normal}.haplotypecaller.filtered.annotated.vcf.gz.tbi",
    )
    return [path for pattern in patterns for path in expand(pattern, normal=matched_normal_samples())]

def somalier_outputs():
    if VCF_ONLY_MODE or not config.get("somalier", {}).get("enabled", False):
        return []
    outputs = ["results/qc/somalier/somalier_flags.tsv"]
    if config.get("somalier", {}).get("ancestry", {}).get("enabled", False):
        outputs.extend(("results/qc/somalier/ancestry.tsv", "results/qc/somalier/ancestry.html"))
    return outputs

def germline_outputs():
    patterns = (
        "results/germline/{normal}.haplotypecaller.g.vcf.gz",
        "results/germline/{normal}.haplotypecaller.g.vcf.gz.tbi",
        "results/germline/{normal}.haplotypecaller.filtered.vcf.gz",
        "results/germline/{normal}.haplotypecaller.filtered.vcf.gz.tbi",
    )
    return [
        path
        for pattern in patterns
        for path in expand(pattern, normal=matched_normal_samples())
    ]

include: "workflow/rules/fastqc.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/validation.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/organoid_catalogs.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/somalier.smk"
include: "workflow/rules/provenance.smk"

rule variants_done:
    input:
        variant_outputs()
    output:
        touch("results/variants/variants.done")


rule qc_done:
    input:
        "results/qc/variants/final_variant_counts_mqc.tsv" if VCF_ONLY_MODE else "results/aggregate/qc_summary.html",
        somalier_outputs()
    output:
        touch("results/aggregate/qc.done")


rule provenance_done:
    input:
        rules.provenance.output.tsv
    output:
        touch("results/aggregate/provenance.done")


rule all:
    default_target: True
    input:
        rules.analysis_manifest.output,
        expand("results/catalogs/{tumour}.stringent.vcf.gz", tumour=tumour_samples()),
        expand("results/catalogs/{tumour}.stringent.tsv", tumour=tumour_samples()),
        expand("results/signatures/{tumour}.sbs96.tsv", tumour=tumour_samples()),
        rules.qc_done.output,
        rules.provenance_done.output
