#!/usr/bin/env bash
# Run one isolated batch through the Snakemake SLURM profile.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PROFILE="$PIPELINE_DIR/config/slurm"
SETUP_SCRIPT="$SCRIPT_DIR/setup.sh"
PYTHON="$PIPELINE_DIR/.pixi/envs/default/bin/python"

BATCH=""
SAMPLES=""
RESUME=false
PREPARED_LAUNCH=""
IS_DRYRUN=false
SNAKEMAKE_ARGS=()
TARGETS=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        --batch) BATCH="${2:?--batch requires a value}"; shift 2 ;;
        --samples) SAMPLES="${2:?--samples requires a value}"; shift 2 ;;
        --resume) RESUME=true; shift ;;
        --prepared-launch) PREPARED_LAUNCH="${2:?--prepared-launch requires a value}"; shift 2 ;;
        -n|--dryrun|--dry-run) IS_DRYRUN=true; SNAKEMAKE_ARGS+=("$1"); shift ;;
        --)
            shift
            while [ "$#" -gt 0 ]; do TARGETS+=("$1"); shift; done
            ;;
        *)
            if [[ "$1" == results/* ]]; then
                TARGETS+=("$1")
            else
                SNAKEMAKE_ARGS+=("$1")
            fi
            shift
            ;;
    esac
done

if [ -z "$BATCH" ]; then
    echo "ERROR: --batch is required. Use snakemake directly for development runs." >&2
    exit 2
fi
if [ ! -x "$PYTHON" ]; then
    echo "ERROR: Pixi environment is missing; run 'pixi install'." >&2
    exit 1
fi

source_setup_script() {
    local restore_nounset=0
    if [[ $- == *u* ]]; then restore_nounset=1; set +u; fi
    # shellcheck disable=SC1090
    source "$SETUP_SCRIPT"
    if [ "$restore_nounset" -eq 1 ]; then set -u; fi
}

if [ -f "$SETUP_SCRIPT" ]; then
    echo "Loading setup script: $SETUP_SCRIPT"
    source_setup_script
elif [ "$IS_DRYRUN" = true ]; then
    echo "WARNING: setup.sh not found; dry run will proceed."
else
    echo "ERROR: setup.sh not found; copy scripts/setup.sh.example and customize it." >&2
    exit 1
fi

if [ -n "$PREPARED_LAUNCH" ]; then
    LAUNCH="$PREPARED_LAUNCH"
    RUN_DIR="$PIPELINE_DIR/runs/$BATCH"
    CONFIG="$RUN_DIR/config/current/config.yaml"
else
    PREPARE=("$PYTHON" "$SCRIPT_DIR/run_manager.py" prepare --batch "$BATCH")
    if [ -n "$SAMPLES" ]; then PREPARE+=(--samples "$SAMPLES"); fi
    if [ "$RESUME" = true ]; then PREPARE+=(--resume); fi
    if [ "$IS_DRYRUN" = true ]; then PREPARE+=(--dry-run); fi
    for target in "${TARGETS[@]+"${TARGETS[@]}"}"; do PREPARE+=(--target "$target"); done
    PREPARE+=(--command "$0" --batch "$BATCH")
    if [ -n "$SAMPLES" ]; then PREPARE+=(--samples "$SAMPLES"); fi
    if [ "$RESUME" = true ]; then PREPARE+=(--resume); fi
    PREPARE+=("${SNAKEMAKE_ARGS[@]}")
    PREPARED_JSON="$("${PREPARE[@]}")"
    LAUNCH="$($PYTHON -c 'import json,sys; print(json.load(sys.stdin)["launch"])' <<<"$PREPARED_JSON")"
    RUN_DIR="$($PYTHON -c 'import json,sys; print(json.load(sys.stdin)["run_dir"])' <<<"$PREPARED_JSON")"
    CONFIG="$($PYTHON -c 'import json,sys; print(json.load(sys.stdin)["config"])' <<<"$PREPARED_JSON")"
fi

export XDG_CACHE_HOME="$RUN_DIR/tmp/xdg-cache"
mkdir -p "$XDG_CACHE_HOME"

CONFIGURED_MUTECT2_SHARDS="$($PYTHON -c 'import sys,yaml; c=yaml.safe_load(open(sys.argv[1])); print(c.get("analysis",{}).get("wgs",{}).get("max_concurrent_mutect2_shards",32))' "$CONFIG")"
MAX_MUTECT2_SHARDS="${MUTECT2_MAX_CONCURRENT_SHARDS:-$CONFIGURED_MUTECT2_SHARDS}"
CONFIGURED_HAPLOTYPECALLER_SHARDS="$($PYTHON -c 'import sys,yaml; c=yaml.safe_load(open(sys.argv[1])); print(c.get("germline",{}).get("max_concurrent_haplotypecaller_shards",16))' "$CONFIG")"
MAX_HAPLOTYPECALLER_SHARDS="${HAPLOTYPECALLER_MAX_CONCURRENT_SHARDS:-$CONFIGURED_HAPLOTYPECALLER_SHARDS}"
if ! [[ "$MAX_MUTECT2_SHARDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: Mutect2 shard concurrency must be a positive integer: $MAX_MUTECT2_SHARDS" >&2
    exit 1
fi
if ! [[ "$MAX_HAPLOTYPECALLER_SHARDS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: HaplotypeCaller shard concurrency must be a positive integer: $MAX_HAPLOTYPECALLER_SHARDS" >&2
    exit 1
fi

CMD=("$PIPELINE_DIR/.pixi/envs/default/bin/snakemake" --snakefile "$PIPELINE_DIR/Snakefile" --directory "$RUN_DIR" --configfile "$CONFIG")
CMD+=(--profile "$PROFILE" --jobs "${SNAKEMAKE_MAX_JOBS:-256}" --rerun-incomplete)
CMD+=(--resources "mutect2_shards=$MAX_MUTECT2_SHARDS" "haplotypecaller_shards=$MAX_HAPLOTYPECALLER_SHARDS")
DEFAULT_RESOURCES=("runtime=60" "mem_mb=4096")
if [ -n "${SLURM_ACCOUNT:-}" ]; then DEFAULT_RESOURCES+=("slurm_account=$SLURM_ACCOUNT"); fi
if [ -n "${SLURM_PARTITION:-}" ]; then DEFAULT_RESOURCES+=("slurm_partition=$SLURM_PARTITION"); fi
if [ -n "${SLURM_EXTRA:-}" ]; then DEFAULT_RESOURCES+=("slurm_extra=$SLURM_EXTRA"); fi
CMD+=(--default-resources "${DEFAULT_RESOURCES[@]}")
CMD+=("${SNAKEMAKE_ARGS[@]}")
if [ "${#TARGETS[@]}" -gt 0 ]; then CMD+=(-- "${TARGETS[@]}"); fi

echo "=== Somatic Variant Calling Pipeline: $BATCH ==="
echo "Run directory: $RUN_DIR"
echo "Launch: $LAUNCH"
echo "Max jobs: ${SNAKEMAKE_MAX_JOBS:-256}"
echo "Max concurrent Mutect2/PoN shards: $MAX_MUTECT2_SHARDS"
echo "Max concurrent HaplotypeCaller shards: $MAX_HAPLOTYPECALLER_SHARDS"

exec "$PYTHON" "$SCRIPT_DIR/run_manager.py" controller --batch "$BATCH" --launch "$LAUNCH" "${CMD[@]}"
