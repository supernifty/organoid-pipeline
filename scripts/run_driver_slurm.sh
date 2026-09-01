#!/usr/bin/env bash
# Submit a persistent batch controller as one SLURM job.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SETUP_SCRIPT="$SCRIPT_DIR/setup.sh"
PYTHON="$PIPELINE_DIR/.pixi/envs/default/bin/python"

DRIVER_DRY_RUN=false
BATCH=""
SAMPLES=""
RESUME=false
RUN_ARGS=()
TARGETS=()
ORIGINAL_ARGS=("$@")
SNAKEMAKE_DRY_RUN=false
while [ "$#" -gt 0 ]; do
    case "$1" in
        --driver-dry-run) DRIVER_DRY_RUN=true; shift ;;
        --batch) BATCH="${2:?--batch requires a value}"; shift 2 ;;
        --samples) SAMPLES="${2:?--samples requires a value}"; shift 2 ;;
        --resume) RESUME=true; shift ;;
        -n|--dryrun|--dry-run)
            SNAKEMAKE_DRY_RUN=true
            RUN_ARGS+=("$1")
            shift
            ;;
        --)
            shift
            while [ "$#" -gt 0 ]; do TARGETS+=("$1"); shift; done
            ;;
        *)
            if [[ "$1" == results/* ]]; then
                TARGETS+=("$1")
            else
                RUN_ARGS+=("$1")
            fi
            shift
            ;;
    esac
done

if [ -z "$BATCH" ]; then
    echo "ERROR: --batch is required." >&2
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
elif [ "$DRIVER_DRY_RUN" = true ]; then
    echo "WARNING: setup.sh not found; driver dry run will proceed."
else
    echo "ERROR: setup.sh not found; copy scripts/setup.sh.example and customize it." >&2
    exit 1
fi

PREPARE=("$PYTHON" "$SCRIPT_DIR/run_manager.py" prepare --batch "$BATCH")
if [ -n "$SAMPLES" ]; then PREPARE+=(--samples "$SAMPLES"); fi
if [ "$RESUME" = true ]; then PREPARE+=(--resume); fi
if [ "$SNAKEMAKE_DRY_RUN" = true ]; then PREPARE+=(--dry-run); fi
for target in "${TARGETS[@]+"${TARGETS[@]}"}"; do PREPARE+=(--target "$target"); done
PREPARE+=(--command "$0" "${ORIGINAL_ARGS[@]}")
PREPARED_JSON="$("${PREPARE[@]}")"
LAUNCH="$($PYTHON -c 'import json,sys; print(json.load(sys.stdin)["launch"])' <<<"$PREPARED_JSON")"
RUN_DIR="$($PYTHON -c 'import json,sys; print(json.load(sys.stdin)["run_dir"])' <<<"$PREPARED_JSON")"

JOBSCRIPT_DIR="$PIPELINE_DIR/tmp/codex/run-driver"
mkdir -p "$JOBSCRIPT_DIR"
JOBSCRIPT="$(mktemp "$JOBSCRIPT_DIR/${BATCH}.XXXXXX")"
cleanup() { rm -f "$JOBSCRIPT"; }
trap cleanup EXIT

printf '#!/usr/bin/env bash\nset -euo pipefail\ncd %q\nexec ./scripts/run_cluster.sh --batch %q --prepared-launch %q' \
    "$PIPELINE_DIR" "$BATCH" "$LAUNCH" >"$JOBSCRIPT"
for arg in "${RUN_ARGS[@]+"${RUN_ARGS[@]}"}"; do printf ' %q' "$arg" >>"$JOBSCRIPT"; done
if [ "${#TARGETS[@]}" -gt 0 ]; then
    for target in "${TARGETS[@]}"; do printf ' %q' "$target" >>"$JOBSCRIPT"; done
fi
printf '\n' >>"$JOBSCRIPT"
chmod +x "$JOBSCRIPT"

DRIVER_MEM="${SNAKEMAKE_DRIVER_MEM_MB:-32768}"
DRIVER_RUNTIME="${SNAKEMAKE_DRIVER_RUNTIME_MIN:-1440}"
DRIVER_CPUS="${SNAKEMAKE_DRIVER_CPUS:-1}"
DRIVER_JOB_NAME="${SNAKEMAKE_DRIVER_JOB_NAME:-snakemake-$BATCH}"
SBATCH_CMD=(sbatch --job-name "$DRIVER_JOB_NAME" --cpus-per-task "$DRIVER_CPUS" --mem "$DRIVER_MEM" --time "$DRIVER_RUNTIME")
SBATCH_CMD+=(--output "$RUN_DIR/log/slurm-driver-%j.out" --error "$RUN_DIR/log/slurm-driver-%j.out")
if [ -n "${SLURM_ACCOUNT:-}" ]; then SBATCH_CMD+=(-A "$SLURM_ACCOUNT"); fi
if [ -n "${SLURM_PARTITION:-}" ]; then SBATCH_CMD+=(-p "$SLURM_PARTITION"); fi
if [ -n "${SLURM_EXTRA:-}" ]; then
    # shellcheck disable=SC2206
    EXTRA=( ${SLURM_EXTRA} )
    SBATCH_CMD+=("${EXTRA[@]}")
fi
SBATCH_CMD+=("$JOBSCRIPT")

echo "=== Batch controller submission: $BATCH ==="
echo "Run directory: $RUN_DIR"
echo "Launch: $LAUNCH"
if [ "$DRIVER_DRY_RUN" = true ]; then
    echo "Generated jobscript:"
    sed -n '1,80p' "$JOBSCRIPT"
    printf 'sbatch command:'; printf ' %q' "${SBATCH_CMD[@]}"; printf '\n'
    exit 0
fi

"$PYTHON" "$SCRIPT_DIR/run_manager.py" transition --batch "$BATCH" --launch "$LAUNCH" --state submitted
set +e
SBATCH_OUTPUT="$("${SBATCH_CMD[@]}" 2>&1)"
SBATCH_STATUS=$?
set -e
if [ "$SBATCH_STATUS" -ne 0 ]; then
    "$PYTHON" "$SCRIPT_DIR/run_manager.py" transition --batch "$BATCH" --launch "$LAUNCH" --state failed --error "$SBATCH_OUTPUT"
    echo "$SBATCH_OUTPUT" >&2
    exit "$SBATCH_STATUS"
fi
echo "$SBATCH_OUTPUT"
if [[ "$SBATCH_OUTPUT" =~ Submitted[[:space:]]+batch[[:space:]]+job[[:space:]]+([0-9]+) ]]; then
    JOB_ID="${BASH_REMATCH[1]}"
else
    "$PYTHON" "$SCRIPT_DIR/run_manager.py" transition --batch "$BATCH" --launch "$LAUNCH" --state failed --error "could not parse sbatch job ID: $SBATCH_OUTPUT"
    echo "ERROR: could not parse sbatch job ID." >&2
    exit 1
fi
"$PYTHON" "$SCRIPT_DIR/run_manager.py" transition --batch "$BATCH" --launch "$LAUNCH" --state submitted --job-id "$JOB_ID"
echo "Controller log: $RUN_DIR/log/slurm-driver-$JOB_ID.out"
