#!/usr/bin/env bash
#SBATCH --account=cnms
#SBATCH --partition=high_mem
#SBATCH --qos=cnms-high_mem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --signal=B:USR1@600

set -eo pipefail

source /etc/profile
set -u

: "${CAMPAIGN_ROOT:?CAMPAIGN_ROOT must be exported by submit_wave.py}"
: "${TWIST_ID:?TWIST_ID must be exported by submit_wave.py}"
: "${PHIX:?PHIX must be exported by submit_wave.py}"
: "${PHIY:?PHIY must be exported by submit_wave.py}"
: "${SEED:?SEED must be exported by submit_wave.py}"
: "${TARGET_R:?TARGET_R must be exported by submit_wave.py}"

THREADS=${FTLM_THREADS:-16}
MAX_RUNTIME_MINUTES=${MAX_RUNTIME_MINUTES:-210}
LANCZOS_MAX_STEPS=${LANCZOS_MAX_STEPS:-80}
LANCZOS_SAVE_STEPS=${LANCZOS_SAVE_STEPS:-}
# Slurm --export uses commas as separators, so submit_wave transports the list
# with ':' and the wrapper restores the CLI representation.
LANCZOS_SAVE_STEPS=${LANCZOS_SAVE_STEPS//:/,}
CHECKPOINT_SERIES=${CHECKPOINT_SERIES:-legacy}
BIN_DIR=${BIN_DIR:-"$CAMPAIGN_ROOT/bin"}
RUN_DIR="$CAMPAIGN_ROOT/runs/twist_$TWIST_ID"
TARGET_TAG=$(printf '%03d' "$TARGET_R")
M_TAG=$(printf '%03d' "$LANCZOS_MAX_STEPS")
if [[ $CHECKPOINT_SERIES == legacy ]]; then
  CHECKPOINT="$RUN_DIR/twist_${TWIST_ID}.ftlmcp"
  PROGRESS_JSONL="$RUN_DIR/progress.jsonl"
  PROGRESS_STATE="$RUN_DIR/progress_state.json"
  OUTPUT="$RUN_DIR/twist_${TWIST_ID}_thermo_R${TARGET_TAG}.csv"
  FILE_TAG="R${TARGET_TAG}"
else
  if [[ ! $CHECKPOINT_SERIES =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "invalid CHECKPOINT_SERIES=$CHECKPOINT_SERIES" >&2
    exit 2
  fi
  CHECKPOINT="$RUN_DIR/twist_${TWIST_ID}_${CHECKPOINT_SERIES}.ftlmcp"
  PROGRESS_JSONL="$RUN_DIR/progress_${CHECKPOINT_SERIES}.jsonl"
  PROGRESS_STATE="$RUN_DIR/progress_state_${CHECKPOINT_SERIES}.json"
  OUTPUT="$RUN_DIR/twist_${TWIST_ID}_thermo_m${M_TAG}_R${TARGET_TAG}.csv"
  FILE_TAG="${CHECKPOINT_SERIES}_m${M_TAG}_R${TARGET_TAG}"
fi

# Resolve the deployment symlinks once.  A later deployment may update
# CAMPAIGN_ROOT/bin while this job is running; the active FTLM and reducer must
# nevertheless come from one immutable build.
FTLM_EXEC=$(readlink -f "$BIN_DIR/ftlm_n_vs_mu")
REDUCER_EXEC=$(readlink -f "$BIN_DIR/ftlm_reduce_checkpoint")
IMMUTABLE_BIN_DIR=$(dirname "$FTLM_EXEC")
if [[ $(dirname "$REDUCER_EXEC") != "$IMMUTABLE_BIN_DIR" ]]; then
  echo "FTLM and reducer resolved to different deployments" >&2
  exit 2
fi
DEPLOYED_COMMIT=$(cat "$IMMUTABLE_BIN_DIR/GIT_COMMIT")
SOURCE_DIR="$CAMPAIGN_ROOT/sources/$DEPLOYED_COMMIT"
FTLM_SHA256=$(sha256sum "$FTLM_EXEC" | awk '{print $1}')
REDUCER_SHA256=$(sha256sum "$REDUCER_EXEC" | awk '{print $1}')

module purge
module load intel/2024.1.0 mkl/2024.1.0
export OMP_NUM_THREADS="$THREADS"
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=false
export MALLOC_ARENA_MAX=4

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ulimit -c 0

COMMAND=(
  "$FTLM_EXEC"
  --lx 4 --ly 4
  --tx 1 --ty 1 --tp 0
  --u 7 --v 0
  --phix "$PHIX" --phiy "$PHIY"
  --beta-list 2.857142857142857,12.5
  --mu-min -3 --mu-max 4 --mu-count 281
  --samples "$TARGET_R"
  --lanczos-steps "$LANCZOS_MAX_STEPS"
  --exact-block-threshold 256
  --threads "$THREADS"
  --seed "$SEED"
  --max-sector-dim 170000000
  --checkpoint "$CHECKPOINT"
  --max-runtime-minutes "$MAX_RUNTIME_MINUTES"
  --progress-jsonl "$PROGRESS_JSONL"
  --progress-state "$PROGRESS_STATE"
  --progress-interval-seconds 60
  --twist-id "$TWIST_ID"
  --output "$OUTPUT"
)
if [[ -n $LANCZOS_SAVE_STEPS ]]; then
  COMMAND+=(--lanczos-save-steps "$LANCZOS_SAVE_STEPS")
fi
if [[ -n ${ONLY_BLOCK_NUP:-} ]]; then
  ONLY_BLOCK="${ONLY_BLOCK_NUP},${ONLY_BLOCK_NDOWN:?},${ONLY_BLOCK_MX:?},${ONLY_BLOCK_MY:?}"
  COMMAND+=(--only-block "$ONLY_BLOCK")
fi

printf '%q ' "${COMMAND[@]}" > "command_${FILE_TAG}.txt"
printf '\n' >> "command_${FILE_TAG}.txt"
cat > "provenance_job-${SLURM_JOB_ID:-unknown}_${FILE_TAG}.txt" <<EOF
git_commit=$DEPLOYED_COMMIT
source_dir=$SOURCE_DIR
ftlm_executable=$FTLM_EXEC
ftlm_sha256=$FTLM_SHA256
reducer_executable=$REDUCER_EXEC
reducer_sha256=$REDUCER_SHA256
slurm_job_id=${SLURM_JOB_ID:-}
node=$(hostname)
threads=$THREADS
target_R=$TARGET_R
target_m=$LANCZOS_MAX_STEPS
lanczos_save_steps=$LANCZOS_SAVE_STEPS
checkpoint_series=$CHECKPOINT_SERIES
requested_memory_gb=250
started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
echo "FTLM_JOB_START utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) host=$(hostname) job=${SLURM_JOB_ID:-} threads=$THREADS memory_gb=250 target_R=$TARGET_R target_m=$LANCZOS_MAX_STEPS save_m=$LANCZOS_SAVE_STEPS series=$CHECKPOINT_SERIES git_commit=$DEPLOYED_COMMIT ftlm_sha256=$FTLM_SHA256"

child_pid=''
request_stop() {
  echo "FTLM_WRAPPER event=STOP_REQUESTED signal=USR1 utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [[ -n $child_pid ]] && kill -0 "$child_pid" 2>/dev/null; then
    kill -USR1 "$child_pid"
  fi
}
request_term() {
  echo "FTLM_WRAPPER event=STOP_REQUESTED signal=TERM utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  if [[ -n $child_pid ]] && kill -0 "$child_pid" 2>/dev/null; then
    kill -TERM "$child_pid"
  fi
}
trap request_stop USR1
trap request_term TERM INT

set +e
"${COMMAND[@]}" &
child_pid=$!
while true; do
  wait "$child_pid"
  rc=$?
  if kill -0 "$child_pid" 2>/dev/null; then
    continue
  fi
  break
done
set -e

echo "FTLM_JOB_END utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) rc=$rc"
status_rc=1
if [[ -s $CHECKPOINT && -s $CHECKPOINT.meta ]]; then
  set +e
  "$REDUCER_EXEC" \
    --checkpoint "$CHECKPOINT" --samples "$TARGET_R" \
    --lanczos-steps "$LANCZOS_MAX_STEPS" \
    --status-json "status_${FILE_TAG}.json" \
    | tee "status_${FILE_TAG}.txt"
  status_rc=${PIPESTATUS[0]}
  set -e
else
  echo "FTLM_WRAPPER event=STATUS_UNAVAILABLE checkpoint_or_metadata_missing=1"
fi

if [[ $rc -eq 75 ]]; then
  if [[ $status_rc -ne 0 ]]; then
    echo "FTLM_WRAPPER event=FAILED reason=incomplete_checkpoint_unreadable" >&2
    exit 1
  fi
  echo "FTLM_WRAPPER event=CHECKPOINTED_INCOMPLETE target_R=$TARGET_R target_m=$LANCZOS_MAX_STEPS series=$CHECKPOINT_SERIES"
  sha256sum "$CHECKPOINT" "$CHECKPOINT.meta" > "SHA256_${FILE_TAG}.txt"
  exit 0
fi
if [[ $rc -ne 0 ]]; then
  echo "FTLM_WRAPPER event=FAILED executable_rc=$rc" >&2
  exit "$rc"
fi
if [[ $status_rc -ne 0 ]]; then
  echo "FTLM_WRAPPER event=FAILED reason=completed_run_checkpoint_unreadable" >&2
  exit 1
fi

if [[ -z ${ONLY_BLOCK_NUP:-} ]]; then
  REDUCE_R_VALUES=()
  for reduce_r in 8 16 32 64 128 256; do
    if (( reduce_r > TARGET_R )); then
      continue
    fi
    REDUCE_R_VALUES+=("$reduce_r")
  done
  if [[ ! " ${REDUCE_R_VALUES[*]} " =~ " ${TARGET_R} " ]]; then
    REDUCE_R_VALUES+=("$TARGET_R")
  fi
  if [[ $CHECKPOINT_SERIES == legacy ]]; then
    for reduce_r in "${REDUCE_R_VALUES[@]}"; do
      reduce_tag=$(printf '%03d' "$reduce_r")
      "$REDUCER_EXEC" \
        --checkpoint "$CHECKPOINT" --samples "$reduce_r" \
        --lanczos-steps "$LANCZOS_MAX_STEPS" \
        --beta-list 2.857142857142857,12.5 \
        --mu-min -3 --mu-max 4 --mu-count 281 \
        --output "$RUN_DIR/twist_${TWIST_ID}_thermo_R${reduce_tag}.csv"
    done
  else
    REDUCE_R_CSV=$(IFS=,; echo "${REDUCE_R_VALUES[*]}")
    "$REDUCER_EXEC" \
      --checkpoint "$CHECKPOINT" \
      --samples-list "$REDUCE_R_CSV" \
      --lanczos-steps-list "$LANCZOS_SAVE_STEPS" \
      --beta-list 2.857142857142857,12.5 \
      --mu-min -3 --mu-max 4 --mu-count 281 \
      --output-template "$RUN_DIR/twist_${TWIST_ID}_thermo_m{m}_R{R}.csv"
  fi
fi
sha256sum "$CHECKPOINT" "$CHECKPOINT.meta" "$RUN_DIR"/*.csv 2>/dev/null \
  > "SHA256_${FILE_TAG}.txt" || true
echo "FTLM_WRAPPER event=COMPLETE target_R=$TARGET_R target_m=$LANCZOS_MAX_STEPS series=$CHECKPOINT_SERIES"
