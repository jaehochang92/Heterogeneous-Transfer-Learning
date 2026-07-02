#!/bin/bash -l
# Submit parallel SLURM array jobs for CuratedOvarianCancer/pilot/2.main.R
#
# Usage:
#   bash sbmt.sh                 # submit with defaults
#   bash sbmt.sh 20              # submit 20 array tasks
#   N_REPLICATES=100 bash sbmt.sh 10
#   bash sbmt.sh --worker        # internal worker entry (called by sbatch)

set -euo pipefail
module load gcc/12.3.0 R/4.5.2

# --- schedule ---------------------------------------------------------------
ACCOUNT="PAS1316"
PARTITION=""
TIME="24:00:00"
CPUS_PER_TASK=8
MEM_PER_CPU="4G"

# total CV replicates split across array tasks
N_REPLICATES="${N_REPLICATES:-50}"
N_ARRAY_TASKS="${1:-10}"
SEED_BASE="${SEED_BASE:-2175}"
N_FOLDS="${N_FOLDS:-5}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="$SCRIPT_DIR/logs"
CACHE_DIR="$SCRIPT_DIR/cache"
mkdir -p "$LOG_DIR" "$CACHE_DIR"

MODE="submit"
if [[ "${1:-}" == "--worker" ]]; then
  MODE="worker"
elif [[ "${2:-}" == "--worker" ]]; then
  MODE="worker"
fi

run_worker() {
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is required in worker mode."
    exit 1
  fi

  local bundle_id="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-job}}_${SLURM_ARRAY_TASK_ID}"
  export SIM_BUNDLE_ID="${bundle_id}"
  export SUBMITTED_TASKS="${SUBMITTED_TASKS:-1}"
  export N_REPLICATES="${N_REPLICATES}"
  export SEED_BASE="${SEED_BASE}"
  export N_FOLDS="${N_FOLDS}"

  cd "$SCRIPT_DIR"
  echo "Worker task ${SLURM_ARRAY_TASK_ID}/${SUBMITTED_TASKS} | bundle=${bundle_id}"
  Rscript 2.main.R
}

submit_array() {
  local n_tasks="$1"

  if [[ "$n_tasks" -lt 1 ]]; then
    echo "N_ARRAY_TASKS must be >= 1"
    exit 1
  fi
  if [[ "$N_REPLICATES" -lt "$n_tasks" ]]; then
    echo "Warning: N_REPLICATES (${N_REPLICATES}) < N_ARRAY_TASKS (${n_tasks}); some tasks will be idle."
  fi

  local -a sbatch_args=(
    --account="$ACCOUNT"
    --job-name="ovarian_pilot"
    --time="$TIME"
    --cpus-per-task="$CPUS_PER_TASK"
    --mem-per-cpu="$MEM_PER_CPU"
    --output="$LOG_DIR/%A_%a.out"
    --error="$LOG_DIR/%A_%a.err"
    --array="1-${n_tasks}"
    --export="ALL,SUBMITTED_TASKS=${n_tasks},N_REPLICATES=${N_REPLICATES},SEED_BASE=${SEED_BASE},N_FOLDS=${N_FOLDS}"
  )

  if [[ -n "$PARTITION" ]]; then
    sbatch_args+=(--partition="$PARTITION")
  fi

  echo "Submitting ${n_tasks} array tasks"
  echo "  N_REPLICATES=${N_REPLICATES}"
  echo "  SEED_BASE=${SEED_BASE}"
  echo "  N_FOLDS=${N_FOLDS}"
  echo "  LOG_DIR=${LOG_DIR}"
  echo "  CACHE_DIR=${CACHE_DIR}"

  sbatch "${sbatch_args[@]}" "$0" --worker
}

if [[ "$MODE" == "worker" ]]; then
  run_worker
else
  submit_array "$N_ARRAY_TASKS"
fi
