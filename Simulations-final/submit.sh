#!/bin/bash -l
set -euo pipefail

# Simulation design
SEED_BASE=1992
FIXED_K=1
FIXED_NP=400
FIXED_NT=30
FIXED_P1=30
FIXED_P2=30
K_VALUES="2,4,8,16"
NP_VALUES="200,400,800,1600"
NT_VALUES="30,60,120,240"
P1_VALUES="30,50,100,200"
P2_VALUES="30,50,100,200"

REPEATS_PER_SWEEP=200
MAX_ARRAY_SIZE=990
TOTAL_SWEEPS=20
TOTAL_JOBS=4000

ACCOUNT="PAS1316"
PARTITION=""
TIME="24:00:00"
CPUS_PER_TASK=1
MEM_PER_CPU="4G"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"

run_worker() {
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is required in worker mode."
    exit 1
  fi

  local local_task_id="$SLURM_ARRAY_TASK_ID"
  local submitted_tasks="${SUBMITTED_TASKS:-1}"

  local base_load=$(( TOTAL_JOBS / submitted_tasks ))
  local extra=$(( TOTAL_JOBS % submitted_tasks ))

  local my_load
  local start_idx
  if [[ $local_task_id -le $extra ]]; then
    my_load=$(( base_load + 1 ))
    start_idx=$(( (local_task_id - 1) * (base_load + 1) + 1 ))
  else
    my_load=$base_load
    start_idx=$(( extra * (base_load + 1) + (local_task_id - extra - 1) * base_load + 1 ))
  fi

  if [[ $my_load -le 0 ]]; then
    echo "Task ${local_task_id}: no assigned workload, exiting."
    exit 0
  fi

  cd "$SCRIPT_DIR"
  ml load gcc/12.3.0 R/4.4.0

  local offset
  for (( offset = 0; offset < my_load; offset++ )); do
    local global_task_id=$(( start_idx + offset ))
    local sweep_task_id=$(( ((global_task_id - 1) % TOTAL_SWEEPS) + 1 ))
    local repeat_idx=$(( ((global_task_id - 1) / TOTAL_SWEEPS) + 1 ))
    local seed=$(( SEED_BASE + repeat_idx - 1 ))

    export SIM_GLOBAL_TASK_ID="$global_task_id"

    Rscript main.r \
      --seedno "$seed" \
      --K "$FIXED_K" \
      --np "$FIXED_NP" \
      --nt "$FIXED_NT" \
      --p1 "$FIXED_P1" \
      --p2 "$FIXED_P2" \
      --rep 1 \
      --K_values "$K_VALUES" \
      --np_values "$NP_VALUES" \
      --nt_values "$NT_VALUES" \
      --p1_values "$P1_VALUES" \
      --p2_values "$P2_VALUES"
  done
}

submit_single_array() {
  local submitted_tasks="$1"
  echo "Submitting simulation jobs as a single array (max ${MAX_ARRAY_SIZE} tasks)."
  echo "Sweep count: ${TOTAL_SWEEPS}"
  echo "Repeats per sweep: ${REPEATS_PER_SWEEP}"
  echo "Total jobs: ${TOTAL_JOBS}"
  echo "Submitted array tasks: ${submitted_tasks}"

  local output_file="$LOG_DIR/%x_%A_%a.out"
  local error_file="$LOG_DIR/%x_%A_%a.err"

  local -a sbatch_args=(
    --job-name="htl_sim"
    --output="$output_file"
    --error="$error_file"
    --time="$TIME"
    --ntasks=1
    --cpus-per-task="$CPUS_PER_TASK"
    --mem-per-cpu="$MEM_PER_CPU"
    --account="$ACCOUNT"
    --array="1-${submitted_tasks}"
    --export="ALL,SUBMITTED_TASKS=${submitted_tasks},TOTAL_SWEEPS=${TOTAL_SWEEPS},TOTAL_JOBS=${TOTAL_JOBS},SEED_BASE=${SEED_BASE},FIXED_K=${FIXED_K},FIXED_NP=${FIXED_NP},FIXED_NT=${FIXED_NT},FIXED_P1=${FIXED_P1},FIXED_P2=${FIXED_P2},K_VALUES=${K_VALUES},NP_VALUES=${NP_VALUES},NT_VALUES=${NT_VALUES},P1_VALUES=${P1_VALUES},P2_VALUES=${P2_VALUES}"
  )

  if [[ -n "$PARTITION" ]]; then
    sbatch_args+=(--partition="$PARTITION")
  fi

  sbatch "${sbatch_args[@]}" "$0" --worker
}

TOTAL_SWEEPS="${TOTAL_SWEEPS:-20}"
TOTAL_JOBS="${TOTAL_JOBS:-4000}"

if [[ "${1:-}" == "--worker" ]]; then
  run_worker
else
  SUBMITTED_TASKS=$(( TOTAL_JOBS < MAX_ARRAY_SIZE ? TOTAL_JOBS : MAX_ARRAY_SIZE ))
  submit_single_array "$SUBMITTED_TASKS"
fi
