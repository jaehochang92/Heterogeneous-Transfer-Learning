#!/bin/bash -l
set -euo pipefail
module load gcc/12.3.0 R/4.5.2

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
P2_VALUES="10,20,30"

REPEATS_PER_SWEEP=200
MAX_ARRAY_SIZE=990

ACCOUNT="PAS1316"
PARTITION=""
TIME="24:00:00"
CPUS_PER_TASK=4
MEM_PER_CPU="4G"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"

SWEEP_PARAMS_CSV="K,np,nt,p1,p2"

parse_args() {
  local selected=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --worker)
        MODE="worker"
        shift
        ;;
      --K)
        selected+=("K")
        shift
        ;;
      --np)
        selected+=("np")
        shift
        ;;
      --nt)
        selected+=("nt")
        shift
        ;;
      --p1)
        selected+=("p1")
        shift
        ;;
      --p2)
        selected+=("p2")
        shift
        ;;
      --all)
        selected=("K" "np" "nt" "p1" "p2")
        shift
        ;;
      *)
        echo "Unknown argument: $1"
        echo "Usage: bash submit.sh [--K] [--np] [--nt] [--p1] [--p2] [--all]"
        exit 1
        ;;
    esac
  done

  if [[ ${#selected[@]} -eq 0 ]]; then
    selected=("K" "np" "nt" "p1" "p2")
  fi

  local uniq=()
  local item
  for item in "${selected[@]}"; do
    if [[ " ${uniq[*]} " != *" ${item} "* ]]; then
      uniq+=("$item")
    fi
  done
  SWEEP_PARAMS_CSV="$(IFS=,; echo "${uniq[*]}")"
}

count_sweeps() {
  Rscript --vanilla -e '
parse_int_list <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals <- vals[nzchar(vals)]
  as.integer(vals)
}
selected <- trimws(unlist(strsplit(Sys.getenv("SWEEP_PARAMS_CSV"), ",", fixed = TRUE)))
selected <- selected[nzchar(selected)]
if (length(selected) == 0) {
  stop("No sweep parameter selected")
}

fixed_cfg <- list(
  K = as.integer(Sys.getenv("FIXED_K")),
  np = as.integer(Sys.getenv("FIXED_NP")),
  nt = as.integer(Sys.getenv("FIXED_NT")),
  p1 = as.integer(Sys.getenv("FIXED_P1")),
  p2 = as.integer(Sys.getenv("FIXED_P2"))
)
sweep_values <- list(
  K = parse_int_list(Sys.getenv("K_VALUES")),
  np = parse_int_list(Sys.getenv("NP_VALUES")),
  nt = parse_int_list(Sys.getenv("NT_VALUES")),
  p1 = parse_int_list(Sys.getenv("P1_VALUES")),
  p2 = parse_int_list(Sys.getenv("P2_VALUES"))
)
sweep_order <- c("K", "np", "nt", "p1", "p2")
rows <- list()
idx <- 1L
for (param in sweep_order) {
  if (!(param %in% selected)) next
  vals <- unique(sweep_values[[param]])
  for (v in vals) {
    cfg <- fixed_cfg
    cfg[[param]] <- v
    rows[[idx]] <- data.frame(
      sweep_param = param,
      K = as.integer(cfg$K),
      np = as.integer(cfg$np),
      nt = as.integer(cfg$nt),
      p1 = as.integer(cfg$p1),
      p2 = as.integer(cfg$p2)
    )
    idx <- idx + 1L
  }
}
plan <- unique(do.call(rbind, rows))
cat(nrow(plan))
' | tr -d '[:space:]'
}

compute_totals() {
  export FIXED_K FIXED_NP FIXED_NT FIXED_P1 FIXED_P2
  export K_VALUES NP_VALUES NT_VALUES P1_VALUES P2_VALUES SWEEP_PARAMS_CSV
  TOTAL_SWEEPS="$(count_sweeps)"
  TOTAL_JOBS=$(( TOTAL_SWEEPS * REPEATS_PER_SWEEP ))
}

run_worker() {
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is required in worker mode."
    exit 1
  fi

  local local_task_id="$SLURM_ARRAY_TASK_ID"
  local submitted_tasks="${SUBMITTED_TASKS:-1}"

  if [[ -z "${TOTAL_SWEEPS:-}" || -z "${TOTAL_JOBS:-}" ]]; then
    compute_totals
  fi

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
      --p2_values "$P2_VALUES" \
      --sweep_params "$SWEEP_PARAMS_CSV"
  done
}

submit_single_array() {
  local submitted_tasks="$1"
  echo "Submitting simulation jobs as a single array (max ${MAX_ARRAY_SIZE} tasks)."
  echo "Selected sweeps: ${SWEEP_PARAMS_CSV}"
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
    --export="ALL,SUBMITTED_TASKS=${submitted_tasks},TOTAL_SWEEPS=${TOTAL_SWEEPS},TOTAL_JOBS=${TOTAL_JOBS},SEED_BASE=${SEED_BASE},FIXED_K=${FIXED_K},FIXED_NP=${FIXED_NP},FIXED_NT=${FIXED_NT},FIXED_P1=${FIXED_P1},FIXED_P2=${FIXED_P2},K_VALUES=${K_VALUES},NP_VALUES=${NP_VALUES},NT_VALUES=${NT_VALUES},P1_VALUES=${P1_VALUES},P2_VALUES=${P2_VALUES},SWEEP_PARAMS_CSV=${SWEEP_PARAMS_CSV}"
  )

  if [[ -n "$PARTITION" ]]; then
    sbatch_args+=(--partition="$PARTITION")
  fi

  sbatch "${sbatch_args[@]}" "$0" --worker
}

MODE="submit"
parse_args "$@"

if [[ "$MODE" == "submit" ]]; then
  compute_totals
fi

if [[ "$MODE" == "worker" ]]; then
  run_worker
else
  SUBMITTED_TASKS=$(( TOTAL_JOBS < MAX_ARRAY_SIZE ? TOTAL_JOBS : MAX_ARRAY_SIZE ))
  submit_single_array "$SUBMITTED_TASKS"
fi