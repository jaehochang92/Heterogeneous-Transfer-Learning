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

ACCOUNT="PAS1316"
PARTITION=""
TIME="24:00:00"
CPUS_PER_TASK=1
MEM_PER_CPU="4G"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"

count_sweeps() {
  Rscript --vanilla -e '
parse_int_list <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals <- vals[nzchar(vals)]
  as.integer(vals)
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
  vals <- unique(sweep_values[[param]])
  for (v in vals) {
    cfg <- fixed_cfg
    cfg[[param]] <- v
    rows[[idx]] <- data.frame(
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

run_worker() {
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is required in worker mode."
    exit 1
  fi

  local global_task_id="$SLURM_ARRAY_TASK_ID"
  local sweep_task_id=$(( ((global_task_id - 1) % TOTAL_SWEEPS) + 1 ))
  local repeat_idx=$(( ((global_task_id - 1) / TOTAL_SWEEPS) + 1 ))
  local seed=$(( SEED_BASE + repeat_idx - 1 ))

  export SIM_GLOBAL_TASK_ID="$global_task_id"
  export SIM_SWEEP_TASK_ID="$sweep_task_id"

  cd "$SCRIPT_DIR"
  ml load gcc/12.3.0 R/4.4.0

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
}

submit_chunks() {
  echo "Submitting simulation jobs in chunks (max ${MAX_ARRAY_SIZE} per array)."
  echo "Sweep count: ${TOTAL_SWEEPS}"
  echo "Repeats per sweep: ${REPEATS_PER_SWEEP}"
  echo "Total jobs: ${TOTAL_JOBS}"

  local start=1
  local chunk=1
  while [[ $start -le $TOTAL_JOBS ]]; do
    local end=$(( start + MAX_ARRAY_SIZE - 1 ))
    if [[ $end -gt $TOTAL_JOBS ]]; then
      end=$TOTAL_JOBS
    fi

    local job_name="htl_sim_${chunk}"
    local output_file="$LOG_DIR/%x_%A_%a.out"
    local error_file="$LOG_DIR/%x_%A_%a.err"

    local -a sbatch_args=(
      --job-name="$job_name"
      --output="$output_file"
      --error="$error_file"
      --time="$TIME"
      --ntasks=1
      --cpus-per-task="$CPUS_PER_TASK"
      --mem-per-cpu="$MEM_PER_CPU"
      --account="$ACCOUNT"
      --array="${start}-${end}"
      --export="ALL,TOTAL_SWEEPS=${TOTAL_SWEEPS},SEED_BASE=${SEED_BASE},FIXED_K=${FIXED_K},FIXED_NP=${FIXED_NP},FIXED_NT=${FIXED_NT},FIXED_P1=${FIXED_P1},FIXED_P2=${FIXED_P2},K_VALUES=${K_VALUES},NP_VALUES=${NP_VALUES},NT_VALUES=${NT_VALUES},P1_VALUES=${P1_VALUES},P2_VALUES=${P2_VALUES}"
    )

    if [[ -n "$PARTITION" ]]; then
      sbatch_args+=(--partition="$PARTITION")
    fi

    sbatch "${sbatch_args[@]}" "$0" --worker
    echo "Submitted chunk ${chunk}: array ${start}-${end}"

    start=$(( end + 1 ))
    chunk=$(( chunk + 1 ))
  done
}

TOTAL_SWEEPS="${TOTAL_SWEEPS:-}"
if [[ -z "$TOTAL_SWEEPS" ]]; then
  export FIXED_K FIXED_NP FIXED_NT FIXED_P1 FIXED_P2
  export K_VALUES NP_VALUES NT_VALUES P1_VALUES P2_VALUES
  TOTAL_SWEEPS="$(count_sweeps)"
fi
TOTAL_JOBS=$(( TOTAL_SWEEPS * REPEATS_PER_SWEEP ))

if [[ "${1:-}" == "--worker" ]]; then
  run_worker
else
  submit_chunks
fi
