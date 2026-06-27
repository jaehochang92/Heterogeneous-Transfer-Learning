#!/bin/bash -l
set -euo pipefail
module load gcc/12.3.0 R/4.5.2

# Schedule
REPEATS_PER_SWEEP=300
MAX_ARRAY_SIZE=990
ACCOUNT="PAS1316"
PARTITION=""
TIME="02:00:00"
CPUS_PER_TASK=4
MEM_PER_CPU="6G"

# Simulation design
SEED_BASE=1992
GAMMA_SCALE=1
FIXED_H=0.6
FIXED_K=4
FIXED_NP=400
FIXED_NT=30
FIXED_P1=30
FIXED_P2=30
H_VALUES="${H_VALUES:-4,2,1,0.5,0.25,0.12}"
K_VALUES="${K_VALUES:-4,8,16,32}"
NP_VALUES="${NP_VALUES:-200,400,800,1600}"
NT_VALUES="${NT_VALUES:-30,60,120,240}"
P1_VALUES="${P1_VALUES:-20,40,80,160}"
P2_VALUES="${P2_VALUES:-10,20,40,80}"

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
LOG_DIR="$SCRIPT_DIR/logs_sand"
mkdir -p "$LOG_DIR"

SWEEP_PARAMS_CSV="${SWEEP_PARAMS_CSV:-K,np,nt,p1,p2}"

parse_args() {
  local selected=()
  MODE="submit"
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --worker)
        MODE="worker"
        shift
        ;;
      --h)
        selected+=("h")
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
        selected=("h" "K" "np" "nt" "p1" "p2")
        shift
        ;;
      *)
        echo "Unknown argument: $1"
        echo "Usage: bash submit_sand.sh [--h] [--K] [--np] [--nt] [--p1] [--p2] [--all]"
        exit 1
        ;;
    esac
  done

  if [[ "$MODE" == "worker" ]]; then
    SWEEP_PARAMS_CSV="${SWEEP_PARAMS_CSV:-h,K,np,nt,p1,p2}"
  else
    if [[ ${#selected[@]} -eq 0 ]]; then
      selected=("h" "K" "np" "nt" "p1" "p2")
    fi

    local uniq=()
    local item
    for item in "${selected[@]}"; do
      if [[ " ${uniq[*]} " != *" ${item} "* ]]; then
        uniq+=("$item")
      fi
    done
    SWEEP_PARAMS_CSV="$(IFS=,; echo "${uniq[*]}")"
  fi
}

count_sweeps() {
  Rscript --vanilla -e '
    parse_int_list <- function(x) {
      vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
      vals <- vals[nzchar(vals)]
      as.integer(vals)
    }
    parse_dbl_list <- function(x) {
      vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
      vals <- vals[nzchar(vals)]
      as.double(vals)
    }
    selected <- trimws(unlist(strsplit(Sys.getenv("SWEEP_PARAMS_CSV"), ",", fixed = TRUE)))
    selected <- selected[nzchar(selected)]
    if (length(selected) == 0) {
      stop("No sweep parameter selected")
    }
    
    fixed_cfg <- list(
      h = as.double(Sys.getenv("FIXED_H")),
      K = as.integer(Sys.getenv("FIXED_K")),
      np = as.integer(Sys.getenv("FIXED_NP")),
      nt = as.integer(Sys.getenv("FIXED_NT")),
      p1 = as.integer(Sys.getenv("FIXED_P1")),
      p2 = as.integer(Sys.getenv("FIXED_P2"))
    )
    sweep_values <- list(
      h = parse_dbl_list(Sys.getenv("H_VALUES")),
      K = parse_int_list(Sys.getenv("K_VALUES")),
      np = parse_int_list(Sys.getenv("NP_VALUES")),
      nt = parse_int_list(Sys.getenv("NT_VALUES")),
      p1 = parse_int_list(Sys.getenv("P1_VALUES")),
      p2 = parse_int_list(Sys.getenv("P2_VALUES"))
    )
    sweep_order <- c("h", "K", "np", "nt", "p1", "p2")
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
          h = as.double(cfg$h),
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
  export FIXED_H FIXED_K GAMMA_SCALE FIXED_NP FIXED_NT FIXED_P1 FIXED_P2
  export H_VALUES K_VALUES NP_VALUES NT_VALUES P1_VALUES P2_VALUES SWEEP_PARAMS_CSV
  TOTAL_SWEEPS="$(count_sweeps)"
  TOTAL_JOBS=$(( TOTAL_SWEEPS * REPEATS_PER_SWEEP ))
}

clean_batch_directories() {
  local bundle_dir="$SCRIPT_DIR/results_sand/batch_bundles"
  if [[ -d "$bundle_dir" ]]; then
    echo "Cleaning previous batch bundles in: $bundle_dir"
    rm -rf "$bundle_dir"
  fi
}

run_worker() {
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is required in worker mode."
    exit 1
  fi
  local array_bundle_base="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-job}}_${SLURM_ARRAY_TASK_ID}"
  
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

  export SIM_BUNDLE_ID="${array_bundle_base}"

  Rscript main_sand.R \
    --gamma_scale "$GAMMA_SCALE" \
    --h "$FIXED_H" \
    --K "$FIXED_K" \
    --np "$FIXED_NP" \
    --nt "$FIXED_NT" \
    --p1 "$FIXED_P1" \
    --p2 "$FIXED_P2" \
    --h_values "$H_VALUES" \
    --K_values "$K_VALUES" \
    --np_values "$NP_VALUES" \
    --nt_values "$NT_VALUES" \
    --p1_values "$P1_VALUES" \
    --p2_values "$P2_VALUES" \
    --sweep_params "$SWEEP_PARAMS_CSV" \
    --out_dir "results_sand" \
    --log_dir "rlogs_sand" \
    --task_start "$start_idx" \
    --task_count "$my_load"
}

submit_single_array() {
  local submitted_tasks="$1"
  echo "Submitting simulation jobs as a single array (max ${MAX_ARRAY_SIZE} tasks)."
  echo "Selected sweeps: ${SWEEP_PARAMS_CSV}"
  echo "Sweep count: ${TOTAL_SWEEPS}"
  echo "Repeats per sweep: ${REPEATS_PER_SWEEP}"
  echo "Total jobs: ${TOTAL_JOBS}"
  echo "Submitted array tasks: ${submitted_tasks}"

  local output_file="$LOG_DIR/%x_%A_%a.log"

  local current_h_vals="$FIXED_H"
  local current_k_vals="$FIXED_K"
  local current_np_vals="$FIXED_NP"
  local current_nt_vals="$FIXED_NT"
  local current_p1_vals="$FIXED_P1"
  local current_p2_vals="$FIXED_P2"
  
  if [[ ",${SWEEP_PARAMS_CSV}," == *,h,* ]];  then current_h_vals="$H_VALUES"; fi
  if [[ ",${SWEEP_PARAMS_CSV}," == *,K,* ]];  then current_k_vals="$K_VALUES"; fi
  if [[ ",${SWEEP_PARAMS_CSV}," == *,np,* ]]; then current_np_vals="$NP_VALUES"; fi
  if [[ ",${SWEEP_PARAMS_CSV}," == *,nt,* ]]; then current_nt_vals="$NT_VALUES"; fi
  if [[ ",${SWEEP_PARAMS_CSV}," == *,p1,* ]]; then current_p1_vals="$P1_VALUES"; fi
  if [[ ",${SWEEP_PARAMS_CSV}," == *,p2,* ]]; then current_p2_vals="$P2_VALUES"; fi
  
  export SUBMITTED_TASKS="${submitted_tasks}"
  export TOTAL_SWEEPS="${TOTAL_SWEEPS}"
  export TOTAL_JOBS="${TOTAL_JOBS}"
  export SEED_BASE="${SEED_BASE}"
  export GAMMA_SCALE="${GAMMA_SCALE}"
  export FIXED_H="${FIXED_H}"
  export FIXED_K="${FIXED_K}"
  export FIXED_NP="${FIXED_NP}"
  export FIXED_NT="${FIXED_NT}"
  export FIXED_P1="${FIXED_P1}"
  export FIXED_P2="${FIXED_P2}"
  export H_VALUES="${current_h_vals}"
  export K_VALUES="${current_k_vals}"
  export NP_VALUES="${current_np_vals}"
  export NT_VALUES="${current_nt_vals}"
  export P1_VALUES="${current_p1_vals}"
  export P2_VALUES="${current_p2_vals}"
  export SWEEP_PARAMS_CSV="${SWEEP_PARAMS_CSV}"
  
  local -a sbatch_args=(
    --job-name="sand_sim"
    --output="$output_file"
    --time="$TIME"
    --ntasks=1
    --cpus-per-task="$CPUS_PER_TASK"
    --mem-per-cpu="$MEM_PER_CPU"
    --account="$ACCOUNT"
    --array="1-${submitted_tasks}"
    --export=ALL
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
  clean_batch_directories
fi

if [[ "$MODE" == "worker" ]]; then
  run_worker
else
  SUBMITTED_TASKS=$(( TOTAL_JOBS < MAX_ARRAY_SIZE ? TOTAL_JOBS : MAX_ARRAY_SIZE ))
  submit_single_array "$SUBMITTED_TASKS"
fi