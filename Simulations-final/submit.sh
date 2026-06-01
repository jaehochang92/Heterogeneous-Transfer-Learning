#!/bin/bash
#SBATCH --job-name=HTL_Sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-2000%1000
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard # update this

mkdir -p logs
mkdir -p results

module load R

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TASK_LOG="logs/task_${SLURM_ARRAY_JOB_ID}_${TASK_ID}.log"
exec > >(tee -a "$TASK_LOG") 2>&1

# 3. Env vars and debugging info (from the 1st task's log)
if [ "$TASK_ID" == "1" ]; then
    echo "=========================================================="
    echo "Starting SLURM Array Job: $SLURM_ARRAY_JOB_ID"
  echo "Array Limit: 2000, Max Concurrent: 1000"
    echo "R Version: $(R --version | head -n 1)"
    echo "=========================================================="
fi

# 4. Run R script
# Pass params to the main script. rep=10 makes 9,900 reps in total.
echo "[$(date '+%F %T')] Task $TASK_ID 시작..."

Rscript main.r \
  --seedno 1992 \
  --K 1 \
  --np 400 \
  --nt 30 \
  --p1 30 \
  --p2 30 \
  --rep 10 \
  --K_values "1,2,4,8" \
  --np_values "200,400,800" \
  --nt_values "30,60,120" \
  --p1_values "30,50,100" \
  --p2_values "10,30,50"

exit_code=$?

if [ "$exit_code" -eq 0 ]; then
    echo "[$(date '+%F %T')] Task $TASK_ID succeeded."
else
    echo "[$(date '+%F %T')] Task $TASK_ID failed (exit=$exit_code)."
fi

exit "$exit_code"