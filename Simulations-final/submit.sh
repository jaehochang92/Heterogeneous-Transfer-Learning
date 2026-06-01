#!/bin/bash
#SBATCH --job-name=HTL_Sim
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --array=1-2000%1000
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard # update this

mkdir -p logs
mkdir -p results

module load R

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TASK_LOG="logs/task_${SLURM_ARRAY_JOB_ID}_${TASK_ID}.log"
exec > >(tee -a "$TASK_LOG") 2>&1

# 3. 환경 변수 및 디버깅 정보 출력 (첫 번째 Task에서만 출력하여 로그 남김)
if [ "$TASK_ID" == "1" ]; then
    echo "=========================================================="
    echo "Starting SLURM Array Job: $SLURM_ARRAY_JOB_ID"
  echo "Array Limit: 2000, Max Concurrent: 1000"
    echo "R Version: $(R --version | head -n 1)"
    echo "=========================================================="
fi

# 4. R 스크립트 실행
# 메인 스크립트에 파라미터를 전달합니다. rep 10으로 설정 시 총 9900번의 반복이 수행됩니다.
echo "[$(date '+%F %T')] Task $TASK_ID 시작..."

Rscript main.r \
  --seedno 1992 \
  --K 1 \
  --np 400 \
  --nt 30 \
  --p1 30 \
  --p2 30 \
  --dgp_fmap "nlinear" \
  --rep 10 \
  --K_values "1,2,4,8" \
  --np_values "200,400,800" \
  --nt_values "30,60,120" \
  --p1_values "30,50,100" \
  --p2_values "10,30,50" \
  --dgp_fmap_values "linear,nlinear"

exit_code=$?

if [ "$exit_code" -eq 0 ]; then
    echo "[$(date '+%F %T')] Task $TASK_ID 완료."
else
    echo "[$(date '+%F %T')] Task $TASK_ID 실패 (exit=$exit_code)."
fi

exit "$exit_code"