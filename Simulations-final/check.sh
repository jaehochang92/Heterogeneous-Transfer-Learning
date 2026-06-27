#!/bin/bash
# =============================================================================
# SAND Simulation Real-time Progress Monitor (Dynamic Config Binding)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_SCRIPT="$SCRIPT_DIR/submit_sand.sh"

if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "[-] ERROR: $SUBMIT_SCRIPT 파일을 찾을 수 없습니다."
    exit 1
fi

RESULTS_DIR="$SCRIPT_DIR/results_sand/batch_bundles"
START_TIME=$(date +%s)

while true; do
    if [ ! -d "$RESULTS_DIR" ]; then
        CURRENT_FILES=0
        ACTUAL_ROWS=0
    else
        CURRENT_FILES=$(ls -1 "$RESULTS_DIR"/sand_res_*.rds 2>/dev/null | wc -l)
        ACTUAL_ROWS=$(Rscript --vanilla -e '
            files <- list.files("results_sand/batch_bundles", pattern = "^sand_res_.*[.]rds$", full.names = TRUE)
            if(length(files) == 0) { cat(0) } else { cat(sum(sapply(files, function(f) nrow(readRDS(f))))) }
        ' 2>/dev/null || echo 0)
    fi

    CURRENT_TIME=$(date +%s)
    ELAPSED=$(( CURRENT_TIME - START_TIME ))
    format_time() {
        local t=$1
        echo "$((t/3600))h $(( (t%3600)/60 ))m $((t%60))s"
    }
    tput cup 7 0
    echo -e " 현재 시각 : $(date '+%Y-%m-%d %H:%M:%S')"
    echo -e " 경과 시간 : $(format_time $ELAPSED)"
    echo -e "-------------------------------------------------------------------------"
    echo -e " 수집 현황 : $(printf "%'d" $ACTUAL_ROWS) Rows 완료"
    echo -e " 파일 적재 : 총 ${CURRENT_FILES} 개의 배치 번들 rds 파일 확보 완료"
    echo -e "-------------------------------------------------------------------------"
    RUNNING_JOBS=$(squeue -u $USER -h | wc -l)
    echo -e " 슬럼 상태 : 현재 스케줄러에서 ${RUNNING_JOBS} 개의 작업 노드가 구동 중입니다."
    sleep 5
done
EOF

chmod +x monitor_progress.sh