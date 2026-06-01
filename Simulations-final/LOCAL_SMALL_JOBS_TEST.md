# Local Small Jobs Submit Test Guide

Here we list some example for submitting jobs in SLURM.

## 1) 단일 태스크 빠른 스모크 테스트

한 번만 돌려서 파이프라인이 깨지지 않는지 확인합니다.

    cd ~/R/Heterogeneous-Transfer-Learning
    SLURM_ARRAY_TASK_ID=1 Rscript main.r \
      --seedno 100 \
      --rep 1 \
      --K 1 \
      --np 80 \
      --nt 20 \
      --p1 10 \
      --p2 5 \
      --dgp_fmap nlinear \
      --K_values 1 \
      --np_values 80 \
      --nt_values 20 \
      --p1_values 10 \
      --p2_values 5 \
      --dgp_fmap_values nlinear

확인 포인트:
- logs/progress_task_0001.log 생성
- results/sim_res_task_0001.rds 생성
- 결과에 Sp 컬럼으로 sp, nsp 두 케이스가 모두 포함

## 2) 로컬에서 배열 태스크 흉내내기 (순차)

작은 배열을 순차 실행해 sweep 분기와 파일 저장을 점검합니다.

    cd ~/R/Heterogeneous-Transfer-Learning
    for tid in 1 2 3 4; do
      echo "Running local task ${tid}"
      SLURM_ARRAY_TASK_ID=${tid} Rscript main.r \
        --seedno 100 \
        --rep 1 \
        --K 1 \
        --np 80 \
        --nt 20 \
        --p1 10 \
        --p2 5 \
        --dgp_fmap nlinear \
        --K_values 1,2 \
        --np_values 80 \
        --nt_values 20 \
        --p1_values 10 \
        --p2_values 5 \
        --dgp_fmap_values linear,nlinear
    done

확인 포인트:
- results/sim_res_task_0001.rds, 0002.rds ... 생성
- task id가 sweep 계획을 올바르게 선택하는지 확인

## 3) 로컬 병렬 소규모 테스트

CPU를 과도하게 쓰지 않도록 2개 동시 실행만 허용합니다.

권장 (xargs 없이 안정적으로 병렬 제어):

        cd ~/R/Heterogeneous-Transfer-Learning
        max_jobs=2
        for tid in 1 2 3 4; do
            SLURM_ARRAY_TASK_ID="$tid" Rscript main.r \
                --seedno 100 \
                --rep 1 \
                --K 1 \
                --np 80 \
                --nt 20 \
                --p1 10 \
                --p2 5 \
                --dgp_fmap nlinear \
                --K_values 1,2 \
                --np_values 80 \
                --nt_values 20 \
                --p1_values 10 \
                --p2_values 5 &
            while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$max_jobs" ]; do
                wait -n 2>/dev/null || wait
            done
        done
        wait

이건 간소화 버전:

        cd ~/R/Heterogeneous-Transfer-Learning
        max_jobs=2
        for tid in 1 2 3 4; do
            SLURM_ARRAY_TASK_ID="$tid" Rscript main.r &
            while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$max_jobs" ]; do
                wait -n 2>/dev/null || wait
            done
        done
        wait

xargs를 꼭 쓰고 싶다면 한 줄 버전:

        cd ~/R/Heterogeneous-Transfer-Learning
        printf "%s\n" 1 2 3 4 | xargs -P2 -I{} bash -lc 'SLURM_ARRAY_TASK_ID={} Rscript main.r --seedno 100 --rep 1 --K 1 --np 80 --nt 20 --p1 10 --p2 5 --dgp_fmap nlinear --K_values 1,2 --np_values 80 --nt_values 20 --p1_values 10 --p2_values 5 --dgp_fmap_values linear,nlinear'

확인 포인트:
- 로그 파일이 서로 덮어쓰지 않는지 확인
- 결과 파일이 태스크 번호별로 분리되는지 확인

## 4) submit.sh를 로컬에서 직접 실행해보기

sbatch 없이 쉘 스크립트 동작만 확인하고 싶을 때 사용합니다.

    cd ~/R/Heterogeneous-Transfer-Learning
    export SLURM_ARRAY_JOB_ID=9999
    export SLURM_ARRAY_TASK_ID=1
    bash submit.sh

참고:
- module load R 명령이 로컬 환경에 없으면 실패할 수 있습니다.
- 그 경우에는 위 1), 2), 3) 방식처럼 Rscript main.r을 직접 호출하세요.

## 5) 결과 종합 시각화 실행

small jobs 테스트 후 그래프 생성 확인용입니다.

    cd ~/R/Heterogeneous-Transfer-Learning
    Rscript plot.R results results/plots

생성 파일:
- results/plots/summary_est.png
- results/plots/summary_pe.png
