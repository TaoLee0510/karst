#!/usr/bin/env bash
# Fig. 2d truth study (scripts/fig2d_truth.R) on the HPC.
#   bash submit_f2t.sh setup                 cells, tasks, engine check (one job)
#   bash submit_f2t.sh fit SIMS [N]          step 2 fits (tasks_SIMS.tsv, N tasks)
#   bash submit_f2t.sh rows STEP SIMS N      rows of N simulations (DEPEND=<job> waits)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
# rc freeze: the study root is RC_RESULTS/part4_rate_choice/fig2d_truth; the round's
# fits and inputs come from RC_UPSTREAM (results_rc for a rerun of this study
# alone); the code is Code/paper (the as-run staged copy is ROOT/code in results_rc).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"
R6=${RC_UPSTREAM}/part1_paper/round_6ae2a24
ROOT=${RC_RESULTS}/part4_rate_choice/fig2d_truth
STAGE=${RC_CODE}/paper
EX=${RC_EXCLUDE}
monitor() { rc_monitor "$1"; }
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date --iso-8601=seconds)" >> ${RC_RESULTS}/logs/jobs.tsv; }
dep() { [[ -z "${DEPEND:-}" ]] || echo "--dependency=${DEPEND_TYPE:-afterok}:${DEPEND}"; }
mkdir -p ${ROOT}/logs ${RC_RESULTS}/logs
case "${1:?}" in
  setup)
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=16G --time=04:00:00 --exclude=${EX} --job-name=f2t_setup \
      --output=${ROOT}/logs/%x_%j.out --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${R6} \
      ${STAGE}/hpc/f2t_job.sh cells tasks engine)
    record ${j} f2t_setup "cells tasks engine"; echo "setup ${j} monitor $(monitor ${j})" ;;
  fit)
    sims=${2:?}; n=${3:-$(($(wc -l < ${ROOT}/tasks/tasks_${sims}.tsv) - 1))}
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${FIT_MEM:-8G} --time=12:00:00 --array=1-${n} $(dep) \
      --exclude=${EX} --job-name=f2t_fit_${sims} --output=${ROOT}/logs/%x_%A_%a.out \
      --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${ROOT}/step2,ABM_CELLS=${ROOT}/salehi_cells.tsv,ABM_FORECAST_HORIZON=3,TASK_FILE=${ROOT}/tasks/tasks_${sims}.tsv,BENCH_TASKS=none \
      ${STAGE}/hpc/f2t_job.sh fit)
    record ${j} f2t_fit_${sims} "1-${n}"; echo "fit ${j} monitor $(monitor ${j})" ;;
  rows)
    step=${2:?}; sims=${3:?}; n=${4:?}
    if [[ ${step} == 1 ]]; then res=${R6}; cells=""; else res=${ROOT}/step2; cells=${ROOT}/salehi_cells.tsv; fi
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${ROWS_MEM:-16G} --time=12:00:00 --array=1-${n} $(dep) \
      --exclude=${EX} --job-name=f2t_rows${step}_${sims} --output=${ROOT}/logs/%x_%A_%a.out \
      --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${res},ABM_CELLS=${cells},F2T_STEP=${step},F2T_SIMS=${sims} \
      ${STAGE}/hpc/f2t_job.sh rows)
    record ${j} f2t_rows${step}_${sims} "1-${n}"; echo "rows ${j} monitor $(monitor ${j})" ;;
esac
