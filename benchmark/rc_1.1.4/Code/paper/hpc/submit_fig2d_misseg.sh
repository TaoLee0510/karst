#!/usr/bin/env bash
# Fig. 2d on ABM_varyMisseg (scripts/fig2d_truth.R --source=misseg) on the HPC.
#   bash submit_f2m.sh setup                 cells, tasks, dispersion (one job)
#   bash submit_f2m.sh fit ARM SIMS [N]      fits + forecasts (tasks_ARM_SIMS.tsv)
#   bash submit_f2m.sh rows ARM SIMS N       rows of N simulations (DEPEND=<job> waits)
#   ARRAY=<indices> runs only those array tasks (reruns)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
# rc freeze: the study root is RC_RESULTS/part4_rate_choice/fig2d_misseg; the round's
# fits and inputs come from RC_UPSTREAM (results_rc for a rerun of this study
# alone); the code is Code/paper (the as-run staged copy is ROOT/code in results_rc).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"
R6=${RC_UPSTREAM}/part1_paper/round_6ae2a24
ROOT=${RC_RESULTS}/part4_rate_choice/fig2d_misseg
STAGE=${RC_CODE}/paper
EX=${RC_EXCLUDE}
monitor() { rc_monitor "$1"; }
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date --iso-8601=seconds)" >> ${RC_RESULTS}/logs/jobs.tsv; }
dep() { [[ -z "${DEPEND:-}" ]] || echo "--dependency=${DEPEND_TYPE:-afterok}:${DEPEND}"; }
mkdir -p ${ROOT}/logs ${RC_RESULTS}/logs
case "${1:?}" in
  setup)
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=8G --time=02:00:00 --exclude=${EX} --job-name=f2m_setup \
      --output=${ROOT}/logs/%x_%j.out --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${R6},F2T_SOURCE=misseg \
      ${STAGE}/hpc/f2t_job.sh cells tasks dispersion)
    record ${j} f2m_setup "cells tasks dispersion"; echo "setup ${j} monitor $(monitor ${j})" ;;
  fit)
    arm=${2:?}; sims=${3:?}; tf=${ROOT}/tasks/tasks_${arm}_${sims}.tsv; n=${4:-$(($(wc -l < ${tf}) - 1))}
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${FIT_MEM:-12G} --time=12:00:00 --array=${ARRAY:-1-${n}} $(dep) \
      --exclude=${EX} --job-name=f2m_fit_${arm}_${sims} --output=${ROOT}/logs/%x_%A_%a.out \
      --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${ROOT}/${arm},ABM_CELLS=${ROOT}/salehi_cells.tsv,ABM_FORECAST_HORIZON=3,MISSEG_FORECAST=TRUE,TASK_FILE=${tf},BENCH_TASKS=none \
      ${STAGE}/hpc/f2t_job.sh fit)
    record ${j} f2m_fit_${arm}_${sims} "${ARRAY:-1-${n}}"; echo "fit ${j} monitor $(monitor ${j})" ;;
  rows)
    arm=${2:?}; sims=${3:?}; n=${4:?}
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${ROWS_MEM:-16G} --time=12:00:00 --array=${ARRAY:-1-${n}} $(dep) \
      --exclude=${EX} --job-name=f2m_rows_${arm}_${sims} --output=${ROOT}/logs/%x_%A_%a.out \
      --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},STAGE=${STAGE},F2T_ROOT=${ROOT},BENCH_RESULTS=${ROOT}/${arm},ABM_CELLS=${ROOT}/salehi_cells.tsv,ABM_FORECAST_HORIZON=3,F2T_SOURCE=misseg,F2T_ARM=${arm},F2T_SIMS=${sims} \
      ${STAGE}/hpc/f2t_job.sh rows)
    record ${j} f2m_rows_${arm}_${sims} "${ARRAY:-1-${n}}"; echo "rows ${j} monitor $(monitor ${j})" ;;
esac
