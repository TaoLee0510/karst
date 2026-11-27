#!/usr/bin/env bash
# The Salehi missegregation-rate grid (A2; as-run: results_rc/part1_paper/
# salehi_rate_grid/submit_bg.sh), rc freeze version: every rate's folder under
# RC_RESULTS/part1_paper/salehi_rate_grid/x<mm>, its tasks in
# Code/paper/tasks/salehi_rate_grid/x<mm>.
#   bash Code/paper/hpc/submit_rate_grid.sh fit M [ARRAY]     fits (one task per row of tasks_salehi.tsv)
#   bash Code/paper/hpc/submit_rate_grid.sh score M [ARRAY]   scripts/score_salehi.R, one task per trajectory (277)
#   bash Code/paper/hpc/submit_rate_grid.sh noise M [ARRAY]   scripts/fig2d_noise.R chunks (28)
#   bash Code/paper/hpc/submit_rate_grid.sh collect M         check_ksets, collect, analyse, fig2d_noise analyse
# M is the multiple (1, 2, 4, ..., 64). x1's fits are the round's
# (BENCH_FITS_FROM: the round folder; link or copy them into x01/fits first).
# DEPEND=<jobs, colon separated> waits (DEPEND_TYPE, default afterany); MEM overrides the memory.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"
ROOT="${RC_RESULTS}/part1_paper/salehi_rate_grid"
step=${1:?}; m=$(printf '%02d' "$((10#${2:?}))")
RES="${ROOT}/x${m}"; TASKS="${RC_CODE}/paper/tasks/salehi_rate_grid/x${m}"
mkdir -p "${RC_RESULTS}/logs/rate_grid" "${RES}"
dep() { [[ -z "${DEPEND:-}" ]] || echo "--dependency=${DEPEND_TYPE:-afterany}:${DEPEND}"; }
common="--qos=${RC_QOS} --exclude=${RC_EXCLUDE} --output=${RC_RESULTS}/logs/rate_grid/%x_%A_%a.out"
envs="ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},BENCH_RESULTS=${RES},BENCH_TASKS=${TASKS}"
case "${step}" in
  fit)
    n=$(( $(wc -l < "${TASKS}/tasks_salehi.tsv") - 1 )); array=${3:-1-${n}}
    j=$(sbatch --parsable ${common} --cpus-per-task=1 --mem=${MEM:-16G} --time=12:00:00 --array=${array} $(dep) \
      --job-name=bg_fit_x${m} --export=${envs},TASK_FILE=${TASKS}/tasks_salehi.tsv "${RC_CODE}/paper/hpc/bg_job.sh" fit) ;;
  score)
    array=${3:-1-277}
    j=$(sbatch --parsable ${common} --cpus-per-task=1 --mem=${MEM:-8G} --time=06:00:00 --array=${array} $(dep) \
      --job-name=bg_score_x${m} --export=${envs} "${RC_CODE}/paper/hpc/bg_job.sh" score) ;;
  noise)
    array=${3:-1-28}
    j=$(sbatch --parsable ${common} --cpus-per-task=1 --mem=${MEM:-8G} --time=06:00:00 --array=${array} $(dep) \
      --job-name=bg_noise_x${m} --export=${envs} "${RC_CODE}/paper/hpc/bg_job.sh" noise_chunk) ;;
  collect)
    j=$(sbatch --parsable ${common} --cpus-per-task=4 --mem=${MEM:-48G} --time=08:00:00 $(dep) \
      --job-name=bg_collect_x${m} --export=${envs},SLURM_CPUS_PER_TASK=4 "${RC_CODE}/paper/hpc/bg_job.sh" collect) ;;
  *) echo "unknown step" >&2; exit 2 ;;
esac
mon=$(rc_monitor "${j}")
printf '%s\t%s\t%s\t%s\n' "${j}" "bg_${step}_x${m}" "${3:-all}" "$(date --iso-8601=seconds)" >> "${RC_RESULTS}/logs/jobs.tsv"
echo "${step} x${m} ${j} monitor ${mon}"
