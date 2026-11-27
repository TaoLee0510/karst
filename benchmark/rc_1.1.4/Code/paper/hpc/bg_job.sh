#!/usr/bin/env bash
# Steps of the Salehi missegregation-rate grid (docs/beta_grid_salehi_6ae2a24.md)
# for one rate's results folder (rc freeze: run from Code/paper):
#   sbatch ... --export=ALL,BENCH_RESULTS=<rate folder>,BENCH_TASKS=<rate task folder> hpc/bg_job.sh STEP
# STEP: fit (TASK_FILE, array task) | score (scripts/score_salehi.R --index,
# array task) | noise_chunk (scripts/fig2d_noise.R --stage=chunk, array task)
# | collect (check_ksets salehi, collect_salehi, analyse_salehi,
# fig2d_noise --stage=analyse).
set -euo pipefail
cd "${RC_HOME:?}/Code/paper"
source hpc/env.sh
load_hpc_environment
export CONTAINER_PWD="${BENCH_DIR}"
echo "Host $(hostname)  job ${SLURM_JOB_ID:-NA} task ${SLURM_ARRAY_TASK_ID:-}  results ${BENCH_RESULTS}  tasks ${BENCH_TASKS:-}  $(date --iso-8601=seconds)"
run() {
  local a s started
  for a in 1 2 3 4 5; do
    s=0; started=$(date +%s)
    container_rscript "$@" || s=$?
    [[ ${s} -ne 0 && ( ${s} -eq 127 || $(( $(date +%s) - started )) -lt 30 ) && ${a} -lt 5 ]] || break
    echo "image did not start (attempt ${a}); trying again" >&2; sleep 30
  done
  return "${s}"
}
for step in "$@"; do
  echo "== ${step} $(date --iso-8601=seconds)"
  case "${step}" in
    fit) run scripts/run_task.R "${TASK_FILE:?}" "${SLURM_ARRAY_TASK_ID}" ;;
    score) run scripts/score_salehi.R --index "${SLURM_ARRAY_TASK_ID}" ;;
    noise_chunk) run scripts/fig2d_noise.R --stage=chunk --index="${SLURM_ARRAY_TASK_ID}" ;;
    collect)
      run scripts/check_ksets.R salehi
      run scripts/collect_salehi.R
      run scripts/analyse_salehi.R
      run scripts/fig2d_noise.R --stage=analyse ;;
    *) echo "unknown step ${step}" >&2; exit 2 ;;
  esac
done
echo "Finished $(date --iso-8601=seconds)"
