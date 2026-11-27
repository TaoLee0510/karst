#!/usr/bin/env bash
# Steps of the cross-validation study (docs/xv_target_6ae2a24.md) run from a
# staging copy of benchmark/paper (R/, scripts/, hpc/) inside the clone, so
# that the clone's tracked files stay as they are until the code is
# committed. BENCH_DIR stays the clone (libraries, binds); the image runs in
# the staging folder (CONTAINER_PWD), whose env.sh passes XV_SOURCE.
#   sbatch ... --export=ALL,STAGE=<staging>,BENCH_RESULTS=<round>,XV_SOURCE=... xv_job.sh STEP ...
# STEP: an R script of scripts/ (without .R); fig4 / fig6 (array task
# SLURM_ARRAY_TASK_ID); fig5_wgd (both FIG456_EXCLUDE settings).
set -euo pipefail
# rc freeze: run from Code/paper (STAGE is no longer used)
cd "${RC_HOME:?}/Code/paper"
source hpc/env.sh
load_hpc_environment
export CONTAINER_PWD="${BENCH_DIR}"
echo "Host $(hostname)  job ${SLURM_JOB_ID:-NA} task ${SLURM_ARRAY_TASK_ID:-}  XV_SOURCE=${XV_SOURCE:-}  results ${BENCH_RESULTS}  $(date --iso-8601=seconds)"
run() {
  local a s
  for a in 1 2 3 4 5; do
    s=0; container_rscript "$@" || s=$?
    [[ ${s} -eq 127 && ${a} -lt 5 ]] || break
    echo "image did not start (attempt ${a}); trying again" >&2; sleep 30
  done
  return "${s}"
}
for step in "$@"; do
  echo "== ${step} $(date --iso-8601=seconds)"
  case "${step}" in
    fig4) run scripts/fig4_cna_effects.R --index "${SLURM_ARRAY_TASK_ID}" ;;
    fig6) run scripts/fig6_steady_states.R --index "${SLURM_ARRAY_TASK_ID}" ;;
    fig5_wgd) FIG456_EXCLUDE=1 run scripts/fig5_wgd.R; FIG456_EXCLUDE=0 run scripts/fig5_wgd.R ;;
    own_chunk) run scripts/xv_own_target.R --index "${SLURM_ARRAY_TASK_ID}" --chunks "${OWN_CHUNKS:?}" ;;
    own_combine) run scripts/xv_own_target.R --combine ;;
    align_make_tasks) run scripts/align_task.R --make-tasks ${ALIGN_PILOT_N:-} ;;
    align_check) run scripts/align_check.R ${ALIGN_CHECK_N:-3} ;;
    align_task) run scripts/align_task.R "${ALIGN_TASKS:?}" "${SLURM_ARRAY_TASK_ID}" ;;
    score_salehi) run scripts/score_salehi.R --index "${SLURM_ARRAY_TASK_ID}" ;;
    check_ksets_abm) run scripts/check_ksets.R abm misseg ;;
    check_ksets_salehi) run scripts/check_ksets.R salehi ;;
    align_summary) run scripts/align_summary.R ;;
    *) run "scripts/${step}.R" ;;
  esac
done
echo "Finished $(date --iso-8601=seconds)"
