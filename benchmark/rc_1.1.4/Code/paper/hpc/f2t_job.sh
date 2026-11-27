#!/usr/bin/env bash
# Steps of the Fig. 2d truth study (scripts/fig2d_truth.R,
# docs/fig2d_truth_6ae2a24.md), run from a staging copy of benchmark/paper
# (R/, scripts/, hpc/) inside the clone, as hpc/xv_job.sh of the
# cross-validation study: BENCH_DIR stays the clone (libraries, binds; KARST
# from the round's lock, lib/karst_6ae2a24), the image runs in the staging
# folder.
#   sbatch ... --export=ALL,STAGE=<staging>,F2T_ROOT=<root>,BENCH_RESULTS=<fits> f2t_job.sh STEP ...
# STEP: cells | tasks | engine | dispersion | fit (TASK_FILE, array task;
# ABM_CELLS, ABM_FORECAST_HORIZON and MISSEG_FORECAST from the environment) |
# rows (F2T_STEP or, with F2T_SOURCE=misseg, F2T_ARM; F2T_SIMS, array task).
# F2T_SOURCE: abm (default) or misseg (docs/fig2d_misseg_6ae2a24.md).
set -euo pipefail
# rc freeze: run from Code/paper (STAGE is no longer used)
cd "${RC_HOME:?}/Code/paper"
source hpc/env.sh
load_hpc_environment
export CONTAINER_PWD="${BENCH_DIR}"
echo "Host $(hostname)  job ${SLURM_JOB_ID:-NA} task ${SLURM_ARRAY_TASK_ID:-}  results ${BENCH_RESULTS}  ABM_CELLS=${ABM_CELLS:-}  $(date --iso-8601=seconds)"
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
    cells|tasks|engine|dispersion) run scripts/fig2d_truth.R --stage="${step}" --source="${F2T_SOURCE:-abm}" ;;
    fit) run scripts/run_task.R "${TASK_FILE:?}" "${SLURM_ARRAY_TASK_ID}" ;;
    rows)
      if [[ "${F2T_SOURCE:-abm}" == misseg ]]; then
        run scripts/fig2d_truth.R --stage=rows --source=misseg --arm="${F2T_ARM:?}" --sims="${F2T_SIMS:?}" \
          --index="${SLURM_ARRAY_TASK_ID}"
      else
        run scripts/fig2d_truth.R --stage=rows --step="${F2T_STEP:?}" --sims="${F2T_SIMS:?}" \
          --index="${SLURM_ARRAY_TASK_ID}"
      fi ;;
    *) echo "unknown step ${step}" >&2; exit 2 ;;
  esac
done
echo "Finished $(date --iso-8601=seconds)"
