#!/usr/bin/env bash
# One array task of the dispersion gate (PLAN.md) on the HPC: Rscript
# GATE_SCRIPT TASK_FILE TASK_ID [GATE_ARGS...] inside the image, with the
# round library (karst 1.1.3, Code/env/lib/karst_6ae2a24), run in Code/dev_repo
# (the repository tree the as-run copy results/stage1 was). rc freeze version:
# paths from Code/rc_env.sh. An image that does not start ("unknown userid",
# status 127, or a failure within 30 s) is retried; the scripts skip finished work.
#   sbatch ... --export=ALL,RC_HOME=...,GATE_SCRIPT=...,TASK_FILE=...,GATE_ARGS="..." gate_job.sh
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "${RC_HOME:?}/Code/rc_env.sh"
export KARST_LIB="${KARST_LIB:-${RC_CODE}/env/lib/karst_6ae2a24}" EXTRA_LIB="${RC_CODE}/env/lib/extra_r"
id="${SLURM_ARRAY_TASK_ID}"
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  script ${GATE_SCRIPT}  karst lib ${KARST_LIB}  $(date --iso-8601=seconds)"
read -r -a extra <<<"${GATE_ARGS:-}"
for attempt in 1 2 3 4 5; do
  status=0
  started=$(date +%s)
  CONTAINER_PWD="${RC_CODE}/dev_repo" RC_PASS_ENV="V7_UPSTREAM" \
    rc_container_rscript "${GATE_SCRIPT:?}" "${TASK_FILE:?}" "${id}" ${extra[@]+"${extra[@]}"} || status=$?
  quick=$(( $(date +%s) - started < 30 ))
  [[ ${status} -ne 0 && ( ${status} -eq 127 || ${quick} -eq 1 ) && ${attempt} -lt 5 ]] || break
  echo "The image did not start on $(hostname) (attempt ${attempt}); trying again" >&2
  sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
