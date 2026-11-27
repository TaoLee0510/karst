#!/usr/bin/env bash
# One array task: Rscript SCRIPT TASK_FILE TASK_ID [SCRIPT_EXTRA...] inside the
# image, with lib/${DEV_LIB}. Submitted by dev/hpc/submit.sh. A start failure
# of the image on some nodes ("unknown userid", status 127) is retried.
set -euo pipefail
source "${RC_HOME:?}/Code/dev_repo/dev/hpc/env.sh"   # rc freeze: RC_HOME from Code/rc_env.sh (exported)
load_hpc_environment
id="${SLURM_ARRAY_TASK_ID}"
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  script ${SCRIPT}  lib ${DEV_LIB}  $(date --iso-8601=seconds)"
read -r -a extra <<<"${SCRIPT_EXTRA:-}"
# A failure within 30 s is taken as the image not starting (status 127, or
# other codes with "unknown userid" on some nodes) and retried; the scripts skip
# finished work, so a quick R error only repeats harmlessly.
for attempt in 1 2 3 4 5; do
  status=0
  started=$(date +%s)
  container_rscript "${SCRIPT}" "${TASK_FILE}" "${id}" ${extra[@]+"${extra[@]}"} || status=$?
  quick=$(( $(date +%s) - started < 30 ))
  [[ ${status} -ne 0 && ( ${status} -eq 127 || ${quick} -eq 1 ) && ${attempt} -lt 5 ]] || break
  echo "The image did not start on $(hostname) (attempt ${attempt}); trying again" >&2
  sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
