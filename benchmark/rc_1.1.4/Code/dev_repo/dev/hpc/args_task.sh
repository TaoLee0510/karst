#!/usr/bin/env bash
# One job of a dev script with free-form arguments (SCRIPT_ARGS, space
# separated; often none, e.g. the evaluation scripts) inside the image, with
# Code/env/lib/${DEV_LIB}. rc freeze addition (the as-run evaluations were run
# with a local Rscript); submitted by submit_args.sh.
set -euo pipefail
source "${RC_HOME:?}/Code/dev_repo/dev/hpc/env.sh"
load_hpc_environment
echo "Host $(hostname)  job ${SLURM_JOB_ID:-NA}  script ${SCRIPT}  lib ${DEV_LIB}  $(date --iso-8601=seconds)"
read -r -a extra <<<"${SCRIPT_ARGS:-}"
# as array_task.sh: a failure within 30 s is taken as the image not starting
# (status 127, or "unknown userid" on some nodes) and retried; a quick R error
# only repeats, since a rerun overwrites the same outputs
for attempt in 1 2 3 4 5; do
  status=0
  started=$(date +%s)
  container_rscript "${SCRIPT}" ${extra[@]+"${extra[@]}"} || status=$?
  quick=$(( $(date +%s) - started < 30 ))
  [[ ${status} -ne 0 && ( ${status} -eq 127 || ${quick} -eq 1 ) && ${attempt} -lt 5 ]] || break
  echo "The image did not start on $(hostname) (attempt ${attempt}); trying again" >&2
  sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
