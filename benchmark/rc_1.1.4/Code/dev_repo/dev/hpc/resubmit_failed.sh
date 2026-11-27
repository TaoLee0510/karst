#!/usr/bin/env bash
# Resubmit the tasks of a finished array that ended OUT_OF_MEMORY, TIMEOUT or
# FAILED, with new resources:
#   dev/hpc/resubmit_failed.sh JOB SCRIPT TASK_FILE MEM TIME
# (DEV_LIB as for submit.sh). The scripts skip tasks whose output exists.
set -euo pipefail
source "${RC_HOME:?}/Code/dev_repo/dev/hpc/env.sh"   # rc freeze: RC_HOME from Code/rc_env.sh (exported)
load_hpc_environment
job="$1"; script="$2"; task_file="$3"; mem="$4"; limit="$5"
ids=$(sacct -j "${job}" -X -n -P --format=JobID,State | awk -F'|' '$2 ~ /OUT_OF_ME|TIMEOUT|FAILED/ {split($1, a, "_"); print a[2]}' | sort -n -u | paste -sd, -)
if [[ -z "${ids}" ]]; then echo "no failed tasks in ${job}"; exit 0; fi
echo "resubmitting tasks ${ids} of ${job} with ${mem} ${limit}"
"${DEV_REPO}/dev/hpc/submit.sh" "${script}" "${task_file}" "${mem}" "${limit}" "${ids}"
