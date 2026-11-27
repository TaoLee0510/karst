#!/usr/bin/env bash
# Submit one array: one task per row of TASK_FILE (header excluded), qos
# xxlarge, one CPU, no throttle, plus the add-hours monitor of the array.
#   dev/hpc/submit.sh SCRIPT TASK_FILE MEM TIME [ARRAY]
# SCRIPT is relative to DEV_REPO. Env: DEV_LIB (default karst_base),
# SCRIPT_EXTRA (more arguments), DEPEND (job ids, afterok), JOB_NAME.
# Every job is recorded in ${DEV_RESULTS}/jobs.tsv.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/env.sh"   # rc freeze
load_hpc_environment
script="$1"; task_file="$(realpath "$2")"; mem="$3"; limit="$4"
n=$(($(wc -l <"${task_file}") - 1))
array="${5:-1-${n}}"
name="${JOB_NAME:-dev_$(basename "${task_file}" .tsv)}"
mkdir -p "${DEV_RESULTS}/logs"
extra=(--exclude="${DEV_EXCLUDE}")
[[ -n "${DEPEND:-}" ]] && extra+=(--dependency="afterok:${DEPEND}")
job=$(sbatch --parsable "${extra[@]}" --qos="${RC_QOS}" --ntasks=1 --cpus-per-task=1 --mem="${mem}" --time="${limit}" \
  --array="${array}" --job-name="${name}" --output="${DEV_RESULTS}/logs/%x_%A_%a.out" \
  --export="ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},DEV_RESULTS=${DEV_RESULTS},DEV_REPO=${DEV_REPO},DEV_LIB=${DEV_LIB},SCRIPT=${script},TASK_FILE=${task_file},SCRIPT_EXTRA=${SCRIPT_EXTRA:-}" \
  "${DEV_REPO}/dev/hpc/array_task.sh")
job="${job%%;*}"
# The monitor must outlive the tasks it extends: qos small (MaxWall 90 days), 3 days.
monitor=$(rc_monitor "${job}" 2>/dev/null || echo NA)
monitor="${monitor%%;*}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${job}" "${monitor}" "${name}" "${script}" "${task_file}" "${n}" \
  "${mem}" "${limit}" "$(date --iso-8601=seconds)" >>"${DEV_RESULTS}/jobs.tsv"
echo "array ${job} (${n} tasks, ${mem}, ${limit}, lib ${DEV_LIB}); monitor ${monitor}"
