#!/usr/bin/env bash
# Submit one job of a dev script with free-form arguments, plus its monitor.
#   dev/hpc/submit_args.sh SCRIPT MEM TIME ["ARG ..."]      SCRIPT relative to DEV_REPO
# Env: DEV_LIB, DEV_RESULTS, DEPEND (job ids, afterany), JOB_NAME.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/env.sh"
load_hpc_environment
script="$1"; mem="$2"; limit="$3"; args="${4:-}"
name="${JOB_NAME:-dev_$(basename "${script}" .R)}"
mkdir -p "${DEV_RESULTS}/logs"
dep=(); [[ -z "${DEPEND:-}" ]] || dep=(--dependency="afterany:${DEPEND}")
job=$(sbatch --parsable ${dep[@]+"${dep[@]}"} --exclude="${DEV_EXCLUDE}" --qos="${RC_QOS}" --ntasks=1 --cpus-per-task=1 --mem="${mem}" \
  --time="${limit}" --job-name="${name}" --output="${DEV_RESULTS}/logs/%x_%j.out" \
  --export="ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},DEV_RESULTS=${DEV_RESULTS},DEV_LIB=${DEV_LIB},SCRIPT=${script},SCRIPT_ARGS=${args}" \
  "${DEV_REPO}/dev/hpc/args_task.sh")
monitor=$(rc_monitor "${job}" 2>/dev/null || echo NA)
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${job}" "${monitor}" "${name}" "${script}" "${args}" "$(date --iso-8601=seconds)" >>"${DEV_RESULTS}/jobs.tsv"
echo "job ${job} (${script}, lib ${DEV_LIB}); monitor ${monitor}"
