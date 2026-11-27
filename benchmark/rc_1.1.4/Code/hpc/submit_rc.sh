#!/usr/bin/env bash
# Submit an array of an rc or Kester script (one task per row of TASK_FILE, no
# throttle, qos xxlarge, the usual excluded nodes) with its add-hours monitor;
# recorded in RC_RESULTS/logs/jobs.tsv.
#   bash Code/hpc/submit_rc.sh NAME SCRIPT TASK_FILE OUT [ARRAY]
#     SCRIPT relative to Code/ (e.g. rc/own_task.R), OUT absolute.
#     CPUS, MEM (16G), TIME (08:00:00), DEPEND (job ids, colon separated,
#     afterany), KARST_LIB, OWN_ROUND_FITS override.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/rc_env.sh"
name=${1:?}; script=${2:?}; tf=$(readlink -f "${3:?}"); out=${4:?}
n=$(( $(wc -l < "${tf}") - 1 )); array=${5:-1-${n}}
mkdir -p "${RC_RESULTS}/logs" "${out}"
dep=""; [[ -z "${DEPEND:-}" ]] || dep="--dependency=afterany:${DEPEND}"
j=$(sbatch --parsable --qos="${RC_QOS}" --cpus-per-task="${CPUS:-1}" --mem="${MEM:-16G}" --time="${TIME:-08:00:00}" \
  --array="${array}" ${dep} --exclude="${RC_EXCLUDE}" --job-name="rc_${name}" \
  --output="${RC_RESULTS}/logs/%x_%A_%a.out" \
  --export=ALL,RC_HOME="${RC_HOME}",RC_RESULTS="${RC_RESULTS}",RC_SCRIPT="${script}",TASK_FILE="${tf}",RC_OUT="${out}" \
  "${RC_CODE}/hpc/rc_job.sh")
m=$(rc_monitor "${j}")
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${j}" "${m}" "rc_${name}" "${script}" "${array}" "$(date --iso-8601=seconds)" >> "${RC_RESULTS}/logs/jobs.tsv"
echo "${name} ${j} monitor ${m}"
