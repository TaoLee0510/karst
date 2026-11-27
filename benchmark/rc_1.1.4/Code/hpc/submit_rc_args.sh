#!/usr/bin/env bash
# Submit one job of an rc / Kester / patch script with free-form arguments
# (Code/hpc/rc_job.sh with RC_ARGS), with its add-hours monitor; recorded in
# RC_RESULTS/logs/jobs.tsv.
#   bash Code/hpc/submit_rc_args.sh NAME SCRIPT "ARG ..."      SCRIPT relative to Code/
#     CPUS, MEM (16G), TIME (08:00:00), DEPEND, KARST_LIB override.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/rc_env.sh"
name=${1:?}; script=${2:?}; args=${3:?}
mkdir -p "${RC_RESULTS}/logs"
dep=""; [[ -z "${DEPEND:-}" ]] || dep="--dependency=afterany:${DEPEND}"
j=$(sbatch --parsable --qos="${RC_QOS}" --cpus-per-task="${CPUS:-1}" --mem="${MEM:-16G}" --time="${TIME:-08:00:00}" \
  ${dep} --exclude="${RC_EXCLUDE}" --job-name="rc_${name}" --output="${RC_RESULTS}/logs/%x_%j.out" \
  --export=ALL,RC_HOME="${RC_HOME}",RC_RESULTS="${RC_RESULTS}",RC_UPSTREAM="${RC_UPSTREAM}",RC_SCRIPT="${script}",RC_ARGS="${args}",OR_CORES="${CPUS:-1}" \
  "${RC_CODE}/hpc/rc_job.sh")
m=$(rc_monitor "${j}")
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${j}" "${m}" "rc_${name}" "${script}" "args" "$(date --iso-8601=seconds)" >> "${RC_RESULTS}/logs/jobs.tsv"
echo "${name} ${j} monitor ${m}"
