#!/usr/bin/env bash
# Submit an array of the 1.1.4 release-candidate work with its monitor.
#   bash submit.sh NAME SCRIPT TASK_FILE OUT [ARRAY]   (MEM, TIME, DEPEND override)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_ROOT=/share/lab_crd/taoli/Project/karst_rc_1.1.4
EX=hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063,hpctpa3pc0062,hpctpa3pc0028,hpctpa3pc0021,hpctpa3pc0024,hpctpa3pc0023,hpctpa3pc0016,hpctpa3pc0022,hpctpa3pc0061,hpctpa3pc0003
name=${1:?}; script=${2:?}; tf=${3:?}; out=${4:?}
n=$(( $(wc -l < ${tf}) - 1 )); array=${5:-1-${n}}
mkdir -p ${RC_ROOT}/logs ${out}
dep=""; [[ -z "${DEPEND:-}" ]] || dep="--dependency=afterany:${DEPEND}"
j=$(sbatch --parsable --qos=xxlarge --cpus-per-task=${CPUS:-1} --mem=${MEM:-16G} --time=${TIME:-08:00:00} --array=${array} ${dep} \
  --exclude=${EX} --job-name=rc_${name} --output=${RC_ROOT}/logs/%x_%A_%a.out \
  --export=ALL,RC_SCRIPT=${script},TASK_FILE=${tf},RC_OUT=${out} ${RC_ROOT}/code/hpc/${JOB_SH:-job.sh})
m=$(sbatch --parsable --qos=small --time=3-00:00:00 --export=ALL,THRESHOLD_HOURS=3,DRY_RUN=FALSE,AUTO_CONFIRM=TRUE ~/submit_hourly_add_hours.sh ${j} 3)
printf '%s\t%s\t%s\t%s\t%s\n' "${j}" "${m}" "rc_${name}" "${array}" "$(date --iso-8601=seconds)" >> ${RC_ROOT}/jobs.tsv
echo "${name} ${j} monitor ${m}"
