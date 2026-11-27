#!/usr/bin/env bash
# Submit a copy batch of the rc freeze as one array (one task per list) with its
# add-hours monitor; the job is recorded in Code/freeze/jobs.tsv.
#   bash submit_copy.sh TASK_FILE [ARRAY]
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC=${RC_HOME:-$(cd "$(dirname "$0")/../.." && pwd)}
tf=$(readlink -f "${1:?}"); n=$(( $(wc -l < "$tf") - 1 )); array=${2:-1-${n}}
EX=hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063,hpctpa3pc0062,hpctpa3pc0028,hpctpa3pc0021,hpctpa3pc0024,hpctpa3pc0023,hpctpa3pc0016,hpctpa3pc0022,hpctpa3pc0061,hpctpa3pc0003
mkdir -p "$RC/Code/freeze/logs"
name=freeze_$(basename "$tf" .tsv)
j=$(sbatch --parsable --qos=xxlarge --cpus-per-task=4 --mem=8G --time=12:00:00 --array="${array}" --exclude=${EX} \
  --job-name="${name}" --output="$RC/Code/freeze/logs/%x_%A_%a.out" \
  --export=ALL,RC_HOME="$RC",TASK_FILE="$tf" "$RC/Code/freeze/copy_task.sh")
m=$(sbatch --parsable --qos=small --time=3-00:00:00 --export=ALL,THRESHOLD_HOURS=3,DRY_RUN=FALSE,AUTO_CONFIRM=TRUE \
  ~/submit_hourly_add_hours.sh "${j}" 3)
printf '%s\t%s\t%s\t%s\t%s\n' "${j}" "${m}" "${name}" "${array}" "$(date --iso-8601=seconds)" >> "$RC/Code/freeze/jobs.tsv"
echo "${name} ${j} monitor ${m}"
