#!/usr/bin/env bash
# History choice study (2026-11-19): every task, one Slurm array per memory
# tier, each with an add-hours monitor. Registry: results/jobs.tsv.
set -euo pipefail
cd /share/lab_crd/taoli/Project/karst/benchmark/paper
HC=/share/lab_crd/taoli/Project/karst/validation/sparse_joint/history_choice
R="${HC}/results"
X=hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063,hpctpa3pc0062
reg="${R}/jobs.tsv"
monitor() { sbatch --parsable --qos=xxlarge --time=12:00:00 --export=ALL,THRESHOLD_HOURS=3 ~/submit_hourly_add_hours.sh "$1" 3; }
log() { printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$(date --iso-8601=seconds)" >>"${reg}"; echo "$2 $1 tasks $3 (monitor $4)"; }
n=$(($(wc -l <"${R}/dev_tasks.tsv") - 1))
j=$(sbatch --parsable --qos=xxlarge --cpus-per-task=1 --mem=4G --time=02:00:00 --exclude="${X}" --array=1-${n} \
  --job-name=hc_dev --output=logs/%x_%A_%a.out "${HC}/hpc_task.sh" dev "${R}/dev_tasks.tsv" "${R}/dev")
log "$j" hc_dev "$n" "$(monitor "$j")"
for mem in 4G 8G 16G 32G; do
  rows=$(awk -F'\t' -v m="${mem}" 'NR > 1 && $6 == m {print NR - 1}' "${R}/abm_tasks.tsv" | paste -sd,)
  [[ -n "${rows}" ]] || continue
  j=$(sbatch --parsable --qos=xxlarge --cpus-per-task=1 --mem="${mem}" --time=12:00:00 --exclude="${X}" --array="${rows}" \
    --job-name="hc_abm_${mem}" --output=logs/%x_%A_%a.out "${HC}/hpc_task.sh" abm "${R}/abm_tasks.tsv" "${R}/abm")
  log "$j" "hc_abm_${mem}" "$(echo "${rows}" | tr ',' '\n' | wc -l)" "$(monitor "$j")"
done
