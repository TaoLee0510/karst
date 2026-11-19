#!/usr/bin/env bash
# One task of the history choice study on HPC (Slurm array task), in the
# ALFA-K paper benchmark's container with the karst library of its lock.
#   sbatch --array=1-N ... hpc_task.sh dev|abm TASK_FILE OUT_DIR
set -euo pipefail
part="$1"
task_file="$2"
out_dir="$3"
cd /share/lab_crd/taoli/Project/karst/benchmark/paper
export BENCH_DIR="$PWD"
source hpc/env.sh
load_hpc_environment
here="${KARST_DIR}/validation/sparse_joint/history_choice"
id="${SLURM_ARRAY_TASK_ID}"
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  ${part}  $(date --iso-8601=seconds)"
echo "Git $(git -C "${KARST_DIR}" rev-parse --short HEAD 2>/dev/null || echo unknown)  lock ${KARST_COMMIT}"
if [[ "${part}" == dev ]]; then
  dataset_id="$(awk -F'\t' -v row=$((id + 1)) 'NR == row {print $1}' "${task_file}")"
  CONTAINER_PWD="${KARST_DIR}" container_rscript "${here}/dev_task.R" "${dataset_id}" "${out_dir}"
else
  container_rscript "${here}/abm_task.R" "${task_file}" "${id}" "${out_dir}" "${here}"
fi
echo "Finished $(date --iso-8601=seconds)"
