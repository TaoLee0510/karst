#!/usr/bin/env bash
# One array task of scripts/fig2d_check.R inside the image: STAGE (salehi |
# abm) and the array index as the chunk. Submitted with sbatch from the
# repository root, e.g.
#   sbatch --qos=xxlarge --exclude=... --array=1-28 --export=ALL,BENCH_DIR=$PWD,STAGE=salehi hpc/fig2d_check_task.sh
set -euo pipefail
cd "${BENCH_DIR}"
source hpc/env.sh
load_hpc_environment
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${SLURM_ARRAY_TASK_ID}  stage ${STAGE}  $(date --iso-8601=seconds)"
container_rscript scripts/fig2d_check.R --stage="${STAGE}" --index="${SLURM_ARRAY_TASK_ID}"
echo "Finished $(date --iso-8601=seconds)"
