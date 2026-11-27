#!/usr/bin/env bash
# One array task of the 1.1.4 release-candidate work (PLAN.md): Rscript SCRIPT
# TASK_FILE TASK_ID OUT inside the benchmark image, with karst 50e5889
# (lib/karst_50e5889 of the HPC clone, read only). Writable: RC_ROOT and, for
# Kester, its output folder. An image that does not start is retried.
#   sbatch ... --export=ALL,RC_SCRIPT=rc_check.R,TASK_FILE=...,RC_OUT=... job.sh
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
KARST_CLONE=/share/lab_crd/taoli/Project/karst
RC_ROOT=/share/lab_crd/taoli/Project/karst_rc_1.1.4
KESTER=/share/lab_crd/taoli/Project/Kester
KESTER_OUT=${KESTER}/results/forecast_50e5889
PAPER_DATA=/share/lab_crd/M010_ALFAK_2023/ALFA-K/data
SIF=/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif
RSCRIPT=/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript
KARST_LIB=${KARST_CLONE}/benchmark/paper/lib/karst_50e5889
EXTRA_LIB=${KARST_CLONE}/benchmark/paper/lib/extra_r
mkdir -p "${KESTER_OUT}"
id=${SLURM_ARRAY_TASK_ID:?}
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  ${RC_SCRIPT}  $(date --iso-8601=seconds)"
for attempt in 1 2 3 4 5; do
  status=0; started=$(date +%s)
  apptainer exec --cleanenv --pwd "${RC_ROOT}/code" \
    --bind "${KARST_CLONE}:${KARST_CLONE}:ro" --bind "${RC_ROOT}:${RC_ROOT}" --bind "${PAPER_DATA}:${PAPER_DATA}:ro" \
    --bind "${KESTER}:${KESTER}:ro" --bind "${KESTER_OUT}:${KESTER_OUT}" \
    --env "KARST_LIB=${KARST_LIB}" --env "KARST_CLONE=${KARST_CLONE}" --env "EXTRA_LIB=${EXTRA_LIB}" --env "PAPER_DATA=${PAPER_DATA}" \
    --env "KESTER=${KESTER}" --env "KESTER_OUT=${KESTER_OUT}" --env "SLURM_ARRAY_TASK_ID=${id}" \
    --env "OMP_NUM_THREADS=1" --env "OPENBLAS_NUM_THREADS=1" --env "MKL_NUM_THREADS=1" \
    "${SIF}" /bin/bash -c '
      source /usr/local/lib/pancankfls-runtime-env.sh
      export R_LIBS_USER="${KARST_LIB}:${R_LIBS_USER:-}:${EXTRA_LIB}"
      exec '"${RSCRIPT}"' --vanilla "$@"' rscript "${RC_ROOT}/code/${RC_SCRIPT:?}" "${TASK_FILE:?}" "${id}" "${RC_OUT:?}" || status=$?
  [[ ${status} -ne 0 && ( ${status} -eq 127 || $(( $(date +%s) - started )) -lt 30 ) && ${attempt} -lt 5 ]] || break
  echo "image did not start on $(hostname) (attempt ${attempt}); trying again" >&2; sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
