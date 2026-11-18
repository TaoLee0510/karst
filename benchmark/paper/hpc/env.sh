# Settings shared by the HPC scripts; sourced from the repository root.

BENCH_DIR="${BENCH_DIR:-/share/lab_crd/taoli/Project/karst/benchmark/paper}"
# The karst clone this folder belongs to (benchmark/paper of it): the
# development benchmark (benchmark/sparse_two_timepoint) and the Salehi
# horizon test (validation/sparse_joint/salehi_horizon) run from it.
KARST_DIR="${KARST_DIR:-$(cd "${BENCH_DIR}/../.." && pwd)}"
# The benchmark of the alfak2 era (KARST's name up to 0.2.0), read only:
# alfak's fits are not refitted but linked from its 8327282 round (same alfak
# dev2@f93b2d2, settings and inputs; hpc/link_alfak.sh), and its later rounds
# are what scripts/compare_rounds.R compares a KARST round with. It is bound
# into the image, so that the links resolve there.
ALFAK2_ERA_DIR="${ALFAK2_ERA_DIR:-/share/lab_crd/taoli/Project/alfak2-paper-benchmark}"
ALFAK_REFERENCE="${ALFAK_REFERENCE:-${ALFAK2_ERA_DIR}/results}"
PAPER_DATA="${PAPER_DATA:-/share/lab_crd/M010_ALFAK_2023/ALFA-K/data}"
# alfakR_dev2 f93b2d2 image used by PANcanKFLs
BENCH_SIF="${BENCH_SIF:-/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif}"
# KARST_COMMIT and KARST_SHA256 of the karst source archive. KARST_LOCK
# names another lock file (e.g. hpc/locks/karst_<commit>.lock), so that a
# round can be prepared while jobs of the round in hpc/karst.lock are still
# queued: they read the lock when they start. Jobs inherit KARST_LOCK from
# the submitting shell (hpc/submit.sh exports ALL).
KARST_LOCK="${KARST_LOCK:-${BENCH_DIR}/hpc/karst.lock}"
source "${KARST_LOCK}"
KARST_LIB="${BENCH_DIR}/lib/karst_${KARST_COMMIT}"
# CRAN packages the image lacks, for the analyses of the paper's Fig. 4-6
# (glmmTMB; hpc/build_extra_r.sh), after karst's library on the path.
EXTRA_LIB="${BENCH_DIR}/lib/extra_r"
# Each karst commit has its own results folder. alfak's fits are linked in
# from ALFAK_REFERENCE (hpc/link_alfak.sh).
BENCH_RESULTS="${BENCH_RESULTS:-${BENCH_DIR}/results_${KARST_COMMIT}}"
SIF_RSCRIPT=/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript

# sbatch and apptainer come from /etc/profile, which does not run under
# `set -eu`; the caller's options are restored afterwards.
load_hpc_environment() {
  if ! command -v sbatch >/dev/null 2>&1 || ! command -v apptainer >/dev/null 2>&1; then
    local options="$-"
    set +eu
    source /etc/profile >/dev/null 2>&1
    [[ "${options}" == *e* ]] && set -e
    [[ "${options}" == *u* ]] && set -u
  fi
  command -v apptainer >/dev/null 2>&1 || { echo "apptainer is unavailable" >&2; return 1; }
}

# Rscript inside the image, with the karst side library first on the
# library path and the extra packages last. Arguments are passed to Rscript;
# it runs in CONTAINER_PWD (default BENCH_DIR).
container_rscript() {
  apptainer exec --cleanenv \
    --pwd "${CONTAINER_PWD:-${BENCH_DIR}}" \
    --bind "${KARST_DIR}:${KARST_DIR}" \
    --bind "${BENCH_DIR}:${BENCH_DIR}" \
    --bind "${PAPER_DATA}:${PAPER_DATA}:ro" \
    --bind "${ALFAK2_ERA_DIR}:${ALFAK2_ERA_DIR}:ro" \
    --env "KARST_LIB=${KARST_LIB}" \
    --env "EXTRA_LIB=${EXTRA_LIB}" \
    --env "PAPER_DATA=${PAPER_DATA}" \
    --env "BENCH_RESULTS=${BENCH_RESULTS}" \
    --env "BENCH_SIF=${BENCH_SIF}" \
    --env "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}" \
    --env "FIG456_EXCLUDE=${FIG456_EXCLUDE:-1}" \
    --env "BENCH_TASKS=${BENCH_TASKS:-}" \
    --env "ALLOW_MISSING=${ALLOW_MISSING:-0}" \
    --env "SCORE_METHODS=${SCORE_METHODS:-}" \
    --env "COMPARE_METHOD=${COMPARE_METHOD:-karst}" \
    --env "COMPARE_OLD_METHOD=${COMPARE_OLD_METHOD:-${COMPARE_METHOD:-karst}}" \
    --env "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}" \
    --env "OMP_NUM_THREADS=1" \
    --env "OPENBLAS_NUM_THREADS=1" \
    --env "MKL_NUM_THREADS=1" \
    "${BENCH_SIF}" /bin/bash -c '
      source /usr/local/lib/pancankfls-runtime-env.sh
      export R_LIBS_USER="${KARST_LIB}:${R_LIBS_USER:-}:${EXTRA_LIB}"
      exec '"${SIF_RSCRIPT}"' --vanilla "$@"' rscript "$@"
}
