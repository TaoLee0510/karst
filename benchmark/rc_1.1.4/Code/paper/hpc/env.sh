# Settings shared by the paper benchmark's HPC scripts (rc freeze version).
# Every path comes from Code/rc_env.sh (RC_HOME); nothing points outside it.
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"

BENCH_DIR="${RC_CODE}/paper"
# The inputs (raw ABM / varyMisseg simulations, Salehi count matrices,
# lineages, metadata) and the paper's published ALFA-K outputs (reference).
PAPER_DATA="${RC_DATA}/paper"
PAPER_REF="${RC_DATA}/reference/paper"
BENCH_SIF="${RC_SIF}"
# KARST_COMMIT and KARST_SHA256 of the karst source archive: 1.1.4 rc
# (50e5889) by default; Code/env/locks/karst_6ae2a24.lock is the 1.1.3
# library the 6ae2a24 rounds were fitted with (default fits identical).
KARST_LOCK="${KARST_LOCK:-${RC_CODE}/env/locks/karst_50e5889.lock}"
source "${KARST_LOCK}"
KARST_LIB="${KARST_LIB:-${RC_CODE}/env/lib/karst_${KARST_COMMIT}}"
# CRAN packages the image lacks, for the analyses of the paper's Fig. 4-6.
EXTRA_LIB="${RC_CODE}/env/lib/extra_r"
# A run's results folder: the round under RC_RESULTS by default. For a rate
# of the Salehi grid, BENCH_RESULTS=${RC_RESULTS}/part1_paper/salehi_rate_grid/x<mm>.
BENCH_RESULTS="${BENCH_RESULTS:-${RC_RESULTS}/part1_paper/round_6ae2a24}"
BENCH_LOGS="${BENCH_LOGS:-${RC_RESULTS}/logs/paper}"
# The round's task tables (Code/paper/tasks/6ae2a24: KARST; the rate grid's in
# tasks/salehi_rate_grid/x<mm>, version B's in tasks/align_karst).
BENCH_TASKS="${BENCH_TASKS:-${BENCH_DIR}/tasks/6ae2a24}"
SIF_RSCRIPT="${RC_RSCRIPT}"

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
  mkdir -p "${BENCH_LOGS}"
}

# Rscript inside the image (rc_container_rscript: Code, data and results_rc
# read only, RC_RESULTS writable), run in CONTAINER_PWD (default BENCH_DIR).
container_rscript() {
  export PAPER_DATA PAPER_REF BENCH_RESULTS BENCH_SIF FIG456_EXCLUDE="${FIG456_EXCLUDE:-1}" \
    COMPARE_METHOD="${COMPARE_METHOD:-karst}" COMPARE_OLD_METHOD="${COMPARE_OLD_METHOD:-${COMPARE_METHOD:-karst}}" \
    ABM_FORECAST_HORIZON="${ABM_FORECAST_HORIZON:-10}" ALLOW_MISSING="${ALLOW_MISSING:-0}"
  KARST_LIB="${KARST_LIB}" EXTRA_LIB="${EXTRA_LIB}" CONTAINER_PWD="${CONTAINER_PWD:-${BENCH_DIR}}" \
  RC_PASS_ENV="PAPER_DATA PAPER_REF BENCH_RESULTS BENCH_SIF FIG456_EXCLUDE BENCH_TASKS ALLOW_MISSING SCORE_METHODS \
COMPARE_METHOD COMPARE_OLD_METHOD STT_OUTPUT_DIR STT_DATA_DIR XV_SOURCE ALIGN ALIGN_RESULTS DEV_EXTEND_ALFAK OWN_CHUNK_DIR \
CHECK_ALL ABM_CELLS ABM_FORECAST_HORIZON F2T_ROOT MISSEG_FORECAST FIG2C_CENSUS_OUT WDIST_TIMEOUT WASSERSTEIN_TIMEOUT" \
    RC_PASS_ENV_IF_SET="F2T_SOURCE F2T_BOOT F2T_CORES F2T_SALEHI_TSV FIG2D_OUT FIG2D_METHODS FIG2D_DESIGN SALEHI_META" \
    rc_container_rscript "$@"
}
