# Settings of the dev studies on HPC (rc freeze version of dev/hpc/env.sh;
# the as-run copy is karst_dev/repo/dev/hpc/env.sh, recorded in
# results_rc/part4_rate_choice/dev_results). Every path comes from
# Code/rc_env.sh (RC_HOME); nothing points outside it.
#   DEV_REPO     this repository copy (Code/dev_repo), read only
#   DEV_RESULTS  where the dev studies write: RC_RESULTS/part4_rate_choice/dev_results
#                (realdata/beta = C5, realdata/noise = C8a, as in results_rc)
#   DEV_BACKUP   the alfak2-era e1-e4 inputs (data/dev/backup)
#   PANCAN_DATA  the NUMBAT-derived karyotype counts (data/pancan)
#   DEV_LIB      the karst build of a job: Code/env/lib/<DEV_LIB>
#                (karst_cd342db = 1.1.2, karst_v4, karst_6ae2a24 = 1.1.3)
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)/rc_env.sh"
DEV_REPO="${RC_CODE}/dev_repo"
DEV_RESULTS="${DEV_RESULTS:-${RC_RESULTS}/part4_rate_choice/dev_results}"
DEV_BACKUP="${RC_DATA}/dev/backup"
PAPER_DATA="${RC_DATA}/paper"
PANCAN_DATA="${RC_DATA}/pancan"
KARST_CLONE="${RC_CODE}/dev_repo"
DEV_SIF="${RC_SIF}"
SIF_RSCRIPT="${RC_RSCRIPT}"
EXTRA_LIB="${RC_CODE}/env/lib/extra_r"
DEV_LIB="${DEV_LIB:-karst_cd342db}"
DEV_EXCLUDE="${DEV_EXCLUDE:-${RC_EXCLUDE}}"
# The v7 simulations' stored truth paths name the as-run results root; they
# are read from V7_UPSTREAM instead (common_v7.R, v7_check_read): the rerun's
# own DEV_RESULTS, or the frozen one for a rerun of one study alone.
V7_UPSTREAM="${V7_UPSTREAM:-${DEV_RESULTS}}"
export DEV_REPO DEV_RESULTS DEV_BACKUP PAPER_DATA PANCAN_DATA KARST_CLONE DEV_LIB V7_UPSTREAM

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

# The inputs the dev studies read from their results root: the CE c2 grid and
# the dev/realdata/A Salehi table (data/dev/realdata), and the results of a
# study that is not being rerun (STUDY ... from RC_UPSTREAM, e.g. beta for C8a
# alone), linked into DEV_RESULTS (a rerun folder, never results_rc).
dev_link_inputs() {
  mkdir -p "${DEV_RESULTS}/realdata"
  local x
  for x in A CE; do [[ -e "${DEV_RESULTS}/realdata/${x}" ]] || ln -s "${RC_DATA}/dev/realdata/${x}" "${DEV_RESULTS}/realdata/${x}"; done
  for x in "$@"; do
    [[ -e "${DEV_RESULTS}/realdata/${x}" ]] || ln -s "${RC_UPSTREAM}/part4_rate_choice/dev_results/realdata/${x}" "${DEV_RESULTS}/realdata/${x}"
  done
}

# Rscript inside the image (rc_container_rscript) with Code/env/lib/${DEV_LIB}
# first on the library path, run in DEV_REPO.
container_rscript() {
  DEV_LIBDIR="${RC_CODE}/env/lib/${DEV_LIB}"
  KARST_LIB="${DEV_LIBDIR}" EXTRA_LIB="${EXTRA_LIB}" CONTAINER_PWD="${DEV_REPO}" DEV_LIBDIR="${DEV_LIBDIR}" \
  RC_PASS_ENV="DEV_REPO DEV_RESULTS DEV_BACKUP DEV_LIBDIR PAPER_DATA KARST_CLONE PANCAN_DATA V7_UPSTREAM" \
    rc_container_rscript "$@"
}
