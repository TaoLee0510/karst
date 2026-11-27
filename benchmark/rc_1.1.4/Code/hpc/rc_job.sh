#!/usr/bin/env bash
# One array task of the rc and Kester scripts (Code/rc, Code/kester) inside the
# image: Rscript Code/<RC_SCRIPT> TASK_FILE TASK_ID RC_OUT, run in the script's
# folder, with karst 1.1.4 rc (Code/env/lib/karst_50e5889) unless KARST_LIB
# says otherwise. Paths come from Code/rc_env.sh; an image that does not start
# is retried (up to 5 times).
#   sbatch ... --export=ALL,RC_HOME=...,RC_SCRIPT=rc/own_task.R,TASK_FILE=...,RC_OUT=... rc_job.sh
# With RC_ARGS set (submit_rc_args.sh), the script gets those arguments
# instead (space separated), e.g. the analysis scripts of the patch studies.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "${RC_HOME:?}/Code/rc_env.sh"
KARST_LIB="${KARST_LIB:-${RC_CODE}/env/lib/karst_50e5889}"
EXTRA_LIB="${RC_CODE}/env/lib/extra_r"
export PAPER_DATA="${RC_DATA}/paper" PAPER_REF="${RC_DATA}/reference/paper"
export STT_DATA_DIR="${STT_DATA_DIR:-${RC_DATA}/own_benchmark}" KESTER="${RC_DATA}/kester"
# the own benchmark's x1 fits (round 6ae2a24): the rerun's by default; the
# frozen ones with OWN_ROUND_FITS=${RC_FROZEN}/part2_own/round_6ae2a24/fits
export OWN_ROUND_FITS="${OWN_ROUND_FITS:-${RC_RESULTS}/part2_own/round_6ae2a24/fits}"
export OR_CORES="${OR_CORES:-${SLURM_CPUS_PER_TASK:-1}}"
id=${SLURM_ARRAY_TASK_ID:-1}
script="${RC_CODE}/${RC_SCRIPT:?}"
if [[ -n "${RC_ARGS:-}" ]]; then read -r -a script_args <<< "${RC_ARGS}"
else mkdir -p "${RC_OUT:?}"; script_args=("${TASK_FILE:?}" "${id}" "${RC_OUT}"); fi
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  ${RC_SCRIPT}  karst lib ${KARST_LIB}  $(date --iso-8601=seconds)"
for attempt in 1 2 3 4 5; do
  status=0; started=$(date +%s)
  KARST_LIB="${KARST_LIB}" EXTRA_LIB="${EXTRA_LIB}" CONTAINER_PWD="$(dirname "${script}")" \
  RC_PASS_ENV="PAPER_DATA PAPER_REF STT_DATA_DIR KESTER OWN_ROUND_FITS OR_CORES RATE_GRID" RC_PASS_ENV_IF_SET="GATE_RESULTS F2T_ROOT F2T_SOURCE F2T_BOOT F2T_CORES" \
    rc_container_rscript "${script}" "${script_args[@]}" || status=$?
  [[ ${status} -ne 0 && ( ${status} -eq 127 || $(( $(date +%s) - started )) -lt 30 ) && ${attempt} -lt 5 ]] || break
  echo "image did not start on $(hostname) (attempt ${attempt}); trying again" >&2; sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
