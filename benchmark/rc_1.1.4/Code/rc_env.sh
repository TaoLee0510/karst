# Paths of the KARST 1.1.4 rc freeze, sourced by every script in Code/.
#
# RC_HOME is the rc_1.1.4 folder (this file's parent's parent). Nothing here
# points outside it:
#   Code/       this code (read only inside the image)
#   Code/env/   the image (SIF), karst sources, locks and built libraries
#   data/       the inputs (read only inside the image)
#   results_rc/ the frozen results (read only inside the image; never written)
#   RC_RESULTS  where a run writes, default RC_HOME/rerun (same layout as
#               results_rc/, so a rerun is checked by comparing the two)
RC_HOME="${RC_HOME:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RC_CODE="${RC_HOME}/Code"
RC_DATA="${RC_HOME}/data"
RC_FROZEN="${RC_HOME}/results_rc"
RC_RESULTS="${RC_RESULTS:-${RC_HOME}/rerun}"
# Where a study reads another study's results (e.g. the paper round's fits
# for fig2d_truth step 1): the run's own output root by default (a full rerun
# in order); RC_UPSTREAM=${RC_FROZEN} reruns one study on the frozen results.
RC_UPSTREAM="${RC_UPSTREAM:-${RC_RESULTS}}"
RC_SIF="${RC_SIF:-${RC_CODE}/env/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif}"
RC_RSCRIPT=/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript
# Slurm: qos and the nodes left out (container start failures seen there)
RC_QOS="${RC_QOS:-xxlarge}"
RC_EXCLUDE="${RC_EXCLUDE:-hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063,hpctpa3pc0062,hpctpa3pc0028,hpctpa3pc0021,hpctpa3pc0024,hpctpa3pc0023,hpctpa3pc0016,hpctpa3pc0022,hpctpa3pc0061,hpctpa3pc0003}"
export RC_HOME RC_CODE RC_DATA RC_FROZEN RC_RESULTS RC_UPSTREAM RC_SIF RC_RSCRIPT RC_QOS RC_EXCLUDE
case "$(realpath -m "${RC_RESULTS}")/" in
  "$(realpath -m "${RC_FROZEN}")"/*) echo "RC_RESULTS must not be inside results_rc (the frozen results)" >&2
                                      return 1 2>/dev/null || exit 1 ;;
esac

# The add-hours monitor of a Slurm job (Code/monitor), qos small, 3 days.
rc_monitor() {
  mkdir -p "${RC_RESULTS}/logs"
  sbatch --parsable --qos=small --time=3-00:00:00 --output="${RC_RESULTS}/logs/monitor_%j.out" --error="${RC_RESULTS}/logs/monitor_%j.err" \
    --export=ALL,THRESHOLD_HOURS=3,DRY_RUN=FALSE,AUTO_CONFIRM=TRUE,ADD_HOURS_SCRIPT="${RC_CODE}/monitor/add_hours_to_running_array_jobs.sh" \
    "${RC_CODE}/monitor/submit_hourly_add_hours.sh" "$1" 3
}

# Rscript inside the image: Code, data and results_rc read only, RC_RESULTS
# writable; KARST_LIB first on the library path, EXTRA_LIB last. Runs in
# CONTAINER_PWD. Variables named in RC_PASS_ENV (space separated) are passed in
# (empty when unset, as the as-run scripts did); those in RC_PASS_ENV_IF_SET
# only when they are set.
# A CONTAINER_PWD inside Code/ is replaced by the same folder of a mirror of
# Code/ (RC_RESULTS/work/Code: folders real, files symbolic links to Code/,
# Code/env left out), so that a script reads Code/ through relative paths but
# whatever it writes into its working folder (e.g. a plot device's
# Rplots.pdf) lands in RC_RESULTS. Remove RC_RESULTS/work after adding files to Code/.
rc_container_rscript() {
  mkdir -p "${RC_RESULTS}/work"
  local envs=() v pwd="${CONTAINER_PWD:?}" farm="${RC_RESULTS}/work/Code"
  case "${pwd}/" in
    "${RC_CODE}/"*)
      ( flock 9
        if [[ ! -e "${farm}/.built" ]]; then
          rm -rf "${farm}"; mkdir -p "${farm}"
          for v in "${RC_CODE}"/*; do [[ "$(basename "${v}")" == env ]] || cp -rs "${v}" "${farm}/"; done
          touch "${farm}/.built"
        fi ) 9> "${RC_RESULTS}/work/.lock"
      pwd="${farm}${pwd#"${RC_CODE}"}" ;;
  esac
  for v in RC_HOME RC_RESULTS RC_UPSTREAM ${RC_PASS_ENV:-}; do envs+=(--env "${v}=${!v:-}"); done
  for v in ${RC_PASS_ENV_IF_SET:-}; do [[ -z "${!v+x}" ]] || envs+=(--env "${v}=${!v}"); done
  apptainer exec --cleanenv --pwd "${pwd}" \
    --bind "${RC_CODE}:${RC_CODE}:ro" --bind "${RC_DATA}:${RC_DATA}:ro" \
    --bind "${RC_FROZEN}:${RC_FROZEN}:ro" --bind "${RC_RESULTS}:${RC_RESULTS}" \
    "${envs[@]}" --env "KARST_LIB=${KARST_LIB:?}" --env "EXTRA_LIB=${EXTRA_LIB:-}" \
    --env "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}" --env "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-}" \
    --env "OMP_NUM_THREADS=1" --env "OPENBLAS_NUM_THREADS=1" --env "MKL_NUM_THREADS=1" \
    "${RC_SIF}" /bin/bash -c '
      source /usr/local/lib/pancankfls-runtime-env.sh
      export R_LIBS_USER="${KARST_LIB}:${R_LIBS_USER:-}${EXTRA_LIB:+:${EXTRA_LIB}}"
      exec '"${RC_RSCRIPT}"' --vanilla "$@"' rscript "$@"
}
