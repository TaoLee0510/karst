#!/usr/bin/env bash
# One array task inside the image. Submitted by hpc/submit.sh, which sets
# MODE (run | reproduce | score_salehi | reproduce_salehi | fig4 | fig6 |
# script | dev | far) and, for MODE=run, dev and far,
# TASK_FILE;
# MODE=script runs SCRIPT_ARGS (a script and its arguments) once per task,
# e.g. a single analysis job. MODE=dev: one fit of the development benchmark
# (benchmark/sparse_two_timepoint/run_task.R), run from the karst clone.
set -euo pipefail
cd "${BENCH_DIR}"
source hpc/env.sh
load_hpc_environment
id="${SLURM_ARRAY_TASK_ID}"
echo "Host $(hostname)  job ${SLURM_ARRAY_JOB_ID:-NA}_${id}  mode ${MODE}  $(date --iso-8601=seconds)"
# Some nodes fail now and then to look up the user when the image starts
# ("FATAL: user: unknown userid", exit status 127), before anything runs; the
# task is then tried again (up to 5 times, 30 s apart). R exits with status 1
# on an error, so a failed script is not repeated, and the scripts skip the
# steps already done.
dispatch() {
  case "${MODE}" in
    run) container_rscript scripts/run_task.R "${TASK_FILE}" "${id}" ;;
    reproduce) container_rscript scripts/check_reproduction.R --out "${BENCH_RESULTS}/reproduction" --index "${id}" ;;
    score_salehi) container_rscript scripts/score_salehi.R --index "${id}" ;;
    reproduce_salehi)
      container_rscript scripts/check_reproduction_salehi.R --out "${BENCH_RESULTS}/reproduction_salehi" --index "${id}" ;;
    fig4) container_rscript scripts/fig4_cna_effects.R --index "${id}" ;;
    fig6) container_rscript scripts/fig6_steady_states.R --index "${id}" ;;
    script) read -r -a script_args <<<"${SCRIPT_ARGS}"; container_rscript "${script_args[@]}" ;;
    dev) CONTAINER_PWD="${RC_CODE}/own" container_rscript benchmark/sparse_two_timepoint/run_task.R "${TASK_FILE}" "${id}" ;;
    far)
      read -r split dataset method < <(awk -F'\t' -v i="${id}" 'NR == i + 1 {print $1, $2, $3}' "${TASK_FILE}")
      CONTAINER_PWD="${RC_CODE}/own" container_rscript validation/sparse_joint/far_distance.R --split="${split}" \
        --ids="${dataset}" --methods="${method}" --fit-only=yes --workers=1 ;;
    *) echo "Unknown MODE ${MODE}" >&2; exit 2 ;;
  esac
}
for attempt in 1 2 3 4 5; do
  status=0
  dispatch || status=$?
  [[ ${status} -eq 127 && ${attempt} -lt 5 ]] || break
  echo "The image did not start on $(hostname) (attempt ${attempt}); trying again" >&2
  sleep 30
done
[[ ${status} -eq 0 ]] || exit "${status}"
echo "Finished $(date --iso-8601=seconds)"
