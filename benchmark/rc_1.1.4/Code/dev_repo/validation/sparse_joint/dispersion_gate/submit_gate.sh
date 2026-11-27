#!/usr/bin/env bash
# Submit the dispersion gate's HPC jobs (PLAN.md); rc freeze version (as-run:
# results_rc/part4_rate_choice/dispersion_gate/results/submit_gate.sh).
#   bash submit_gate.sh v7 ARRAY            KARST fits of results/step1/tasks/fits_v7.tsv
#   bash submit_gate.sh vm M ARRAY          KARST fits of fits_vm_x<M>.tsv (benchmark run_task.R, no forecasts)
#   bash submit_gate.sh pick ARRAY          cell-half hold-out fits of pick_v7.tsv
# The gate writes to G = RC_RESULTS/part4_rate_choice/dispersion_gate/results;
# its task tables are read from GATE_TASKS (default the frozen
# results_rc/.../results/step1/tasks, as made by task_generators/); the v7
# simulations come from V7_RESULTS (the dev results root holding
# realdata/noise/v7: RC_UPSTREAM/part4_rate_choice/dev_results); fig2d_misseg's
# salehi_cells.tsv from RC_UPSTREAM. DEPEND=<job> waits (afterok, or DEPEND_TYPE).
# Every job gets an add-hours monitor and a line in RC_RESULTS/logs/jobs.tsv.
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/rc_env.sh"
G=${RC_RESULTS}/part4_rate_choice/dispersion_gate/results
TASKS=${GATE_TASKS:-${RC_FROZEN}/part4_rate_choice/dispersion_gate/results/step1/tasks}
V7_RESULTS=${V7_RESULTS:-${RC_UPSTREAM}/part4_rate_choice/dev_results}
M6=${RC_UPSTREAM}/part4_rate_choice/fig2d_misseg
EX=${RC_EXCLUDE}
monitor() { rc_monitor "$1"; }
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$(date --iso-8601=seconds)" >> ${RC_RESULTS}/logs/jobs.tsv; }
dep() { [[ -z "${DEPEND:-}" ]] || echo "--dependency=${DEPEND_TYPE:-afterok}:${DEPEND}"; }
mkdir -p ${G}/logs ${RC_RESULTS}/logs
common="--export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},V7_UPSTREAM=${V7_RESULTS}"
case "${1:?}" in
  v7)
    array=${2:?}
    export GATE_ARGS="${RC_CODE}/dev_repo ${V7_RESULTS} ${G}/step1/v7"
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${MEM:-8G} --time=06:00:00 --array=${array} $(dep) \
      --exclude=${EX} --job-name=gate_v7 --output=${G}/logs/%x_%A_%a.out \
      ${common},GATE_SCRIPT=validation/sparse_joint/dispersion_gate/step1_fit_v7.R,TASK_FILE=${TASKS}/fits_v7.tsv \
      ${RC_CODE}/dev_repo/validation/sparse_joint/dispersion_gate/gate_job.sh)
    record ${j} gate_v7 "${array}"; echo "v7 ${j} monitor $(monitor ${j})" ;;
  pick)
    array=${2:?}
    export GATE_ARGS="${RC_CODE}/dev_repo ${V7_RESULTS} ${G}/step1/pick"
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${MEM:-8G} --time=06:00:00 --array=${array} $(dep) \
      --exclude=${EX} --job-name=gate_pick --output=${G}/logs/%x_%A_%a.out \
      ${common},GATE_SCRIPT=validation/sparse_joint/dispersion_gate/step1_pick_v7.R,TASK_FILE=${TASKS}/pick_v7.tsv \
      ${RC_CODE}/dev_repo/validation/sparse_joint/dispersion_gate/gate_job.sh)
    record ${j} gate_pick "${array}"; echo "pick ${j} monitor $(monitor ${j})" ;;
  vm)
    # the paper benchmark's run_task.R with the round's 1.1.3 library, as fig2d_misseg's staged code ran it
    m=${2:?}; array=${3:?}; mm=$(printf '%02d' ${m})
    j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${MEM:-12G} --time=06:00:00 --array=${array} $(dep) \
      --exclude=${EX} --job-name=gate_vm_x${mm} --output=${G}/logs/%x_%A_%a.out \
      ${common},KARST_LOCK=${RC_CODE}/env/locks/karst_6ae2a24.lock,F2T_ROOT=${M6},BENCH_RESULTS=${G}/step1/vm/x${mm},ABM_CELLS=${M6}/salehi_cells.tsv,ABM_FORECAST_HORIZON=3,TASK_FILE=${TASKS}/fits_vm_x${mm}.tsv,BENCH_TASKS=none \
      ${RC_CODE}/paper/hpc/f2t_job.sh fit)
    record ${j} gate_vm_x${mm} "${array}"; echo "vm x${mm} ${j} monitor $(monitor ${j})" ;;
esac
