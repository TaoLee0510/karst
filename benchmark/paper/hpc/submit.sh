#!/usr/bin/env bash
# Submit an array of tasks on qos xxlarge, one CPU per task, no limit on the
# number of tasks running at once. Run from the repository root.
#
#   hpc/submit.sh run TASK_FILE MEM TIME [ARRAY]   tasks of a task table
#   hpc/submit.sh reproduce N MEM TIME [ARRAY]     reproduction check of the
#                                                  first N paper simulations
#   hpc/submit.sh score_salehi N MEM TIME [ARRAY]  scores of the first N
#                                                  Salehi trajectories
#   hpc/submit.sh fig4 N MEM TIME [ARRAY]          Fig. 4 bootstrap chunks
#                                                  (arm x chunk, 70 tasks)
#   hpc/submit.sh fig6 N MEM TIME [ARRAY]          Fig. 6 steady states, one
#                                                  task per arm (7)
#   SCRIPT_ARGS="scripts/x.R ..." hpc/submit.sh script N MEM TIME
#                                                  any script, N times
#   hpc/submit.sh reproduce_salehi N MEM TIME      the same check for the
#                 [ARRAY]                          paper's Salehi fits
#   hpc/submit.sh horizon TASK_FILE MEM TIME [ARRAY]  Salehi horizon test
#                                                  (HORIZON_OUT: output folder)
#   hpc/submit.sh dev TASK_FILE MEM TIME [ARRAY]   development benchmark fits
#   hpc/submit.sh far TASK_FILE MEM TIME [ARRAY]   its distance 3-4 sweep fits
#   hpc/submit.sh build MEM TIME                   install karst (one job)
#   hpc/submit.sh build_extra MEM TIME             install src/extra/ into
#                                                  lib/extra_r (one job)
#
# ARRAY defaults to every task (1-N). DEPEND=<job id>[:<job id> ...] starts
# the job only after those jobs succeeded (DEPEND_TYPE=afterany: ended). QOS=<qos> replaces xxlarge (12 h at most), e.g.
# QOS=large for up to 24 h. DRY_RUN=1 prints the sbatch call only.
# EXCLUDE=<nodes> keeps the jobs off those nodes (sbatch has no SBATCH_EXCLUDE
# variable, so that name is read here too; default: the three nodes below).
# Every submitted job is recorded in BENCH_RESULTS/jobs.tsv (job id, mode,
# task file, memory, time, array, time of submission), from which
# hpc/sacct_jobs.sh reads each task's peak memory.
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
load_hpc_environment
mode="$1"
mkdir -p "${BENCH_DIR}/logs"
qos="${QOS:-xxlarge}"
exclude="${EXCLUDE:-${SBATCH_EXCLUDE:-hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063}}"
extra=(--exclude="${exclude}")
if [[ -n "${DEPEND:-}" ]]; then
  extra+=(--dependency="${DEPEND_TYPE:-afterok}:${DEPEND}")
fi

if [[ "${mode}" == build || "${mode}" == build_extra ]]; then
  if [[ "${mode}" == build ]]; then job_name=karst_build; script=hpc/build_karst.sh
  else job_name=extra_r_build; script=hpc/build_extra_r.sh; fi
  cmd=(sbatch --parsable "${extra[@]}" --qos="${qos}" --ntasks=1 --cpus-per-task=1 --mem="$2" --time="$3"
       --job-name="${job_name}" --output="${BENCH_DIR}/logs/%x_%j.out" --export="ALL,BENCH_DIR=${BENCH_DIR}"
       "${script}")
  task_file=""
  mem="$2"
  limit="$3"
  array=""
else
  if [[ "${mode}" == run || "${mode}" == horizon || "${mode}" == dev || "${mode}" == far ]]; then
    task_file="$(realpath "$2")"
    n=$(($(wc -l <"${task_file}") - 1))
    name="bench_$(basename "${task_file}" .tsv)"
  else
    task_file=""
    n="$2"
    name="bench_${mode}"
  fi
  array="${5:-1-${n}}"
  mem="$3"
  limit="$4"
  cmd=(sbatch --parsable ${extra[@]+"${extra[@]}"} --qos="${qos}" --ntasks=1 --cpus-per-task=1 --mem="${mem}"
       --time="${limit}" --array="${array}" --job-name="${name}" --output="${BENCH_DIR}/logs/%x_%A_%a.out"
       --export="ALL,BENCH_DIR=${BENCH_DIR},MODE=${mode},TASK_FILE=${task_file}" hpc/array_task.sh)
fi
printf '%q ' "${cmd[@]}"
printf '\n'
if [[ "${DRY_RUN:-0}" != 1 ]]; then
  job="$("${cmd[@]}")"
  job="${job%%;*}"
  echo "Submitted batch job ${job}"
  registry="${BENCH_RESULTS}/jobs.tsv"
  mkdir -p "${BENCH_RESULTS}"
  [[ -s "${registry}" ]] || printf 'job_id\tmode\ttask_file\tmem\ttime\tarray\tsubmitted\n' >"${registry}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${job}" "${mode}" "${task_file}" "${mem}" "${limit}" "${array}" \
    "$(date --iso-8601=seconds)" >>"${registry}"
fi
