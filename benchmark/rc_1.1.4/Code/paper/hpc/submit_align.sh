#!/usr/bin/env bash
# Version B (ALIGN=karst) of round 6ae2a24, after the pilot: the full round's
# inputs, the Salehi scores, the collecting and the analyses, into
# results_6ae2a24/align_karst/. alfak's and KARST's Fig. 4 and 6 arms are
# rerun (their landscapes are on E); the paper's own arm is copied.
#   bash submit_align.sh tasks      the inputs table and the align_task array
#   bash submit_align.sh downstream Salehi scores, collect, analyses, Fig. 4-6
#                                   (the ABM collect and analysis run in a job of their own)
#                                   (DEPEND=<job id> to wait for the array)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
# rc freeze: the round is the rerun round under RC_RESULTS (run A1 first); the
# library is karst_wip, the staged 1.1.4 work tree this study ran with
# (Code/env/lib/karst_wip; source Code/env/src/karst_wip).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"
STAGE=${RC_CODE}/paper
R6=${RC_RESULTS}/part1_paper/round_6ae2a24
LOGS=${RC_RESULTS}/logs/paper; mkdir -p ${LOGS}
B=${R6}/align_karst
EX=${RC_EXCLUDE}
COMMON="ALL,STAGE=${STAGE},BENCH_RESULTS=${R6},RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},KARST_LIB=${RC_CODE}/env/lib/karst_wip,ALIGN=karst,BENCH_TASKS=${RC_CODE}/paper/tasks/6ae2a24"
monitor() { rc_monitor "$1"; }
record() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$(date --iso-8601=seconds)" >> ${RC_RESULTS}/logs/jobs.tsv; }
dep() { [[ -z "${1:-}" ]] || echo "--dependency=afterok:$1"; }
case "${1:?}" in
  tasks)
    mk=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=4G --time=00:30:00 --exclude=${EX} --job-name=align_make \
      --output=${LOGS}/%x_%j.out --export=${COMMON} ${STAGE}/hpc/xv_job.sh align_make_tasks)
    record ${mk} align_make make-tasks 4G 00:30:00 ""; echo "make ${mk} monitor $(monitor ${mk})"
    n=${N_INPUTS:?}
    a=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${ALIGN_MEM:-16G} --time=12:00:00 --array=1-${n} \
      --dependency=afterok:${mk} --exclude=${EX} --job-name=align --output=${LOGS}/%x_%A_%a.out \
      --export=${COMMON},ALIGN_TASKS=${B}/tasks/inputs.tsv ${STAGE}/hpc/xv_job.sh align_task)
    record ${a} align "align_task 1-${n}" ${ALIGN_MEM:-16G} 12:00:00 1-${n}; echo "align ${a} monitor $(monitor ${a})" ;;
  downstream)
    d=$(dep "${DEPEND:-}")
    sc=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=16G --time=12:00:00 --array=1-${N_TRAJ:?} ${d} \
      --exclude=${EX} --job-name=align_score --output=${LOGS}/%x_%A_%a.out --export=${COMMON} \
      ${STAGE}/hpc/xv_job.sh score_salehi)
    record ${sc} align_score "score_salehi 1-${N_TRAJ}" 16G 12:00:00 1-${N_TRAJ}; echo "score ${sc} monitor $(monitor ${sc})"
    co=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=4 --mem=96G --time=12:00:00 --dependency=afterok:${sc} \
      --exclude=${EX} --job-name=align_collect --output=${LOGS}/%x_%j.out --export=${COMMON},SLURM_CPUS_PER_TASK=4 \
      ${STAGE}/hpc/xv_job.sh check_ksets_salehi collect_salehi analyse_salehi fig5_wgd)
    record ${co} align_collect "check_ksets_salehi collect_salehi analyse_salehi fig5_wgd" 96G 12:00:00 ""
    echo "collect ${co} monitor $(monitor ${co})"
    f4=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=16G --time=06:00:00 --array=1-60 --dependency=afterok:${co} \
      --exclude=${EX} --job-name=align_fig4 --output=${LOGS}/%x_%A_%a.out --export=${COMMON} ${STAGE}/hpc/xv_job.sh fig4)
    record ${f4} align_fig4 fig4 16G 06:00:00 1-60; echo "fig4 ${f4} monitor $(monitor ${f4})"
    f6=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=32G --time=12:00:00 --array=1-6 --dependency=afterok:${co} \
      --exclude=${EX} --job-name=align_fig6 --output=${LOGS}/%x_%A_%a.out --export=${COMMON} ${STAGE}/hpc/xv_job.sh fig6)
    record ${f6} align_fig6 fig6 32G 12:00:00 1-6; echo "fig6 ${f6} monitor $(monitor ${f6})"
    # the paper's own arm (its fits are the paper's, not aligned): copied
    mkdir -p ${B}/analysis/fig456/fig4 ${B}/analysis/fig456/fig6 ${B}/analysis/fig456/fig4_all_inputs ${B}/analysis/fig456/fig6_all_inputs
    cp -p ${R6}/analysis/fig456/fig4/paper_* ${B}/analysis/fig456/fig4/
    cp -p ${R6}/analysis/fig456/fig6/paper.Rds ${B}/analysis/fig456/fig6/
    cp -p ${R6}/analysis/fig456/fig4_all_inputs/* ${B}/analysis/fig456/fig4_all_inputs/
    cp -p ${R6}/analysis/fig456/fig6_all_inputs/* ${B}/analysis/fig456/fig6_all_inputs/
    su=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=16G --time=02:00:00 --dependency=afterok:${f4}:${f6} \
      --exclude=${EX} --job-name=align_fig456 --output=${LOGS}/%x_%j.out --export=${COMMON} \
      ${STAGE}/hpc/xv_job.sh fig456_summary)
    record ${su} align_fig456 fig456_summary 16G 02:00:00 ""; echo "fig456 ${su} monitor $(monitor ${su})" ;;
esac
