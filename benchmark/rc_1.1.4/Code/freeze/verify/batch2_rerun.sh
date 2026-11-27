#!/usr/bin/env bash
# Verification of batch 2: the steps of the first run of batch2_run.sh
# (2026-11-26 18:04) that failed or ran twice, rerun with the corrected
# scripts (REPORT_batch2.md, "First run"):
#   C1 fig2d_truth_analyse.R steps 1 and 2 (Code/paper fix: no pmis column on
#   the truth simulations); C4, C6, C7 analyses (their arguments had been
#   given in the wrong order); the dev studies' downstream and evaluation steps
#   (the login-node run and batch2_run_dev.sh had both submitted them, two
#   writers per folder; batch2_run_dev.sh now starts its folders empty, runs
#   salehi_fit.R on the mirrored screen output and the evaluations in order).
#   bash Code/freeze/verify/batch2_rerun.sh        (HPC login node)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch2
source "${RC_HOME}/Code/rc_env.sh"
export RC_UPSTREAM=${RC_FROZEN}
V=${RC_RESULTS}; F=${RC_FROZEN}; P4=part4_rate_choice
mon() { echo "$1 $2 monitor $(rc_monitor "$2")"; }
A1=${V}/links/fig2d_truth; rm -rf ${A1}/analysis
for s in 1 2; do j=$(cd ${RC_CODE}/paper && F2T_ROOT=${A1} F2T_BOOT=1000 F2T_CORES=4 SCRIPT_ARGS="scripts/fig2d_truth_analyse.R ${s}" \
  bash hpc/submit.sh script 1 32G 06:00:00 | awk '/Submitted/ {print $4}' | tail -1); mon c1_analyse_${s} ${j}; done
K6=${RC_CODE}/env/lib/karst_6ae2a24; T6=${RC_CODE}/patch/outside_domain/tasks; T7=${RC_CODE}/patch/marginal_pl/tasks
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_fr_analyse patch/forecast_rate/scripts/analyse.R \
  "${F}/${P4}/forecast_rate/results/out ${V}/${P4}/forecast_rate/results/final" 2>&1 | tail -1
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_od_step1 patch/outside_domain/scripts/step1_analyse.R \
  "${F}/${P4}/outside_domain/results/sim_out ${T6}/sim_tasks.tsv ${V}/${P4}/outside_domain/results/step1" 2>&1 | tail -1
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_mp_step01 patch/marginal_pl/scripts/step01_analyse.R \
  "${F}/${P4}/marginal_pl/results/out ${T7}/sim_tasks.tsv ${V}/${P4}/marginal_pl/results/step01" 2>&1 | tail -1
j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=8G --time=12:00:00 --exclude=${RC_EXCLUDE} --job-name=v2_dev_driver \
  --output=${V}/logs/%x_%j.out --export=ALL,RC_HOME=${RC_HOME} ${RC_CODE}/freeze/verify/batch2_run_dev.sh)
mon dev_driver ${j}
