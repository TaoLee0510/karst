#!/usr/bin/env bash
# Verification of batch 2, after the first comparison (REPORT_batch2.md):
#   C2's beta_true analysis again with F2T_SALEHI_TSV (the round's Salehi
#   reference table, as the as-run figure analysis/figures/fig2d_misseg.pdf
#   had it; the first rerun left it unset);
#   C5's evaluation on the records as they were when it ran (batch2_beta_asof.sh).
#   bash Code/freeze/verify/batch2_rerun2.sh       (HPC login node)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch2
source "${RC_HOME}/Code/rc_env.sh"
export RC_UPSTREAM=${RC_FROZEN}
V=${RC_RESULTS}; R6=${RC_FROZEN}/part1_paper/round_6ae2a24; A2=${V}/links/fig2d_misseg
mon() { echo "$1 $2 monitor $(rc_monitor "$2")"; }
j=$(cd ${RC_CODE}/paper && F2T_ROOT=${A2} F2T_SOURCE=misseg F2T_BOOT=1000 F2T_CORES=4 F2T_SALEHI_TSV=${R6}/analysis/fig2d_noise/fig2d_noise.tsv \
  SCRIPT_ARGS="scripts/fig2d_truth_analyse.R beta_true" bash hpc/submit.sh script 1 32G 06:00:00 | awk '/Submitted/ {print $4}' | tail -1)
mon c2_analyse_true ${j}
j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=8G --time=06:00:00 --exclude=${RC_EXCLUDE} --job-name=v2_beta_asof \
  --output=${V}/logs/%x_%j.out --export=ALL,RC_HOME=${RC_HOME} ${RC_CODE}/freeze/verify/batch2_beta_asof.sh)
mon beta_asof ${j}
