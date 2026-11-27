#!/usr/bin/env bash
# Verification of batch 2, step 2: compare every file the reruns of
# batch2_run.sh wrote (links of the mirrors left out) with the same path in
# results_rc (compare.R --pairs, one job in the image).
#   bash Code/freeze/verify/batch2_compare.sh      (HPC login node, after batch2_run.sh's jobs)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch2
source "${RC_HOME}/Code/rc_env.sh"
V=${RC_RESULTS}; F=${RC_FROZEN}/part4_rate_choice; P=${V}/compare/pairs.tsv; mkdir -p "${V}/compare"
printf 'label\trerun\tfrozen\n' > "${P}"
pair() { [[ -e "$2" ]] && printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "${P}" || echo "missing rerun: $2" >&2; }
pair c1_setup ${V}/fresh/fig2d_truth ${F}/fig2d_truth
pair c1_fit_rows ${V}/part4_rate_choice/fig2d_truth ${F}/fig2d_truth
pair c1_analysis ${V}/links/fig2d_truth ${F}/fig2d_truth
pair c2_setup ${V}/fresh/fig2d_misseg ${F}/fig2d_misseg
pair c2_fit ${V}/part4_rate_choice/fig2d_misseg ${F}/fig2d_misseg
pair c2_rows ${V}/links/fig2d_misseg_rows ${F}/fig2d_misseg
pair c2_analysis ${V}/links/fig2d_misseg ${F}/fig2d_misseg
pair c3_gate ${V}/part4_rate_choice/dispersion_gate ${F}/dispersion_gate
pair c3_tasks ${V}/links/gate_tasks ${F}/dispersion_gate/results
pair c4_forecast_rate ${V}/part4_rate_choice/forecast_rate ${F}/forecast_rate
pair c6_outside_domain ${V}/part4_rate_choice/outside_domain ${F}/outside_domain
pair c7_marginal_pl ${V}/part4_rate_choice/marginal_pl ${F}/marginal_pl
for r in dev_fit dev_down dev_screen dev_sfit dev_eval; do pair dev_${r} ${V}/${r}/realdata ${F}/dev_results/realdata; done
pair dev_eval_asof ${V}/dev_eval_asof/realdata ${F}/dev_results/realdata   # batch2_beta_asof.sh
echo "$(( $(wc -l < "${P}") - 1 )) pairs"
cd "${RC_CODE}/paper"
SCRIPT_ARGS="../freeze/verify/compare.R --pairs ${P} ${V}/compare" bash hpc/submit.sh script 1 32G 04:00:00 | tail -1
