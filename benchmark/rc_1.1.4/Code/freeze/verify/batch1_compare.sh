#!/usr/bin/env bash
# Verification of batch 1, step 3: pair every rerun output under
# RC_HOME/verify/batch1 with its frozen counterpart in results_rc and compare
# (compare.R --pairs, one job in the image); the reproduction check and the
# own-benchmark subset are summarised by their own scripts.
#   bash Code/freeze/verify/batch1_compare.sh      (HPC login node, after steps 1 and 2)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch1
source "${RC_HOME}/Code/rc_env.sh"
V=${RC_RESULTS}; F=${RC_FROZEN}; L=${V}/links; P=${V}/compare/pairs.tsv; mkdir -p "${V}/compare"
printf 'label\trerun\tfrozen\n' > "${P}"
pair() { [[ -e "$2" ]] && printf '%s\t%s\t%s\n' "$1" "$2" "$3" >> "${P}" || echo "missing rerun: $2" >&2; }
# step 1: sampled fits (each method folder of a rerun input)
for d in $(cd "${V}" && find part1_paper/round_6ae2a24/fits part1_paper/salehi_rate_grid -mindepth 4 -maxdepth 6 -type d \
             \( -name karst -o -name alfak -o -name karst_pairwise -o -name karst_state_space \) 2>/dev/null); do
  pair "fit_$(echo "${d}" | tr '/' '_')" "${V}/${d}" "${F}/${d}"
done
for f in $(cd "${V}" && find part2_own -path '*/fits/*.rds' -o -path '*/objects/*.rds' 2>/dev/null); do
  pair "own_$(basename "${f}" .rds)" "${V}/${f}" "${F}/${f}"
done
pair kester_inputs "${V}/part3_kester/inputs" "${F}/part3_kester/inputs"
pair kester_manifest "${V}/part3_kester/fold_manifest.tsv" "${F}/part3_kester/fold_manifest.tsv"
pair kester_tasks "${V}/part3_kester/tasks.tsv" "${F}/part3_kester/tasks.tsv"
for d in $(cd "${V}" && find part3_kester/fits -mindepth 3 -maxdepth 3 -type d 2>/dev/null); do
  pair "kester_$(echo "${d}" | tr '/' '_')" "${V}/${d}" "${F}/${d}"
done
# step 2: analyses on the frozen intermediate files
pair a1_analysis_abm "${L}/round/analysis/abm" "${F}/part1_paper/round_6ae2a24/analysis/abm"
pair a1_analysis_salehi "${L}/round/analysis/salehi" "${F}/part1_paper/round_6ae2a24/analysis/salehi"
pair a3_xv_alfak_target "${L}/round/tables/xv_alfak_target.Rds" "${F}/part1_paper/round_6ae2a24/tables/xv_alfak_target.Rds"
for d in $(cd "${L}/round" && find align_karst/fits -mindepth 4 -maxdepth 4 -type d 2>/dev/null); do
  pair "a4_$(echo "${d}" | tr '/' '_')" "${L}/round/${d}" "${F}/part1_paper/round_6ae2a24/${d}"
done
pair a2_compare "${L}/salehi_rate_grid/analysis" "${F}/part1_paper/salehi_rate_grid/analysis"
pair a2_paper_grid "${L}/paper_grid" "${F}/part1_paper/paper_grid"
pair b2_summary "${L}/own_grid/summary" "${F}/part2_own/rate_grid/summary"
pair d_kester_scores "${L}/kester/scores.tsv" "${F}/part3_kester/scores.tsv"
pair d_kester_summary "${L}/kester/summary" "${F}/part3_kester/summary"
echo "$(( $(wc -l < "${P}") - 1 )) pairs"
cd "${RC_CODE}/paper"
SCRIPT_ARGS="../freeze/verify/compare.R --pairs ${P} ${V}/compare" bash hpc/submit.sh script 1 32G 04:00:00 | tail -1
SCRIPT_ARGS="../freeze/verify/own_subset_check.R ${L}/own_grid_subset ${F}/part2_own/rate_grid ${V}/compare/own_subset.tsv" \
  bash hpc/submit.sh script 1 8G 01:00:00 | tail -1
SCRIPT_ARGS="../freeze/verify/rc_check_summary.R ${V}/part0_freeze/rc_check2 ${V}/compare/rc_check.tsv" \
  bash hpc/submit.sh script 1 4G 00:30:00 | tail -1
