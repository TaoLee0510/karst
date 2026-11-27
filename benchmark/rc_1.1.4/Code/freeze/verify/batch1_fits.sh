#!/usr/bin/env bash
# Verification of batch 1, step 1: sampled reruns of fits from Code/ and data/
# into RC_HOME/verify/batch1 (never results_rc). Compared afterwards with
# compare.R (batch1_compare.sh).
#   bash Code/freeze/verify/batch1_fits.sh          (HPC login node, from RC_HOME)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch1
source "${RC_HOME}/Code/rc_env.sh"
V=${RC_RESULTS}; T=${V}/tasks; mkdir -p "${T}" "${V}/logs"
rows() {  # rows SRC OUT ROW... : header + the given data rows, the task column renumbered 1..n when present
  local src=$1 out=$2; shift 2
  awk -F'\t' -v OFS='\t' -v want="$*" 'BEGIN {n = split(want, w, " "); for (i = 1; i <= n; i++) keep[w[i] + 1] = i}
    NR == 1 {for (i = 1; i <= NF; i++) if ($i == "task") tc = i; print; next}
    (NR in keep) {if (tc) $tc = keep[NR]; print}' "${src}" > "${out}"
}
P=${RC_CODE}/paper/tasks
# A1 KARST: two ABM, one varyMisseg, two Salehi (ends and all) tasks of the round
rows ${P}/6ae2a24/tasks_abm_mem4G.tsv ${T}/a1_abm.tsv 1 2
rows ${P}/6ae2a24/tasks_misseg_mem4G.tsv ${T}/a1_misseg.tsv 1
rows ${P}/6ae2a24/tasks_salehi_mem8G.tsv ${T}/a1_salehi.tsv 1
rows ${P}/6ae2a24/tasks_salehi_mem16G.tsv ${T}/a1_salehi_all.tsv 1
cd "${RC_CODE}/paper"
for t in a1_abm a1_misseg a1_salehi a1_salehi_all; do
  n=$(( $(wc -l < ${T}/${t}.tsv) - 1 ))
  echo "${t}: $(bash hpc/submit.sh run ${T}/${t}.tsv 16G 12:00:00 1-${n} | tail -1)"
done
# ALFA-K task tables of the round (the round's ALFA-K fits are the alfak2-era round's)
SCRIPT_ARGS="scripts/make_tasks.R ${V}/part1_paper/tasks_alfak alfak" bash hpc/submit.sh script 1 4G 00:30:00 | tail -1
# A2 rate grid: x08 first task (KARST), x64 first alfak task
first_alfak=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "method") m = i; if ($i == "task") t = i}; next} $m == "alfak" {print $t; exit}' ${P}/salehi_rate_grid/x64/tasks_salehi.tsv)
bash "${RC_CODE}/paper/hpc/submit_rate_grid.sh" fit 8 1
bash "${RC_CODE}/paper/hpc/submit_rate_grid.sh" fit 64 "${first_alfak}"
# B1 own round: KARST default and independent, ALFA-K minobs 10, a far-sweep fit
export STT_OUTPUT_DIR=${V}/part2_own/round_6ae2a24 STT_DATA_DIR=${RC_DATA}/own_benchmark
O=${RC_CODE}/own/tasks
bal=$(awk -F'\t' '$3 == "karst_sparse_joint_balance" {print NR - 1; exit}' ${O}/karst_main.tsv)
ind=$(awk -F'\t' '$3 == "karst_sparse_joint" {print NR - 1; exit}' ${O}/karst_main.tsv)
a10=$(awk -F'\t' '$3 == "alfak_minobs10" {print NR - 1; exit}' ${O}/alfak_main.tsv)
bash hpc/submit.sh dev ${O}/karst_main.tsv 16G 12:00:00 ${bal},${ind} | tail -1
bash hpc/submit.sh dev ${O}/alfak_main.tsv 16G 12:00:00 ${a10} | tail -1
bash hpc/submit.sh far ${O}/karst_far.tsv 16G 12:00:00 1 | tail -1
unset STT_OUTPUT_DIR STT_DATA_DIR
# B2 own multiples: x2 two KARST fits, x64 one, ALFA-K x8 one
R=${RC_CODE}/rc/tasks
bash "${RC_CODE}/hpc/submit_rc.sh" v_own_k2 rc/own_task.R ${R}/own_karst_sparse_joint_balance_x2.tsv ${V}/part2_own/rate_grid 1,2
bash "${RC_CODE}/hpc/submit_rc.sh" v_own_k64 rc/own_task.R ${R}/own_karst_sparse_joint_balance_x64.tsv ${V}/part2_own/rate_grid 1
bash "${RC_CODE}/hpc/submit_rc.sh" v_own_a8 rc/own_task.R ${R}/own_alfak_minobs10_x8.tsv ${V}/part2_own/rate_grid 1
# D Kester: setup, then six fits (KARST and ALFA-K at x1 and x8, a GBM-like fold)
K=${V}/part3_kester
s=$(bash "${RC_CODE}/hpc/submit_rc.sh" v_kester_setup kester/setup.R ${R}/one.tsv ${K} | awk '{print $2}')
echo "kester setup ${s}"
cp "${RC_CODE}/kester/tasks/tasks.tsv" ${T}/kester_tasks.tsv
ids=$(awk -F'\t' 'NR > 1 && (($1 == "APKS1_P10-22_to_P25_all" && ($2 == 1 || $2 == 8) && ($3 == "karst" || $4 == "10")) || ($1 == "A3_cells_1" && $2 == 1 && ($3 == "karst" || $4 == "10"))) {printf "%s%s", sep, $5; sep = ","}' ${T}/kester_tasks.tsv)
echo "kester tasks ${ids}"
DEPEND=${s} bash "${RC_CODE}/hpc/submit_rc.sh" v_kester_fit kester/kester_task.R ${T}/kester_tasks.tsv ${K} "${ids}"
# Part 0: the whole reproduction check (60 fits) against the frozen fits
bash "${RC_CODE}/hpc/submit_rc.sh" v_rc_check rc/rc_check.R ${R}/rc_check.tsv ${V}/part0_freeze/rc_check2
