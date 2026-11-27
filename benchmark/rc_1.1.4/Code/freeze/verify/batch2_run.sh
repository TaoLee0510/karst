#!/usr/bin/env bash
# Verification of batch 2 (results_rc/part4_rate_choice): sampled reruns of
# every study's steps from Code/ and data/ into RC_HOME/verify/batch2, each on
# the frozen upstream results (RC_UPSTREAM = results_rc). A step that reads
# frozen outputs of its own study gets them through a mirror (folders real,
# files symbolic links to results_rc; only in verify/, never in results_rc);
# its own outputs are always computed fresh. Compared afterwards by
# batch2_compare.sh. The dev studies' downstream and evaluation steps are in
# batch2_run_dev.sh, which this script submits at its end as a Slurm job.
#   bash Code/freeze/verify/batch2_run.sh          (HPC login node)
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
export RC_RESULTS=${RC_HOME}/verify/batch2
source "${RC_HOME}/Code/rc_env.sh"
export RC_UPSTREAM=${RC_FROZEN}
V=${RC_RESULTS}; F=${RC_FROZEN}; P4=part4_rate_choice; mkdir -p "${V}/logs" "${V}/tasks"
L6=${RC_CODE}/env/locks/karst_6ae2a24.lock
mirror() { local src=$1 dst=$2; mkdir -p "$(dirname "${dst}")"; cp -rsn "${src}" "$(dirname "${dst}")/"; }   # dst = dirname/basename(src)
mon() { echo "$1 $2 monitor $(rc_monitor "$2")"; }
row() { awk -F'\t' -v c="$2" -v r="${3:-1}" 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == c) k = i; next} NR == r + 1 {print $k}' "$1"; }
rowsel() { awk -F'\t' -v r="$2" 'NR == 1 || NR == r + 1' "$1"; }   # header + data row r
f2t() {   # f2t NAME STEP ARRAY EXPORTS... : one f2t_job.sh array (paper code), 1.1.3 library as run
  local name=$1 step=$2 array=$3; shift 3; local ex; ex=$(IFS=,; echo "$*")
  sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=${MEM:-16G} --time=12:00:00 --array=${array} --exclude=${RC_EXCLUDE} \
    --job-name=v2_${name} --output=${V}/logs/%x_%A_%a.out \
    --export=ALL,RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},RC_UPSTREAM=${RC_UPSTREAM},KARST_LOCK=${L6},${ex} \
    ${RC_CODE}/paper/hpc/f2t_job.sh ${step}
}
R6=${F}/part1_paper/round_6ae2a24

## C1 fig2d_truth -----------------------------------------------------------------------------
T1=${V}/${P4}/fig2d_truth; FT1=${F}/${P4}/fig2d_truth
mon c1_setup "$(f2t c1_setup 'cells tasks engine' 1 F2T_ROOT=${V}/fresh/fig2d_truth BENCH_RESULTS=${R6})"
sim1=$(row ${FT1}/tasks/sims_pilot.tsv sim_id 1)
mkdir -p ${T1}; cp -p ${FT1}/salehi_cells.tsv ${T1}/; mirror ${FT1}/tasks ${T1}/tasks
mkdir -p ${T1}/step2/fits/abm; mirror ${FT1}/step2/fits/abm/${sim1} ${T1}/step2/fits/abm/${sim1}
d1=$(row ${FT1}/tasks/tasks_pilot.tsv design 1); m1=$(row ${FT1}/tasks/tasks_pilot.tsv method 1)
rm -rf ${T1}/step2/fits/abm/${sim1}/${d1}/${m1}   # task 1's fit is refitted
mon c1_fit "$(f2t c1_fit fit 1 F2T_ROOT=${T1} BENCH_RESULTS=${T1}/step2 ABM_CELLS=${T1}/salehi_cells.tsv ABM_FORECAST_HORIZON=3 TASK_FILE=${FT1}/tasks/tasks_pilot.tsv BENCH_TASKS=none)"
mon c1_rows2 "$(MEM=16G f2t c1_rows2 rows 1 F2T_ROOT=${T1} BENCH_RESULTS=${FT1}/step2 ABM_CELLS=${T1}/salehi_cells.tsv F2T_STEP=2 F2T_SIMS=pilot)"
mon c1_rows1 "$(MEM=16G f2t c1_rows1 rows 1 F2T_ROOT=${T1} BENCH_RESULTS=${R6} ABM_CELLS= F2T_STEP=1 F2T_SIMS=pilot)"
A1=${V}/links/fig2d_truth; mkdir -p ${A1}/step1 ${A1}/step2; cp -p ${FT1}/salehi_cells.tsv ${A1}/; mirror ${FT1}/tasks ${A1}/tasks
mirror ${FT1}/step1/rows ${A1}/step1/rows; mirror ${FT1}/step2/rows ${A1}/step2/rows
for s in 1 2; do j=$(cd ${RC_CODE}/paper && F2T_ROOT=${A1} F2T_BOOT=1000 F2T_CORES=4 SCRIPT_ARGS="scripts/fig2d_truth_analyse.R ${s}" \
  bash hpc/submit.sh script 1 32G 06:00:00 | awk '/Submitted/ {print $4}' | tail -1); mon c1_analyse_${s} ${j}; done

## C2 fig2d_misseg ------------------------------------------------------------------------------
T2=${V}/${P4}/fig2d_misseg; FT2=${F}/${P4}/fig2d_misseg
mon c2_setup "$(f2t c2_setup 'cells tasks dispersion' 1 F2T_ROOT=${V}/fresh/fig2d_misseg BENCH_RESULTS=${R6} F2T_SOURCE=misseg)"
sim2=$(row ${FT2}/tasks/sims_beta_fixed_pilot.tsv sim_id 1)
mkdir -p ${T2}/beta_fixed/fits/misseg; cp -p ${FT2}/salehi_cells.tsv ${T2}/; mirror ${FT2}/tasks ${T2}/tasks
mirror ${FT2}/beta_fixed/fits/misseg/${sim2} ${T2}/beta_fixed/fits/misseg/${sim2}
d2=$(row ${FT2}/tasks/tasks_beta_fixed_pilot.tsv design 1); m2=$(row ${FT2}/tasks/tasks_beta_fixed_pilot.tsv method 1)
rm -rf ${T2}/beta_fixed/fits/misseg/${sim2}/${d2}/${m2}
mon c2_fit "$(MEM=12G f2t c2_fit fit 1 F2T_ROOT=${T2} BENCH_RESULTS=${T2}/beta_fixed ABM_CELLS=${T2}/salehi_cells.tsv ABM_FORECAST_HORIZON=3 MISSEG_FORECAST=TRUE TASK_FILE=${FT2}/tasks/tasks_beta_fixed_pilot.tsv BENCH_TASKS=none)"
mkdir -p ${V}/links/fig2d_misseg_rows/beta_fixed; cp -p ${FT2}/salehi_cells.tsv ${V}/links/fig2d_misseg_rows/; mirror ${FT2}/tasks ${V}/links/fig2d_misseg_rows/tasks
mon c2_rows "$(f2t c2_rows rows 1 F2T_ROOT=${V}/links/fig2d_misseg_rows BENCH_RESULTS=${FT2}/beta_fixed ABM_CELLS=${FT2}/salehi_cells.tsv ABM_FORECAST_HORIZON=3 F2T_SOURCE=misseg F2T_ARM=beta_fixed F2T_SIMS=pilot)"
A2=${V}/links/fig2d_misseg; mkdir -p ${A2}; cp -p ${FT2}/salehi_cells.tsv ${FT2}/dispersion.tsv ${A2}/; mirror ${FT2}/tasks ${A2}/tasks
for a in beta_fixed beta_true; do mkdir -p ${A2}/${a}; mirror ${FT2}/${a}/rows ${A2}/${a}/rows; done
j=$(cd ${RC_CODE}/paper && F2T_ROOT=${A2} F2T_SOURCE=misseg F2T_BOOT=1000 F2T_CORES=4 F2T_SALEHI_TSV=${R6}/analysis/fig2d_noise/fig2d_noise.tsv SCRIPT_ARGS="scripts/fig2d_truth_analyse.R beta_fixed" bash hpc/submit.sh script 1 32G 06:00:00 | awk '/Submitted/ {print $4}' | tail -1); mon c2_analyse_fixed ${j}
j=$(cd ${RC_CODE}/paper && DEPEND=${j} DEPEND_TYPE=afterany F2T_ROOT=${A2} F2T_SOURCE=misseg F2T_BOOT=1000 F2T_CORES=4 F2T_SALEHI_TSV=${R6}/analysis/fig2d_noise/fig2d_noise.tsv SCRIPT_ARGS="scripts/fig2d_truth_analyse.R beta_true" bash hpc/submit.sh script 1 32G 06:00:00 | awk '/Submitted/ {print $4}' | tail -1); mon c2_analyse_true ${j}

## C3 dispersion gate ---------------------------------------------------------------------------
G=${V}/${P4}/dispersion_gate/results; FG=${F}/${P4}/dispersion_gate/results; S0=${RC_DATA}/dispersion_gate_step0
bash ${RC_CODE}/hpc/submit_rc_args.sh v2_gate_step0 dev_repo/validation/sparse_joint/dispersion_gate/step0_dispersion.R \
  "--paper-main=${RC_DATA}/paper/raw/ABM --paper-vm=${RC_DATA}/paper/raw/ABM_varyMisseg --v7-e=${S0}/v7_e --v7-c2=${S0}/v7_c2 --own=${S0}/own --out=${G}/step0 --cores=4" 2>&1 | tail -1
GT=${V}/links/gate_tasks; mkdir -p ${GT}/step1/tasks; mirror ${FG}/step0 ${GT}/step0
j=$(GATE_RESULTS=${GT} bash ${RC_CODE}/hpc/submit_rc_args.sh v2_gate_mk_v7_tasks \
  dev_repo/validation/sparse_joint/dispersion_gate/task_generators/mk_v7_tasks.R "-" 2>&1 | tail -1 | awk '{print $2}')
echo "mk_v7 ${j}"
DEPEND=${j} GATE_RESULTS=${GT} bash ${RC_CODE}/hpc/submit_rc_args.sh v2_gate_mk_pick_tasks \
  dev_repo/validation/sparse_joint/dispersion_gate/task_generators/mk_pick_tasks.R "-" 2>&1 | tail -1
GATE_RESULTS=${GT} bash ${RC_CODE}/hpc/submit_rc_args.sh v2_gate_mk_vm dev_repo/validation/sparse_joint/dispersion_gate/task_generators/mk_vm_tasks.R \
  "${F}/${P4}/fig2d_misseg/tasks" 2>&1 | tail -1
GATE_TASKS=${FG}/step1/tasks bash ${RC_CODE}/dev_repo/validation/sparse_joint/dispersion_gate/submit_gate.sh v7 1
GATE_TASKS=${FG}/step1/tasks bash ${RC_CODE}/dev_repo/validation/sparse_joint/dispersion_gate/submit_gate.sh pick 1
GATE_TASKS=${FG}/step1/tasks bash ${RC_CODE}/dev_repo/validation/sparse_joint/dispersion_gate/submit_gate.sh vm 2 1
CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_gate_analyse dev_repo/validation/sparse_joint/dispersion_gate/step1_analyse.R \
  "--v7=${FG}/step1/v7 --pick=${FG}/step1/pick --vm=${FG}/step1/vm --vm-x1=${FG}/step1/vm_x1 --vm-raw=${RC_DATA}/paper/raw/ABM_varyMisseg --cells=${F}/${P4}/fig2d_misseg/salehi_cells.tsv --gate-results=${FG} --out=${G}/step1/analysis" 2>&1 | tail -1

## C4 forecast rate, C6 out-of-domain, C7 marginal PL ---------------------------------------------
K6=${RC_CODE}/env/lib/karst_6ae2a24
O4=${V}/${P4}/forecast_rate/out; T4=${RC_CODE}/patch/forecast_rate/tasks
abm4=$(awk -F'\t' 'NR > 1 && $1 == "abm" {print NR - 1; exit}' ${T4}/sims.tsv); vm4=$(awk -F'\t' 'NR > 1 && $1 == "misseg" {print NR - 1; exit}' ${T4}/sims.tsv)
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc.sh v2_fr_sim patch/forecast_rate/scripts/sim_task.R ${T4}/sims.tsv ${O4} ${abm4},${vm4}
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc.sh v2_fr_alfak patch/forecast_rate/scripts/alfak_task.R ${T4}/sims.tsv ${O4} ${abm4},${vm4}
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc.sh v2_fr_salehi patch/forecast_rate/scripts/salehi_task.R ${T4}/salehi.tsv ${O4} 1
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_fr_analyse patch/forecast_rate/scripts/analyse.R \
  "${F}/${P4}/forecast_rate/results/out ${V}/${P4}/forecast_rate/results/final" 2>&1 | tail -1
O6=${V}/${P4}/outside_domain/out; T6=${RC_CODE}/patch/outside_domain/tasks
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc.sh v2_od_sim patch/outside_domain/scripts/sim_task.R ${T6}/sim_tasks.tsv ${O6} 1
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_od_step1 patch/outside_domain/scripts/step1_analyse.R \
  "${F}/${P4}/outside_domain/results/sim_out ${T6}/sim_tasks.tsv ${V}/${P4}/outside_domain/results/step1" 2>&1 | tail -1
KARST_LIB=${K6} CPUS=3 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_od_step0 patch/outside_domain/scripts/step0_engine_check.R "3" 2>&1 | tail -1
O7=${V}/${P4}/marginal_pl/out; T7=${RC_CODE}/patch/marginal_pl/tasks
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc.sh v2_mp_own patch/marginal_pl/scripts/own_task.R ${T7}/own_dev.tsv ${O7} 1
rowsel ${T7}/sim_tasks.tsv 1 > ${V}/tasks/mp_paper_1.tsv; awk -F'\t' 'NR > 1 && $1 == "misseg" {print; exit}' ${T7}/sim_tasks.tsv >> ${V}/tasks/mp_paper_1.tsv
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc_args.sh v2_mp_paper patch/marginal_pl/scripts/paper_task.R "${V}/tasks/mp_paper_1.tsv ${O7}" 2>&1 | tail -1
KARST_LIB=${K6} bash ${RC_CODE}/hpc/submit_rc_args.sh v2_mp_salehi patch/marginal_pl/scripts/salehi_pl.R "-" 2>&1 | tail -1
KARST_LIB=${K6} CPUS=4 bash ${RC_CODE}/hpc/submit_rc_args.sh v2_mp_step01 patch/marginal_pl/scripts/step01_analyse.R \
  "${F}/${P4}/marginal_pl/results/out ${T7}/sim_tasks.tsv ${V}/${P4}/marginal_pl/results/step01" 2>&1 | tail -1

## C5 real-data beta, C8a noise (dev studies) --------------------------------------------------------
FD=${F}/${P4}/dev_results/realdata; D=${RC_CODE}/dev_repo/dev/realdata
dev() {   # dev ROOT LIB SCRIPT TASK_FILE ARRAY [MEM]
  DEV_RESULTS=$1 DEV_LIB=$2 JOB_NAME=v2_$(basename $3 .R) bash ${RC_CODE}/dev_repo/dev/hpc/submit.sh $3 $4 ${6:-16G} 12:00:00 $5 | tail -1; }
links() { (source ${RC_CODE}/dev_repo/dev/hpc/env.sh; DEV_RESULTS=$1; dev_link_inputs); }
# fits from the inputs only (fresh root)
DF=${V}/dev_fit; mkdir -p ${DF}; links ${DF}
dev ${DF} karst_cd342db dev/realdata/beta/beta_sim_fit.R ${D}/beta/tasks/beta_sim.tsv 1
dev ${DF} karst_cd342db dev/realdata/beta/beta_salehi.R ${D}/beta/tasks/beta_salehi.tsv 1
dev ${DF} karst_cd342db dev/realdata/beta/beta_numbat.R ${D}/beta/tasks/beta_numbat.tsv 1
nb=$(awk -F'\t' 'NR > 1 && $0 ~ /numbat/ {print NR - 1; exit}' ${D}/noise/tasks/n1.tsv)
dev ${DF} karst_cd342db dev/realdata/noise/n1_measure.R ${D}/noise/tasks/n1.tsv ${nb}
dev ${DF} karst_cd342db dev/realdata/noise/v2_task.R ${D}/noise/tasks/v2_sim.tsv 1
dev ${DF} karst_cd342db dev/realdata/noise/v3/sim_eset.R ${D}/noise/v3/tasks/truth.tsv 1
dev ${DF} karst_6ae2a24 dev/realdata/noise/v7/sim_eset.R ${D}/noise/v7/tasks/sims_e.tsv 1
dev ${DF} karst_6ae2a24 dev/realdata/noise/v7/sim_c2.R ${D}/noise/v7/tasks/sims_c2.tsv 1
# the dev studies' downstream and evaluation steps: batch2_run_dev.sh, as a Slurm job
# (its mirrors take too long for the login node)
j=$(sbatch --parsable --qos=${RC_QOS} --cpus-per-task=1 --mem=8G --time=12:00:00 --exclude=${RC_EXCLUDE} --job-name=v2_dev_driver \
  --output=${V}/logs/%x_%j.out --export=ALL,RC_HOME=${RC_HOME} ${RC_CODE}/freeze/verify/batch2_run_dev.sh)
mon dev_driver ${j}
echo "submitted; c1 sim ${sim1} (${d1} ${m1}), c2 sim ${sim2} (${d2} ${m2})"
