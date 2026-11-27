#!/usr/bin/env bash
# Batch 1 of the rc freeze (results_rc/README.md): builds the copy task table
# Code/freeze/tasks/copy_batch1.tsv and one list file per task. Each task copies
# the entries of its list (paths relative to SRC) from SRC to DST with rsync -aL
# (symlinks are followed, so every copy is a real file). Sources are read only.
#   bash make_copy_tasks_batch1.sh            (on the HPC, from any directory)
set -euo pipefail
P=/share/lab_crd/taoli/Project
RC=${RC_HOME:-$P/karst/benchmark/rc_1.1.4}
R6=$P/karst/benchmark/paper/results_6ae2a24
PAPER=$P/karst/benchmark/paper
OWN=$P/karst/benchmark/results/sparse_two_timepoint_6ae2a24
RCW=$P/karst_rc_1.1.4
KES=$P/Kester/results
PD=/share/lab_crd/M010_ALFAK_2023/ALFA-K/data
SIF=/share/lab_crd/taoli/Docker/cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif
RR=$RC/results_rc
T=$RC/Code/freeze/tasks; L=$T/lists_batch1
mkdir -p "$L"
tf=$T/copy_batch1.tsv
printf 'label\tsrc\tdst\tlist\n' > "$tf"
add() {   # add LABEL SRC DST ENTRY...
  local label=$1 src=$2 dst=$3; shift 3
  printf '%s\n' "$@" > "$L/$label.txt"
  printf '%s\t%s\t%s\t%s\n' "$label" "$src" "$dst" "$L/$label.txt" >> "$tf"
}
chunks() {   # chunks LABEL SRC DST SIZE : one task per SIZE entries of SRC
  local label=$1 src=$2 dst=$3 size=$4 i=0
  ls "$src" | split -l "$size" - "$L/.tmp_${label}_"
  for f in "$L"/.tmp_${label}_*; do
    i=$((i + 1)); mv "$f" "$L/${label}_$(printf '%02d' $i).txt"
    printf '%s\t%s\t%s\t%s\n' "${label}_$(printf '%02d' $i)" "$src" "$dst" "$L/${label}_$(printf '%02d' $i).txt" >> "$tf"
  done
}

## E. part 0: freeze and reproduction check ----------------------------------------------------
add p0_results $RCW/results $RR/part0_freeze rc_check2 rc_check_superseded_refit_from_fit_data own_pilot pilotcmp
add p0_work $RCW $RR/logs/rc_work code jobs.tsv logs tasks

## A1 round 6ae2a24 (paper benchmark, every test) ---------------------------------------------
D1=$RR/part1_paper/round_6ae2a24
chunks a1_fits_abm $R6/fits/abm $D1/fits/abm 25
chunks a1_fits_misseg $R6/fits/misseg $D1/fits/misseg 25
chunks a1_fits_salehi $R6/fits/salehi $D1/fits/salehi 25
add a1_round $R6 $D1 analysis figures scores tables jobs.tsv fig6_dominant_source.Rds
add a1_paper_logs $PAPER $RR/part1_paper logs tasks
## A3 CV-target variants
add a3_xv $R6 $D1 analysis_xv_alfak_target analysis_xv_own_target figures_xv_alfak_target figures_xv_own_target xv_target_code
## A4 version B
chunks a4_fits_abm $R6/align_karst/fits/abm $D1/align_karst/fits/abm 25
chunks a4_fits_misseg $R6/align_karst/fits/misseg $D1/align_karst/fits/misseg 25
chunks a4_fits_salehi $R6/align_karst/fits/salehi $D1/align_karst/fits/salehi 25
add a4_align $R6/align_karst $D1/align_karst align_checks_steady.Rds align_tasks.Rds analysis fig6_dominant_source.Rds figures scores tables tasks
## A2 Salehi rate grid x02-x64 (x01 fits are the round's fits/salehi)
D2=$RR/part1_paper/salehi_rate_grid
add a2_grid $R6/beta_grid $D2 code jobs.tsv logs oom_tasks.txt submit_bg.sh tasks
add a2_x01 $R6/beta_grid/x01 $D2/x01 analysis figures scores tables
for m in 02 04 08 16 32 64; do
  add a2_x${m} $R6/beta_grid/x${m} $D2/x${m} analysis figures scores tables
  chunks a2_x${m}_fits $R6/beta_grid/x${m}/fits/salehi $D2/x${m}/fits/salehi 50
done

## B1 own round 6ae2a24 (results; datasets and nodes go to data/) -----------------------------
D3=$RR/part2_own/round_6ae2a24
add b1_round $OWN $D3 fits truth tasks nonpositive_fitness metrics_heldout.csv own_coverage.json own_summary.json \
  jobs.tsv round_jobs.tsv rerun_602.sh rerun_oom.sh rerun_oom2.sh score_job.sh score_split.sh submit_nonpositive.sh \
  submit_round.sh submit_summaries.sh sweeps/far/fits sweeps/far/objects sweeps/far/metrics_far_heldout.csv
## B2 own rate grid x02-x64
D4=$RR/part2_own/rate_grid
add b2_tables $RCW/results/own_grid $D4 metrics_heldout_grid.csv direction_heldout_grid.csv summary
for m in 02 04 08 16 32 64; do add b2_x${m} $RCW/results/own_grid $D4 x${m}; done

## D Kester ---------------------------------------------------------------------------------------
add d_kester $KES/forecast_50e5889 $RR/part3_kester $(ls $KES/forecast_50e5889)

## data/ -------------------------------------------------------------------------------------------
add data_paper $PD $RC/data/paper raw/ABM raw/ABM_varyMisseg raw/salehi/README.md raw/salehi/arm_loci.Rds \
  raw/salehi/metadata.csv processed/salehi/alfak_inputs processed/salehi/lineages.Rds
add data_reference $PD $RC/data/reference/paper processed/ABM processed/ABM_vary_misseg processed/ABM_summaries.Rds \
  processed/ABM_vary_misseg_summaries.Rds processed/salehi/alfak_outputs processed/salehi/alfak_outputs_proc \
  processed/salehi/arm_level_cn.Rds processed/salehi/chrom_level_cn.Rds processed/salehi/landscape_analysis_bootstrap.Rds \
  processed/salehi/novel_kary_predictions.Rds processed/salehi/steadyStatePredictions.Rds
add data_own $OWN $RC/data/own_benchmark datasets nodes dataset_index_heldout.csv
add data_own_far $OWN/sweeps/far $RC/data/own_benchmark/far nodes
add data_kester $KES $RC/data/kester/results cohort_backbone/APKS1_karst_counts.tsv cohort_backbone/APKS2_karst_counts.tsv \
  cohort_backbone/APKS3_karst_counts.tsv cohort_backbone_a3_ap1/A3_karst_counts.tsv cohort_backbone_a3_ap1/AP1_karst_counts.tsv

## Code/env: the image and the karst sources -------------------------------------------------------
add env_sif $(dirname $SIF) $RC/Code/env $(basename $SIF)
add env_karst $PAPER $RC/Code/env src/karst_50e5889.tar.gz src/karst_6ae2a24.tar.gz hpc/locks/karst_50e5889.lock hpc/karst.lock

echo "$(( $(wc -l < "$tf") - 1 )) tasks in $tf"
