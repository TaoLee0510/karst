#!/usr/bin/env bash
# Batch 2 of the rc freeze (results_rc/part4_rate_choice): the result folders
# of the missegregation-rate studies C1-C9, as make_copy_tasks_batch1.sh.
# Their inputs and the upstream results they read are added by
# make_copy_tasks_batch2b.sh once mapped.
#   bash make_copy_tasks_batch2.sh            (on the HPC)
set -euo pipefail
P=/share/lab_crd/taoli/Project
RC=${RC_HOME:-$P/karst/benchmark/rc_1.1.4}
R6=$P/karst/benchmark/paper/results_6ae2a24
DG=$P/karst/validation/sparse_joint/dispersion_gate
VS=$P/karst/validation/sparse_joint
PATCH=$P/karst_patch
DEVR=$P/karst_dev/results/realdata
D=$RC/results_rc/part4_rate_choice
T=$RC/Code/freeze/tasks; L=$T/lists_batch2
mkdir -p "$L"
tf=$T/copy_batch2.tsv
printf 'label\tsrc\tdst\tlist\n' > "$tf"
add() { local label=$1 src=$2 dst=$3; shift 3; printf '%s\n' "$@" > "$L/$label.txt"
        printf '%s\t%s\t%s\t%s\n' "$label" "$src" "$dst" "$L/$label.txt" >> "$tf"; }
chunks() { local label=$1 src=$2 dst=$3 size=$4 i=0
  ls "$src" | split -l "$size" - "$L/.tmp_${label}_"
  for f in "$L"/.tmp_${label}_*; do i=$((i + 1)); mv "$f" "$L/${label}_$(printf '%02d' $i).txt"
    printf '%s\t%s\t%s\t%s\n' "${label}_$(printf '%02d' $i)" "$src" "$dst" "$L/${label}_$(printf '%02d' $i).txt" >> "$tf"; done; }

## C1 fig2d_truth (paper main ABM, true rate 5e-5)
add c1_top $R6/fig2d_truth $D/fig2d_truth code engine_check.tsv jobs.tsv logs salehi_cells.tsv submit_f2t.sh tasks step1 step2/rows step2/truth
chunks c1_fits $R6/fig2d_truth/step2/fits/abm $D/fig2d_truth/step2/fits/abm 25
## C2 fig2d_misseg (varyMisseg, fixed 5e-5 against the true rate)
add c2_top $R6/fig2d_misseg $D/fig2d_misseg code dispersion.tsv jobs.tsv logs salehi_cells.tsv submit_f2m.sh tasks truth beta_fixed/rows beta_true/rows
chunks c2_fixed $R6/fig2d_misseg/beta_fixed/fits/misseg $D/fig2d_misseg/beta_fixed/fits/misseg 25
chunks c2_true $R6/fig2d_misseg/beta_true/fits/misseg $D/fig2d_misseg/beta_true/fits/misseg 25
## C3 dispersion gate
add c3_top $DG/results $D/dispersion_gate/results jobs.tsv logs stage1 submit_gate.sh step1/pick step1/tasks
for x in $(ls $DG/results/step1/v7); do add c3_v7_$x $DG/results/step1/v7 $D/dispersion_gate/results/step1/v7 $x; done
for x in $(ls $DG/results/step1/vm); do add c3_vm_$x $DG/results/step1/vm $D/dispersion_gate/results/step1/vm $x; done
## C4 forecast rate, C6 out-of-domain, C7 marginal PL (patch studies)
add c4_forecast_rate $PATCH/forecast_rate $D/forecast_rate $(ls $PATCH/forecast_rate)
add c6_outside_domain $PATCH/outside_domain $D/outside_domain $(ls $PATCH/outside_domain)
add c7_marginal_pl $PATCH/marginal_pl $D/marginal_pl $(ls $PATCH/marginal_pl)
## C5 real-data beta study (dev)
for x in $(ls $DEVR/beta); do add c5_$x $DEVR/beta $D/realdata_beta $x; done
## C8 beta from the data: noise N1-N3, v2-v7; the alfak2-era cloud-level study
add c8_noise_n1 $DEVR/noise $D/beta_from_data/noise $(cd $DEVR/noise && ls -d n1_*)
for x in $(cd $DEVR/noise && ls -d v*); do add c8_noise_$x $DEVR/noise $D/beta_from_data/noise $x; done
add c8_cloud_level $VS/cloud_level $D/beta_from_data/cloud_level $(ls $VS/cloud_level)
## C9 the alfak2-era beta-choice study
add c9_beta_choice $VS/beta_choice $D/beta_choice $(ls $VS/beta_choice)
echo "$(( $(wc -l < "$tf") - 1 )) tasks in $tf"
