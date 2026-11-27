#!/usr/bin/env bash
# Verification of batch 2, the dev studies' downstream and evaluation part of
# batch2_run.sh, as a Slurm job (building its mirrors took too long for the
# login node, which ended batch2_run.sh there). Otherwise as batch2_run.sh:
# sampled reruns of
# every study's steps from Code/ and data/ into RC_HOME/verify/batch2, each on
# the frozen upstream results (RC_UPSTREAM = results_rc). A step that reads
# frozen outputs of its own study gets them through a mirror (folders real,
# files symbolic links to results_rc; only in verify/, never in results_rc);
# its own outputs are always computed fresh. Compared afterwards by
# batch2_compare.sh. Submitted by batch2_run.sh at its end. Each of its roots
# (dev_down, dev_screen, dev_sfit, dev_eval) starts empty, so every output
# there has one writer.
#   sbatch ... --export=ALL,RC_HOME=... Code/freeze/verify/batch2_run_dev.sh
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=${RC_HOME:?}
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
## C5 real-data beta, C8a noise (dev studies) --------------------------------------------------------
FD=${F}/${P4}/dev_results/realdata; D=${RC_CODE}/dev_repo/dev/realdata
dev() {   # dev ROOT LIB SCRIPT TASK_FILE ARRAY [MEM]
  DEV_RESULTS=$1 DEV_LIB=$2 JOB_NAME=v2_$(basename $3 .R) bash ${RC_CODE}/dev_repo/dev/hpc/submit.sh $3 $4 ${6:-16G} 12:00:00 $5 | tail -1; }
links() { (source ${RC_CODE}/dev_repo/dev/hpc/env.sh; DEV_RESULTS=$1; dev_link_inputs); }
rm -rf ${V}/dev_down ${V}/dev_screen ${V}/dev_sfit ${V}/dev_eval
# downstream steps on the frozen upstream (only their inputs mirrored)
DD=${V}/dev_down; mkdir -p ${DD}/realdata/beta ${DD}/realdata/noise/v3 ${DD}/realdata/noise/v7; links ${DD}
mirror ${FD}/beta/salehi ${DD}/realdata/beta/salehi
mirror ${FD}/noise/v3/sims ${DD}/realdata/noise/v3/sims
mirror ${FD}/noise/v4 ${DD}/realdata/noise/v4
mirror ${FD}/noise/v7/sims ${DD}/realdata/noise/v7/sims
fd=$(row ${D}/noise/v4/tasks/fits_dev.tsv dataset_id 1); tag=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) h[$i] = i; next} NR == 2 {print $h["kind"] "_" $h["set"] "_" $h["variant"]; exit}' ${D}/noise/v4/tasks/fits_dev.tsv)
mt=$(row ${D}/noise/v4/tasks/fits_dev.tsv method 1)
rm -f ${DD}/realdata/noise/v4/fits/${mt}/${tag}/${fd}.rds   # refitted
tag5=$(row ${D}/noise/v5/tasks/step12.tsv tag 1); tag6=$(row ${D}/noise/v6/tasks/stats.tsv tag 1)
dev ${DD} karst_cd342db dev/realdata/beta/beta_horizon.R ${D}/beta/tasks/beta_horizon.tsv 1
dev ${DD} karst_cd342db dev/realdata/beta/decomp/decomp_horizon.R ${D}/beta/decomp/tasks/decomp.tsv 1
dev ${DD} karst_cd342db dev/realdata/beta/disp/disp_forecast.R ${D}/beta/disp/tasks/disp.tsv 1
dev ${DD} karst_cd342db dev/realdata/noise/v3/v3_task.R ${D}/noise/v3/tasks/analysis.tsv 1
DS=${V}/dev_screen; mkdir -p ${DS}/realdata/noise/v3; links ${DS}; mirror ${FD}/noise/v3/sims ${DS}/realdata/noise/v3/sims
dev ${DS} karst_cd342db dev/realdata/noise/v4/screen_task.R ${D}/noise/v4/tasks/screen.tsv 1
DSF=${V}/dev_sfit; mkdir -p ${DSF}/realdata/noise/v4/screen; links ${DSF}; mirror ${FD}/noise/v4/screen/salehi ${DSF}/realdata/noise/v4/screen/salehi
dev ${DSF} karst_v4 dev/realdata/noise/v4/salehi_fit.R ${D}/noise/v4/tasks/salehi_fit.tsv 1
dev ${DD} karst_v4 dev/realdata/noise/v4/fit_task.R ${D}/noise/v4/tasks/fits_dev.tsv 1 8G
dev ${DD} karst_v4 dev/realdata/noise/v5/step12_task.R ${D}/noise/v5/tasks/step12.tsv 1 6G
dev ${DD} karst_v4 dev/realdata/noise/v6/stats_task.R ${D}/noise/v6/tasks/stats.tsv 1 8G
dev ${DD} karst_6ae2a24 dev/realdata/noise/v7/stats_task.R ${D}/noise/v7/tasks/stats.tsv 1
# evaluations on the frozen results (a full mirror, evaluation outputs removed)
# The evaluations ran locally on a partial copy of the dev results: the mirror
# holds exactly the files that copy held (dev_eval_files.txt, written
# from the local dev/results/realdata), each a link to the identical frozen file.
DE=${V}/dev_eval; mkdir -p ${DE}/realdata; links ${DE}
while IFS= read -r f; do mkdir -p "${DE}/realdata/$(dirname "${f}")"; ln -sfn "${FD}/${f}" "${DE}/realdata/${f}"; done < ${RC_CODE}/freeze/verify/dev_eval_files.txt
echo "eval mirror: $(find ${DE}/realdata -type l | wc -l) links"
ev() { DEV_RESULTS=${DE} DEV_LIB=$1 JOB_NAME=v2_$(basename $2 .R) DEPEND=${3:-} bash ${RC_CODE}/dev_repo/dev/hpc/submit_args.sh $2 32G 08:00:00 | tail -1; }
# one after the other (afterany), as they ran: disp_eval reads decomp_eval's
# evaluation.rds, beta_summary beta_eval's
j=
for e in karst_cd342db:beta/beta_eval.R karst_cd342db:beta/beta_summary.R karst_cd342db:beta/decomp/decomp_eval.R \
         karst_cd342db:beta/disp/disp_eval.R karst_cd342db:noise/n23_eval.R karst_cd342db:noise/v2_eval.R \
         karst_cd342db:noise/v3/v3_eval.R karst_cd342db:noise/v4/v4_eval.R karst_v4:noise/v5/step12_eval.R \
         karst_v4:noise/v6/calib_eval.R karst_6ae2a24:noise/v7/calib_eval.R; do
  out=$(ev ${e%%:*} dev/realdata/${e#*:} ${j}); echo "${out}"; j=$(awk '{print $2}' <<< "${out}")
done
echo "submitted; v4 fit ${mt}/${tag}/${fd}, v5 ${tag5}, v6 ${tag6}"
