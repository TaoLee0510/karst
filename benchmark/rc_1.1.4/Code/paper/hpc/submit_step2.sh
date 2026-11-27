#!/usr/bin/env bash
# Step 2 of the cross-validation study: the paper's tests with KARST's
# cross-validation from another source (XVS: alfak_target, own_target, ...;
# XV_SOURCE in R/xv_source.R), into analysis_xv_<XVS>/.
#   XVS=own_target bash submit_step2.sh
# alfak's and the paper's Fig. 4 and 6 arms
# do not depend on KARST's scores and are copied from analysis/fig456/; the
# KARST arms (fig4 tasks 1-40, fig6 tasks 1-4) are rerun.
# /etc/profile (sbatch, lmod) does not run under set -u
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
# rc freeze: the round is the rerun round under RC_RESULTS (run A1 first); the
# library is karst_wip, the staged 1.1.4 work tree this study ran with
# (Code/env/lib/karst_wip; source Code/env/src/karst_wip).
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/rc_env.sh"
STAGE=${RC_CODE}/paper
R6=${RC_RESULTS}/part1_paper/round_6ae2a24
LOGS=${RC_RESULTS}/logs/paper; mkdir -p ${LOGS}
EX=${RC_EXCLUDE}
OLD=${R6}/analysis/fig456
XVS=${XVS:?}
NEW=${R6}/analysis_xv_${XVS}/fig456
[[ ! -e ${NEW} ]] || { echo "${NEW} exists already" >&2; exit 1; }
mkdir -p ${NEW}/fig4 ${NEW}/fig6 ${NEW}/fig4_all_inputs ${NEW}/fig6_all_inputs
cp -p ${OLD}/fig4/alfak_* ${OLD}/fig4/paper_* ${NEW}/fig4/
cp -p ${OLD}/fig6/alfak_2.Rds ${OLD}/fig6/alfak_all.Rds ${OLD}/fig6/paper.Rds ${NEW}/fig6/
cp -p ${OLD}/fig4_all_inputs/* ${NEW}/fig4_all_inputs/
cp -p ${OLD}/fig6_all_inputs/* ${NEW}/fig6_all_inputs/
echo "copied: fig4 $(ls ${NEW}/fig4 | wc -l), fig6 $(ls ${NEW}/fig6 | wc -l)"
common=(--qos=${RC_QOS} --exclude=${EX} --export=ALL,STAGE=${STAGE},BENCH_RESULTS=${R6},XV_SOURCE=${XVS},RC_HOME=${RC_HOME},RC_RESULTS=${RC_RESULTS},KARST_LIB=${RC_CODE}/env/lib/karst_wip)
monitor() { rc_monitor "$1"; }
record() { printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" "$(date --iso-8601=seconds)" >> ${RC_RESULTS}/logs/jobs.tsv; }
f4=$(sbatch --parsable "${common[@]}" --cpus-per-task=1 --mem=4G --time=02:00:00 --array=1-40 --job-name=xv_fig4_${XVS} \
  --output=${LOGS}/%x_%A_%a.out ${STAGE}/hpc/xv_job.sh fig4)
record ${f4} xv_fig4 fig4 4G 02:00:00 1-40; echo "fig4 ${f4} monitor $(monitor ${f4})"
f6=$(sbatch --parsable "${common[@]}" --cpus-per-task=1 --mem=4G --time=02:00:00 --array=1-4 --job-name=xv_fig6_${XVS} \
  --output=${LOGS}/%x_%A_%a.out ${STAGE}/hpc/xv_job.sh fig6)
record ${f6} xv_fig6 fig6 4G 02:00:00 1-4; echo "fig6 ${f6} monitor $(monitor ${f6})"
an=$(sbatch --parsable "${common[@]}" --cpus-per-task=4 --mem=48G --time=06:00:00 --job-name=xv_analyse_${XVS} \
  --output=${LOGS}/%x_%j.out ${STAGE}/hpc/xv_job.sh analyse_abm analyse_salehi fig5_wgd summary_comparison)
record ${an} xv_analyse "analyse_abm analyse_salehi fig5_wgd summary_comparison" 48G 06:00:00 ""; echo "analyse ${an} monitor $(monitor ${an})"
su=$(sbatch --parsable "${common[@]}" --cpus-per-task=1 --mem=16G --time=02:00:00 --job-name=xv_fig456_${XVS} \
  --dependency=afterok:${f4}:${f6}:${an} --output=${LOGS}/%x_%j.out ${STAGE}/hpc/xv_job.sh fig456_summary)
record ${su} xv_fig456 fig456_summary 16G 02:00:00 ""; echo "fig456_summary ${su} monitor $(monitor ${su})"
