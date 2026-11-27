#!/usr/bin/env bash
# Verification of batch 2: C5's evaluation (beta_eval.R, then beta_summary.R)
# on the records as they were when it ran. The frozen beta/evaluation.rds was
# written on 2026-11-20 at 05:06; at 07:27 the decomp study's extension
# (beta_salehi_ext, fit_ext.tsv) added 310 Salehi fits at x128 and x256 to the
# same beta/fits/salehi folders, which a rerun on the final folders reads too.
# The mirror here holds the files of dev_eval_files.txt under beta/ that are
# not newer than the frozen evaluation.rds (links to results_rc), plus the
# linked inputs. Submitted by batch2_rerun2.sh as a Slurm job.
#   sbatch ... --export=ALL,RC_HOME=... Code/freeze/verify/batch2_beta_asof.sh
source /etc/profile >/dev/null 2>&1 || true
set -euo pipefail
RC_HOME=${RC_HOME:?}
export RC_RESULTS=${RC_HOME}/verify/batch2
source "${RC_HOME}/Code/rc_env.sh"
export RC_UPSTREAM=${RC_FROZEN}
V=${RC_RESULTS}; FD=${RC_FROZEN}/part4_rate_choice/dev_results/realdata
DA=${V}/dev_eval_asof; rm -rf ${DA}; mkdir -p ${DA}/realdata ${V}/tasks
(source ${RC_CODE}/dev_repo/dev/hpc/env.sh; DEV_RESULTS=${DA}; dev_link_inputs)
(cd ${FD} && find beta -type f ! -newer beta/evaluation.rds | LC_ALL=C sort) > ${V}/tasks/beta_asof_files.txt
grep '^beta/' ${RC_CODE}/freeze/verify/dev_eval_files.txt | LC_ALL=C sort | LC_ALL=C comm -12 - ${V}/tasks/beta_asof_files.txt > ${V}/tasks/beta_asof_mirror.txt
while IFS= read -r f; do mkdir -p "${DA}/realdata/$(dirname "${f}")"; ln -sfn "${FD}/${f}" "${DA}/realdata/${f}"; done < ${V}/tasks/beta_asof_mirror.txt
echo "as-of mirror: $(wc -l < ${V}/tasks/beta_asof_mirror.txt) links ($(grep -c '^beta/' ${RC_CODE}/freeze/verify/dev_eval_files.txt) in the full eval mirror)"
ev() { DEV_RESULTS=${DA} DEV_LIB=karst_cd342db JOB_NAME=v2_asof_$(basename $1 .R) DEPEND=${2:-} bash ${RC_CODE}/dev_repo/dev/hpc/submit_args.sh $1 32G 08:00:00 | tail -1; }
out=$(ev dev/realdata/beta/beta_eval.R); echo "${out}"
ev dev/realdata/beta/beta_summary.R $(awk '{print $2}' <<< "${out}")
