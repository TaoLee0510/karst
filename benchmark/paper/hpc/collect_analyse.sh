#!/usr/bin/env bash
# Tables and tests of a round, after its fits and the Salehi scores. Run as a
# Slurm job from the repository root, optionally naming the steps to run:
#
#   sbatch --qos=xxlarge --cpus-per-task=4 --mem=64G --time=06:00:00 \
#     --exclude=hpctpa3pc0013,hpctpa3pc0020,hpctpa3pc0033,hpctpa3pc0063 \
#     --output=logs/%x_%j.out --job-name=collect_analyse hpc/collect_analyse.sh [STEP ...]
#
# collect_salehi.R forks SLURM_CPUS_PER_TASK workers for the regressions of
# Fig. 3, each about 4 GB: 8 workers do not fit in 32 GB.
# Step fig5_wgd runs twice: on the analysable trajectories, and on every
# trajectory for the paper's own fits (reference of fig456_summary.R).
# Step "sacct" (on the host) reads the peak memory of every task of the job
# registries REGISTRIES (default: this round's BENCH_RESULTS/jobs.tsv and
# hpc/jobs/8327282.tsv, the round alfak's fits are linked from) into
# BENCH_RESULTS/tables/sacct.tsv; RESOURCES_ARGS goes to
# collect_resources.R (default: alfak=8327282 karst=jobs). With
# COMPARE_WITH (an earlier round's results folder), the default steps end
# with compare_rounds.R: whether karst's two-time-point fits are unchanged
# (COMPARE_PATTERN, a pattern of the designs compared, default "ends$"; "."
# compares every design of method karst, e.g. the default "anchored" fits of
# a round that should equal the earlier one).
set -euo pipefail
cd "${BENCH_DIR:-$PWD}"
source hpc/env.sh
load_hpc_environment
steps=("$@")
if [[ ${#steps[@]} -eq 0 ]]; then
  steps=(sacct collect_abm analyse_abm example_abm check_ksets collect_salehi analyse_salehi cv_score_analysis
         collect_resources summary_comparison)
  # karst's "auto" arm (scripts/auto_picks.R, compose_auto.R): its criterion
  [[ ! -s "${BENCH_RESULTS}/tables/auto_picks.tsv" ]] || steps+=(auto_criterion)
  # karst's missegregation rate picked per input (scripts/beta_picks.R,
  # compose_beta.R; ac36fd1 round): its criterion
  [[ ! -s "${BENCH_RESULTS}/tables/beta_picks.tsv" ]] || steps+=(beta_criterion)
  # the paper's Fig. 5 (one job) and, once the Fig. 4 and 6 arrays have run
  # (hpc/submit.sh fig4 / fig6), the Fig. 4-6 conclusions
  steps+=(fig5_wgd)
  [[ ! -d "${BENCH_RESULTS}/analysis/fig456/fig4" && ! -d "${BENCH_RESULTS}/analysis/fig456/fig6" ]] ||
    steps+=(fig456_summary)
  [[ -z "${COMPARE_WITH:-}" ]] || steps+=(compare_rounds)
fi
read -r -a registries <<<"${REGISTRIES:-${BENCH_RESULTS}/jobs.tsv hpc/jobs/8327282.tsv}"
read -r -a resources_args <<<"${RESOURCES_ARGS:-alfak=8327282 karst=jobs}"
echo "Results ${BENCH_RESULTS}  karst ${KARST_COMMIT}  $(date --iso-8601=seconds)"
for step in "${steps[@]}"; do
  echo "== ${step} $(date --iso-8601=seconds)"
  case "${step}" in
    sacct) hpc/sacct_jobs.sh "${BENCH_RESULTS}/tables/sacct.tsv" "${registries[@]}" ;;
    collect_resources) container_rscript scripts/collect_resources.R "${resources_args[@]}" ;;
    compare_rounds) container_rscript scripts/compare_rounds.R "${COMPARE_WITH}" "${COMPARE_PATTERN:-ends\$}" ;;
    fig5_wgd)
      FIG456_EXCLUDE=1 container_rscript scripts/fig5_wgd.R
      FIG456_EXCLUDE=0 container_rscript scripts/fig5_wgd.R ;;
    *) container_rscript "scripts/${step}.R" ;;
  esac
done
echo "Finished $(date --iso-8601=seconds)"
