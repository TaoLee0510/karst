# KARST with the tests of the ALFA-K paper

This folder of the karst repository evaluates the default model of KARST
(the sparse joint count model, `karst::fit_sparse_joint()`) with the tests of

> Beck, Li, Andor. ALFA-K: local adaptive mapping of karyotype fitness
> landscapes. *Nat Commun* (2025). doi:10.1038/s41467-025-67750-0

and compares every score with alfak. The reference alfak is **alfak
dev2@f93b2d2**: package alfakR, branch `alfakR_dev2` of TaoLee0510/alfakR at
commit f93b2d2, the version PANcanKFLs uses (the code calls it as
`alfakR::`). The paper's own fits are **alfak main@4ad82b0**: they date from
June 2025, and 4ad82b0 (2025-05-14) is the last commit of Richard-Beck/alfakR
main before them (inferred from the dates; the paper does not name a commit).
They are used only to check that this code reproduces the paper's scores and
as a reference line ("paper").

The code is not part of the KARST package (`.Rbuildignore`). Up to version
0.2.0 KARST was called alfak2 and this benchmark was a repository of its own
(alfak2-paper-benchmark); the reports of that time are in `docs/alfak2_era/`.

## Design

Two time points are the minimum input of KARST. Each time window is fitted with
its first and last time point by both methods, and with every time point by alfak
(the paper's design) and, from KARST's multi-timepoint version on, by
KARST:

| Fit | Input |
|---|---|
| KARST, two time points | first and last time point of the window |
| KARST, all time points | every time point of the window (multi-timepoint KARST) |
| alfak, two time points | the same two time points |
| alfak, all time points | every time point of the window (the paper's design) |
| KARST pairwise, all time points | every time point, `fit_sparse_joint(series = "pairwise")` (c7f809c on) |
| KARST state space, all time points | every time point, `series = "state_space"` (c7f809c on) |

"KARST, all time points" is KARST's default, `series = "anchored"`. The
pairwise and state-space methods are fitted to the designs with more than
two time points only (with two time points all three give the same fit; methods
`karst_pairwise`, `karst_state_space`).

Two arms of the alfak2 era are not run any more: "auto" (one of the three
"all" fits picked per input by a rule; it failed its criterion in round
05a8855 and KARST has no `series = "auto"`) and the missegregation rate
picked per input (`beta = "auto"`, removed after round ac36fd1). Their
scripts are kept (`scripts/auto_*.R`, `compose_auto.R`, `beta_*.R`,
`compose_beta.R`, `make_beta_tasks.R`, `prune_beta_candidates.R`); the steps
that use them run only when their tables exist (`hpc/collect_analyse.sh`).

Each method is compared with the other on the same time points (two time points;
all time points). `scripts/summary_comparison.R` puts both groups in one figure
(`figures/summary_comparison`): accuracy on paired inputs, and failures,
fit time and peak memory (`tables/resources.Rds`, from `run.Rds` and
Slurm's accounting) over every input. The CV score is reported as its
median and share > 0 only: each method's score is measured against its own
full-data estimates, so the two are not comparable.

* **ABM data** (paper Supplementary Note 3): the paper's windows are the last
  2, 4 and 8 passages before day 120. Reduced to their first and last
  passage they are 1, 3 and 7 passages apart (about 6, 17 and 41 days). The
  2-passage window is the paper's `ntp = 2`.
* **Salehi data** (paper Fig. 2–3): each of the 277 trajectories reduced to
  its first and last passage, 1 to 9 passages apart (5 to 135 days). The
  trajectories include every run of consecutive passages, so the first +
  middle passage pairs of long trajectories are among these pairs.
* **Missegregation rate** (paper Supplementary Fig. 3): the 125 simulations
  at five rates, same windows, minobs 10.

How the scores change with the time between the two time points is the
stability of each method.

## Tests

| Paper | Here |
|---|---|
| Supp. Fig. 2a: ρ, r, R² overall and for frequent, distance-1 and distance-2 karyotypes, by λ and minobs | three fits × three windows |
| Supp. Fig. 2b, 2f: frequent-karyotype R² vs overall R², Sankey | same |
| Supp. Fig. 2c–e: example simulation (w_1p6_m_5e-05_rep_74, λ = 1.6) | same simulation, every fit (`scripts/example_abm.R`) |
| Supp. Note 3.2: initial frequencies (Pearson, RMSE, late karyotypes) | all fits |
| Supp. Fig. 3: robustness to the missegregation rate (Pearson r) | three fits × three windows |
| Supp. Fig. 4a–c: cross-validation score as a diagnostic | same |
| Supp. Fig. 5a: forecast angle θ by sign of the CV score | same forecaster, only the landscape differs |
| Supp. Fig. 5b, Fig. 1g: fraction of forecasts beating the no-evolution baseline | same |
| Fig. 2a–b: CV score vs training length | by the interval between time points |
| Fig. 2c: angle metric, cluster bootstrap test | paper's selection rule, and paired by gap |
| Fig. 2d: predicted vs observed Wasserstein distance | same |
| Fig. 2e–f: fraction beating the static baseline | same |
| Fig. 2g–h, unrelated lineages, theoretical null | do not depend on the method; reused |
| Fig. 3e–f: logistic regression of novel karyotype emergence | f from each fit |

| Fig. 4c, d, f: fitness effects of CNAs, mixed models over 200 bootstrap draws of non-overlapping lineages | each fit's landscapes, paper's selection (`scripts/fig4_cna_effects.R`) |
| Fig. 5b–c: WGD, landscape shape and fitness effects (SA906) | same (`scripts/fig5_wgd.R`) |
| Fig. 6: steady states under missegregation rates 4e-5 to 5.1e-3 | same (`scripts/fig6_steady_states.R`) |

Fig. 4–6 are not tests of accuracy: they ask whether each method's fits lead
to the biological conclusions the paper drew from its own fits
(`scripts/fig456_summary.R`, `analysis/fig456/conclusions.tsv`; reference:
the paper's stored results, and its own fits through the same code). They
were added for the multi-timepoint round 05a8855. Not repeated, because they
do not depend on the fits or only illustrate: Fig. 4a (lineage sampling),
4b (example heatmap), 4e (passage overlap, counts only), 5a (counts only),
the analytic scenarios of Fig. 6. Fig. 4 needs glmmTMB, which the image
lacks: version 1.1.11 is installed into `lib/extra_r`
(`hpc/build_extra_r.sh`).

Two ways of reporting: the paper's rule (failed fits dropped, test
instances selected by each method's own CV score) and a paired comparison
on the inputs where both methods have a fit, without selection.

## Settings

* alfak algorithm settings as PANcanKFLs: `nboot = 500`,
  `correct_efflux = TRUE`, `passage_times = NULL` (times are column names ×
  `dt`), minobs 5, 10 and 20.
* Data settings as the paper: ABM `n0 = 2e5`, `nb = 2e7`, `pm = 5e-5`;
  missegregation simulations `pm` = the simulated rate; Salehi `n0 = 1e5`,
  `nb = 1e7`, `pm = 5e-5`.
* KARST: default model, the same `pm`, `n0`, `nb`, `nboot` and
  `correct_efflux`, outputs in alfak's format (`karst::as_alfak()`) for
  each minobs. On the Salehi data `min_cn = 0`, because most inputs contain
  karyotypes called with zero copies of a chromosome, which alfak accepts.
* Forecasts: `alfakR::predict_evo()` in ABM mode for both methods, with the
  paper's settings (S03 for the ABM data, `02_run_alfak.R` for the Salehi
  data), five replicates with seeds 1–5.

Both methods run in the PANcanKFLs image
`cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif`; KARST is installed into
a side library with that image's R and compilers (`hpc/build_karst.sh`),
from the commit in `hpc/karst.lock`.

## Inputs that cannot be analysed

KARST models the karyotypes next to the major ones of its time points. When
every cell of one time point lies two or more steps from every major karyotype
(with the first and last passage: 24 Salehi trajectories, all from the SA609
family), the time points are not linked at the karyotype level and KARST cannot
fit them. This is decided from the data before any fit, with the tested
KARST's own rules (`R/coverage.R`), for each design on the time points it fits
(first and last; every time point). A fit is analysed only when its own input is
analysable, for every method, the paper's line included;
`excluded_inputs.tsv` shows what the other fits give on the excluded inputs.
The ABM and missegregation inputs are all analysable, with two time points and
with every time point.

With more than two time points (8512694), KARST takes its majors from the
first and the last time point only and builds its domain around them; the first
and the last time point must have a cell inside it. A time point in between with no
cell on a major says nothing under the model and is left out of the fit
(`diagnostics$dropped_samples`). An input on every time point is therefore
analysable exactly when its first and last time points are: on the Salehi data the
same 24 trajectories are excluded on every passage as on first and last. 53
of the 187 analysable trajectories of more than two passages lose 103
passages in between, 10 of them all of their passages in between (their fit
is the fit of first and last). Most are from the SA000, SA001 and
SA609-derived families, where a clone rises and falls between the first and
the last passage. (With 21aedca, which stopped instead, 61 of these 205
trajectories could not be fitted.)

## Rounds

Each KARST commit has its own results folder, `results_<commit>`, and its
own task tables, `tasks/<commit>/`, on the HPC in
`/share/lab_crd/taoli/Project/karst/benchmark/paper` (a clone of the karst
repository; results, tasks, logs and libraries are not under version
control). A round refits KARST only. alfak's fits are linked from the
benchmark of the alfak2 era (`hpc/link_alfak.sh`, from
`/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results`, round
8327282, where the folders are named after the package, `alfakR`); alfak,
its settings and the inputs do not change. That folder is bound read only
into the image (`hpc/env.sh`), so that the links resolve there.
`hpc/submit.sh` records every job in `results_<commit>/jobs.tsv`; the job
registries of the alfak2-era rounds go in `hpc/jobs/` (not tracked; the copy
of 8327282's points to its task tables by absolute path).
`scripts/compare_rounds.R` compares KARST's fits with an earlier round's
(`COMPARE_WITH=<results folder> hpc/collect_analyse.sh`); for a round of
the alfak2 era, whose KARST folders are named `alfak2`, `alfak2_pairwise`
and `alfak2_state_space`, set `COMPARE_OLD_METHOD` to that name.

| KARST | Results | Report |
|---|---|---|
| b9a275e (1.0.0, the rename; stopped with every fit but one done, superseded by e7e32f8) | `results_b9a275e/` | |
| e7e32f8 (1.0.0 with the fixes of the code review of 2026-11-17; all KARST fits redone) | `results_e7e32f8/` | `docs/report_e7e32f8.md` |

The same round fits KARST with `initial_state = "free"` to the Salehi data
(method `karst_free`, both designs) into `results_<commit>_free/`, for the
horizon test (`../../validation/sparse_joint/salehi_horizon`); the round's
tables do not include it.

A round whose commit should leave the fits unchanged can refit a check subset
only (`tasks/<commit>/`). If `scripts/compare_rounds.R`
(`COMPARE_METHOD=karst | karst_pairwise | karst_state_space`) finds the
subset identical to the earlier round, link every other fit from there
(`SKIP_EXISTING=1 SOURCE_NAME=<method> hpc/link_alfak.sh results_<earlier>
<method>`). Link only after every check task has ended: `run_task.R`
creates a fit folder late, so a running task's folder would be linked and
then written through. Queued jobs read `hpc/karst.lock` when they start; to
prepare a round while another round's jobs are queued, put its lock in
`hpc/locks/` and export `KARST_LOCK` for its submissions (`hpc/env.sh`).

## Data

Read in place, never copied into this repository:

* `PAPER_DATA` (default `/share/lab_crd/M010_ALFAK_2023/ALFA-K/data`):
  `raw/ABM` (the manuscript simulations), `raw/ABM_varyMisseg`,
  `processed/salehi/`, and the paper's fits and summaries in `processed/`.
* Results go to `BENCH_RESULTS` (default `<this folder>/results_<commit>`,
  not under version control).

## Running

On the HPC, from this folder of the karst clone. Before every array, report
its tasks and resources; give each array an add-hours monitor
(`~/submit_hourly_add_hours.sh`).

1. **Stage the commit.** `git archive --format=tar.gz <commit>` of the
   karst repository to `src/karst_<commit>.tar.gz`, and the commit and the
   archive's SHA-256 in `hpc/karst.lock`.
2. **Build.** `hpc/submit.sh build 8G 02:00:00` installs karst into
   `lib/karst_<commit>`; glmmTMB for Fig. 4 goes to `lib/extra_r`
   (`hpc/submit.sh build_extra 8G 02:00:00`, or a copy of the alfak2 era's).
3. **Link alfak.** `hpc/link_alfak.sh "${ALFAK_REFERENCE}"`.
4. **Task tables** (inside the image, `container_rscript` of `hpc/env.sh`):

   ```bash
   Rscript scripts/make_tasks.R tasks/<commit> karst karst_pairwise karst_state_space
   Rscript scripts/split_tasks.R tasks/<commit>/tasks_abm.tsv tasks/<commit>/tasks_misseg.tsv \
     tasks/<commit>/tasks_salehi.tsv
   ```

5. **Fits**, one array per memory tier (`tasks/<commit>/tiers.tsv`):
   `hpc/submit.sh run <tier file> <mem> <time>`.
6. **After the fits**, in this order: `hpc/submit.sh score_salehi 277 8G
   02:00:00`; `hpc/collect_analyse.sh` (sbatch, 4 CPUs, 64 GB, 6 h; with
   `REGISTRIES="hpc/jobs/8327282.tsv ${BENCH_RESULTS}/jobs.tsv"`); the Fig. 4
   and Fig. 6 arrays (`hpc/submit.sh fig4 70 8G 04:00:00` and, for the
   paper's fits on every trajectory, `FIG456_EXCLUDE=0 ... fig4 70 8G
   04:00:00 61-70`; `fig6 7 8G 06:00:00` and `FIG456_EXCLUDE=0 ... fig6 7 8G
   06:00:00 7`); then `hpc/collect_analyse.sh fig456_summary`; the Fig. 2d
   noise check (`SCRIPT_ARGS="scripts/fig2d_noise.R --stage=chunk"
   hpc/submit.sh script 28 8G 03:00:00`, then `--stage=analyse`); and
   `COMPARE_WITH=... hpc/collect_analyse.sh compare_rounds`.

Fig. 4 and Fig. 6 read `tables/salehi_scores.Rds`, so they run after
`collect_analyse`. Slurm requeues a task after NODE_FAIL by itself (check
`sacct -D`); do not resubmit it.

KARST's tiers: 4 GB and 6 h for two time points of the ABM and
missegregation data; 8 GB and 12 h for the Salehi data and the anchored fits
of all time points; 16 GB for pairwise and 32 GB for state space (12 h).

Each array task fits one method to one input, then forecasts and scores
(`scripts/run_task.R`); finished steps are skipped when a task is rerun.
alfak's memory grows with the number of time points it fits and falls with
minobs (from under 1 GB to over 32 GB); its fits are linked, not rerun.

The same scripts run two other checks in this image with the round's
library: the Salehi horizon test (`hpc/submit.sh horizon TASK_FILE ...`,
`../../validation/sparse_joint/salehi_horizon`) and the development
benchmark of simulated tumours (`hpc/submit.sh dev TASK_FILE ...`,
`../sparse_two_timepoint`).

## Files

* `R/config.R` paths, settings and the ABM windows
* `R/abm.R` simulations, true fitness, inputs
* `R/salehi.R` trajectories, inputs, forecasts
* `R/coverage.R` whether KARST can fit an input, decided before fitting
* `R/fit.R` fitting with alfak or KARST
* `R/labels.R` names and colours of the fits compared
* `R/forecast.R` ABM forecasts
* `R/metrics_abm.R` scores of S04 (`abm_preds_summary.Rds`)
* `R/metrics_salehi.R` scores of 03, references of 04, regression data of 05
* `scripts/check_reproduction.R`, `scripts/check_reproduction_salehi.R` the
  scores of the paper's fits, recomputed
* `scripts/make_tasks.R`, `scripts/split_tasks.R`, `scripts/run_task.R` task
  tables, memory tiers and one task
* `scripts/score_salehi.R` Salehi scores of one trajectory
* `scripts/collect_abm.R`, `scripts/collect_salehi.R` tables of all fits
* `scripts/analyse_abm.R`, `scripts/analyse_salehi.R` the paper's tests for
  every fit, and the paired comparison
* `scripts/cv_score_analysis.R` what the CV score of each method measures
* `scripts/collect_resources.R` fit time and peak memory of every input
* `scripts/summary_comparison.R` KARST against alfak on two time points and on
  all time points, in one figure
* `scripts/compare_rounds.R` KARST's fits of two rounds, side by side
* `scripts/initial_state_summary.R`, `scripts/initial_state_vs_truth.R`
  KARST's choice of `initial_state = "auto"` on two time points (alfak2
  round 827c4b5): share switched to "free", and the switched fits against
  the reference round's, each against the truth
* `scripts/fig2d_noise.R` the paper's Fig. 2d on the Salehi passages, raw and
  with the passages' sampling noise held fixed (uses `scripts/fig2d_check.R`)
* `scripts/example_abm.R` the paper's Supp. Fig. 2c–e for every fit
* `scripts/auto_picks.R`, `scripts/compose_auto.R`, `scripts/auto_criterion.R`
  the "auto" arm of the alfak2 era: picks, links, criterion; `scripts/beta_*.R`,
  `compose_beta.R`, `make_beta_tasks.R`, `prune_beta_candidates.R` its
  missegregation rate picked per input
* `R/paper_fig456.R`, `scripts/fig4_cna_effects.R`, `scripts/fig5_wgd.R`,
  `scripts/fig6_steady_states.R`, `scripts/fig456_summary.R` the paper's
  Fig. 4–6 on every fit
* `docs/report_<commit>.md` results of each KARST round;
  `docs/alfak2_era/` the reports of the alfak2 era
* `hpc/` image, side library, Slurm scripts, `link_alfak.sh`,
  `sacct_jobs.sh` (peak memory from Slurm's accounting) and `jobs/` (the
  job registries of the alfak2-era rounds, not tracked)
