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

## Code adapted from the paper

Some files of this folder adapt the analysis code of the ALFA-K paper by
Richard J. Beck (https://github.com/Richard-Beck/ALFA-K), parts of it
unchanged:

| File | Adapted from |
|---|---|
| `scripts/fig4_cna_effects.R` | `scripts/06_cna_fitness_effects.R` |
| `scripts/fig5_wgd.R` | `scripts/07_wgd.R`, `dev/wgd4Noemi.R` |
| `scripts/fig6_steady_states.R` | `scripts/04_analyze_fit_data.R` |
| `R/paper_fig456.R` | `R/utils_karyo.R` (karyotype helpers, unchanged) |
| `R/metrics_salehi.R` | `scripts/04_analyze_fit_data.R`, `R/utils_karyo.R` |
| `R/metrics_abm.R` | `scripts/S04_process_abm_results.R`, `R/utils_karyo.R` |

The ALFA-K repository states no license. The adapted code remains the work of
its author and is not covered by the MIT license of karst (`LICENSE.md`),
which applies to the karst package and the rest of this repository.

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
picked per input (`beta = "auto"`, removed after round ac36fd1). The "auto"
arm is no method of the current code any more (`karst_auto` has no entry in
`R/fit.R` or `R/labels.R`, and `hpc/collect_analyse.sh` no step for it); its
scripts (`scripts/auto_picks.R`, `compose_auto.R`, `auto_criterion.R`) are
kept for the record only and do not run with the current code. The scripts
of the missegregation-rate arm are kept (`beta_*.R`, `compose_beta.R`,
`make_beta_tasks.R`, `prune_beta_candidates.R`); its steps run only when
their tables exist (`hpc/collect_analyse.sh`).

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
| Fig. 5b–c: WGD, landscape shape and fitness effects (SA906) | same; KARST one minobs per trajectory (`scripts/fig5_wgd.R`) |
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

Since 2026-11-18 (`fig456_conclusions()` in `R/paper_fig456.R`):

* A statement agrees with the reference when both have the same
  significance and, where both are significant, the same sign. Two results
  that are not significant agree whatever the signs of their estimates.
  (Before, the signs had to agree too, so that a non-significant estimate
  near 0 could flip the verdict from one round to the next.)
* Intercepts (4f, 5b) are not statements: the paper draws no conclusion
  from them, and the 5b intercept depends on each method's absolute fitness
  level. Each arm has 16 statements: 4d (4), 4f (6), 5b (5), 5c (1).
* Fig. 6 is descriptive (`agrees` NA): the share of landscapes whose
  dominant karyotype changes with the missegregation rate, beside the
  reference's. The former check, whether any landscape switches, holds for
  every arm by construction.
* Fig. 5: the paper's 07 keeps one fit per trajectory and minobs. alfak and
  the paper's fits keep that rule; KARST derives every minobs from one fit,
  so its minobs repeat a trajectory's karyotypes with the same fitness
  (round e7e32f8: KARST 2 time points, 19866 rows of which 10723 distinct
  trajectory-karyotype pairs). KARST keeps one minobs per trajectory, as in
  Fig. 4: the highest CV score, of a tie the lowest minobs
  (`fig5_select_fits()`; `fig5_tests.tsv`, column `minobs_rule`).
* Fig. 4d: the variance of delta f per frequent karyotype is taken without
  `na.rm`, as the paper's 06 does (line 336), so a frequent karyotype with a
  neighbour missing from the landscape drops out of the variance models
  (`PDX_kv`, `cisplatin_kv`). The paper's line is kept to reproduce its
  numbers; each draw records the share dropped
  (`analysis/fig456/fig4_kvar_dropped.tsv`) and the two models with
  `na.rm = TRUE` (`PDX_kv_narm`, `cisplatin_kv_narm` in
  `fig4_estimates.tsv`; no reference, so not a statement). The shares have
  not been computed with the current code yet (the bootstrap runs on the
  HPC; next round).

Two ways of reporting: the paper's rule (failed fits dropped, test
instances selected by each method's own CV score) and a paired comparison
on the inputs where both methods have a fit, without selection.

A fit failed when its run did (`run.Rds` status not "ok"; for the paper's
own Salehi rows, 03's NaN CV score of a fit without `xval.Rds`), in every
table: `fit_failures*.tsv`, `excluded_inputs.tsv` and
`summary_performance.tsv` (`salehi_fit_ok()`, `R/metrics_salehi.R`). A fit
that ran but has no CV score (CV score NA; on the Salehi data of round
e7e32f8, 65 alfak fits and 1 KARST fit, 28 of them with a forecast; on the
analysable inputs 53 and 26) is not a failure: it is left out only where a
fit is selected by its CV score, and both Salehi panels of
`summary_comparison` (forecasts; novel karyotypes) pair every fit that
ran. Until 2026-11-18 the Salehi tables counted a missing CV score as a
failure and dropped these fits from the paired forecasts.

The paper's rule keeps every minobs that ties for an input's highest CV
score (S05: `Rxv == max(Rxv)`). alfak's scores never tie; KARST derives every
minobs from one fit, and two or three of them often have the same
cross-validation predictions and scores that differ only by rounding, so that
a strict comparison counted such an input once, twice or three times
depending on the rounding. Since 2026-11-18 scores within 1e-8 of the best
count as tied and the lowest minobs of a tie is kept (`R/select.R`), in every
analysis that selects by CV score (Supp. Fig. 4c and 5, Fig. 1g, 2b-2f, 3e-f,
4 and 6). For alfak and the paper's own fits the selection is the paper's.

## Fig. 2d and sampling noise

`scripts/fig2d_noise.R` takes one row per trajectory, minobs and daughter
passage, on the rows every arm of a design has. A daughter passage is the
evaluation passage of every trajectory that ends at its parent, for each
minobs, so the same daughter gives many rows. Round e7e32f8: two time
points, every trajectory, 471 rows from 57 daughters (up to 24 rows per
daughter); one interval, 147 rows from 53 daughters (up to 3, the minobs);
all time points, 3–10 passages, 346 rows from 46 daughters (up to 21). The
rows are not independent. The lineage-type bootstrap (6 types) resamples
whole types and so widens the intervals accordingly, but the row counts
overstate the information.

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
and `alfak2_state_space`, set `COMPARE_OLD_METHOD` to that name. A fit that
only one of the two rounds has is listed with status "missing" on the other
side. `hpc/link_alfak.sh` stops when the reference has no `fits` folder or
no folder of the method to link.

| KARST | Results | Report |
|---|---|---|
| b9a275e (1.0.0, the rename; stopped with every fit but one done, superseded by e7e32f8) | `results_b9a275e/` | |
| e7e32f8 (1.0.0 with the fixes of the code review of 2026-11-17; all KARST fits redone) | `results_e7e32f8/` | `docs/report_e7e32f8.md` |
| 2c60af1 (1.1.0, fixes of the code review of 2026-11-18; complete, but 85 KARST fit tasks failed when one candidate history's posterior mode was not found) | `results_2c60af1/` | |
| 50415cf (1.1.1, the posterior mode searched once more; stopped soon after submission, superseded by cd342db) | `results_50415cf/` (partial) | |
| cd342db (1.1.2, candidate histories compared on the objective at the posterior mode; all KARST fits redone) | `results_cd342db/` | `docs/report_cd342db.md` (also compares with 2c60af1) |
| 6ae2a24 (1.1.3, a fit's own posterior rebuilt at its stored mode, for the draws of alfak's files; all KARST fits redone) | `results_6ae2a24/` | `docs/report_6ae2a24.md` (also compares with cd342db) |

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

Before they collect anything, `collect_abm.R` and `collect_salehi.R` check
every fit of the round's task tables (`R/tasks.R`; `BENCH_TASKS`, default
`tasks/<commit>`, `none` to skip; `karst_free`'s tasks are checked in
`results_<commit>_free`, the others in `results_<commit>`): `run.Rds` in every
minobs folder, for a fit that ran `score.Rds` and (ABM) `forecast.Rds`, for
a Salehi fit that ran `forecast.Rds` where the trajectory has a daughter
passage, and Salehi scores for every fit folder with `run.Rds`
(`score_salehi` must run after the last fit). A task killed for time or
memory before KARST's fit ended leaves no folder at all (KARST writes its
minobs folders after the fit), so without this check it would vanish from
the tables, while alfak's killed fits count as failed. Anything missing
stops the step and is listed in `tables/missing_{abm,salehi}.tsv`: rerun
those tasks. `ALLOW_MISSING=1 hpc/collect_analyse.sh` collects anyway, a fit
without `run.Rds` counting as failed (status "missing").

`run_task.R` can be run again on a task after any failure. Its files
(`run.Rds`, forecast replicates, `forecast.Rds`, `score.Rds`) are written
under a temporary name and renamed, so that a killed task leaves no
truncated file that a rerun would take as done. A rerun of a KARST task
fits only the minobs without `run.Rds`, so that finished minobs keep their
bootstrap draws, forecasts and scores; it resumes a forecast that did not
finish (replicates already written are kept), and scores a fit again when
its forecast is rerun. Only forecasts with status "ok" in `forecast.Rds` are
scored (ABM and Salehi); replicates of a failed or unfinished forecast are
never scored. A fit whose score cannot be computed makes the task exit with
status 1 (Slurm: FAILED) after the other steps.

KARST's tiers: 4 GB and 6 h for two time points of the ABM and
missegregation data; 8 GB and 12 h for the Salehi data and the anchored fits
of all time points; 16 GB for pairwise and 32 GB for state space (12 h).
`karst_free` takes the tiers of the default (before 2026-11-18,
`split_tasks.R` gave it alfak's).

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
* `R/select.R` selection by CV score, the paper's rule with ties resolved
  (an option without a score is never chosen; a missing value in a
  grouping column is an error)
* `R/tasks.R` the round's task tables and the check of their outputs before
  collecting
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
  with the passages' sampling noise held fixed (uses `scripts/fig2d_check.R`;
  see "Fig. 2d and sampling noise")
* `scripts/example_abm.R` the paper's Supp. Fig. 2c–e for every fit
* `scripts/auto_picks.R`, `scripts/compose_auto.R`, `scripts/auto_criterion.R`
  the "auto" arm of the alfak2 era: picks, links, criterion (obsolete, kept
  for the record); `scripts/beta_*.R`,
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
