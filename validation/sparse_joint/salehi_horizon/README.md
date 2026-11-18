# Salehi forecasts against later passages

Started 2026-11-16 at the user's request: use the Salehi passages as
two-time-point inputs, (a, b), (a, c), (a, d), ..., and the passages that
follow each input as the truth for its forecast. The ALFA-K paper's Fig. 2d
scores only the next passage, by a rank correlation that sampling noise
dominates (`../cloud_level/README.md`, step 2). Here every later passage is a
test, and the score is a proper one.

## Design

- **Inputs and fits.** The "ends" fits of a round of the paper benchmark
  (`benchmark/paper`), read only: KARST at its default (`initial_state =
  "balance"`), alfak (alfak dev2@f93b2d2) and KARST with `initial_state =
  "free"` (method `karst_free`, fitted by the same round into a folder of its
  own); minobs 5, 10 and 20. The test makes no fit of its own. The first run
  (2026-11-16, alfak2 era) read round 05a8855 and the free fits of the
  alfak2-salehi-variants project; the KARST 1.0.0 run reads round e7e32f8.
- **Forecasts.** From the trajectory's last passage, as the paper benchmark
  forecasts (`forecast_salehi_dir()`: `alfakR::predict_evo`, ABM, 5
  replicates, missegregation 5e-5, population 1e7), but up to the furthest
  descendant of that passage within four levels of the metadata's parent
  links (5 days per passage for SA906 and SA039, 15 for the others).
- **Truth.** Every descendant passage, from the metadata's parent links; its
  counts are the last column of a trajectory that ends with it.
- **Horizon.** Days between the parent and the later passage. A level of
  the parent links (`depth` in the rows) is a sampled passage, not one
  passage: adjacent levels are one passage (15 days) apart in most PDX lines
  (two, 30 days, in 13 rows of SA532 and SA609), but 5 to 21 passages (25 to
  105 days) apart in the p53 lines SA906 and SA039. Up to 2026-11-18 the
  tables were split by level and called it "passages ahead"; 117 of the 471
  rows of "one passage ahead" were 25 to 105 days ahead, and "3–4" held p53
  k.o. rows 75 to 175 days ahead. The tables below are by days: up to 15
  (one PDX passage), 16–30, 31–60, over 60. Beyond 60 days the rows are
  mostly of the p53 k.o. line (108 of 131), and the lineage bootstrap has
  three clusters (two in the all-time-points table), so those intervals are
  rough.
- **Scores**, on the karyotypes of the fit's landscape (the same set for all
  three methods), for each descendant:
  - the mean log-likelihood per cell of its counts under the forecast,
    minus that under the parent's frequencies (persistence), each mixed with
    a uniform distribution over the landscape with weight 0.01 (and 0.001).
    Sampling noise of the later passage enters both alike, so a forecast
    gains only where it puts probability on where the cells go. Persistence
    is the parent's raw frequencies, while the forecast is an ABM run with
    missegregation (5e-5 per copy and division), which moves some
    probability to one-step neighbours of the parent's karyotypes even
    without selection. Part of a forecast's gain over persistence can come
    from that spreading alone, so "beats persistence" is not by itself
    evidence of selection being forecast; the differences between methods
    (KARST − alfak), which share the engine and the spreading, are not
    affected. A persistence baseline with the same spreading (the ABM at
    flat fitness) was not run;
  - whether the forecast is closer to it than the parent in Wasserstein
    distance;
  - the distance between two random halves of it (its sampling noise).
- **Comparison.** On the rows that the compared methods share (trajectory,
  minobs, later passage), paired; 95% intervals from a bootstrap over the six
  lineage types and over the samples.

## Files

- `make_tasks.R`, `forecast_score.R`: the task table and one task (one fit),
  run on HPC in the paper benchmark's image (`benchmark/paper`:
  `hpc/submit.sh horizon TASKS.tsv 8G 02:00:00`, with `HORIZON_OUT` the rows'
  folder).
- `reuse_rows.R`: carries a method's rows over from an earlier run under its
  current name (alfak's rows of the first run, whose fits are the same).
- `analyse.R`, `analyse_all.R`: the tables and figures in `results/report/`
  and `results/report_all/` (not tracked); an optional argument names the
  results folder (default `results/`). The local copies of round e7e32f8
  are in `results/e7e32f8/` (reports only, with the paired rows in
  `rows.Rds`; the per-task rows are on HPC).
- Fits of the two designs ("ends", "all") are scored on landscapes of
  different karyotypes: the "all" landscape is larger, so the uniform floor
  (0.01 / K) and the persistence baseline differ, and the gains of the two
  designs are not comparable. From 2026-11-18 `forecast_score.R` also writes
  scores on a common support (`ll_fc_s_*`, `ll_p0_s_*`: the union of the
  landscapes of every fit of the trajectory and minobs, both designs, under
  the task's `support_roots`, which `make_tasks.R` sets), and
  `analyse_all.R` compares the designs on those only; with rows written
  before that it leaves the comparison out.

## KARST 1.0.0 rerun (2026-11-18)

With round e7e32f8 of the paper benchmark (KARST 1.0.0), 1662 tasks for the
two-time-point fits (KARST, KARST free) and 2460 for the fits of every passage
(anchored, pairwise, state space, free); alfak's rows were carried over from
the first run (`reuse_rows.R`), its fits being the same folders. Every table of
`results/report/` and `results/report_all/` equals that of the first run to
within 1e-14: the numbers below stand for KARST 1.0.0.

## Results (2026-11-16)

All 2493 tasks finished (array 21981908): 396 KARST, 402 alfak and 396
KARST "free" fits had later passages; the others had no fit (inputs the
benchmark excluded, or minobs without a fit) or no descendant. 986 rows (fit ×
later passage) on 143 trajectories are shared by all three methods.

**KARST forecasts the later passages better than alfak at every horizon.**
Log score per cell over persistence (mean; lineage-type bootstrap), by days
ahead of the parent. Recomputed 2026-11-18 with the current `analyse.R` from
the local rows of round e7e32f8 (`results/e7e32f8/report/rows.Rds`); the
pooled numbers are those of the first run:

| Days ahead | Rows (trajectories) | KARST | alfak | KARST free | KARST − alfak | Rows where KARST is better |
| --- | --- | --- | --- | --- | --- | --- |
| ≤ 15 (one PDX passage) | 341 (99) | +0.57 (0.25, 0.74) | +0.09 (−0.46, 0.31) | +0.52 | +0.48 (0.26, 0.91) | 57% |
| 16–30 | 257 (84) | +0.50 (−0.03, 0.86) | −0.55 (−1.64, −0.01) | +0.38 | +1.05 (0.55, 1.91) | 71% |
| 31–60 | 257 (68) | +0.02 (−0.91, 0.55) | −1.32 (−2.39, −0.55) | −0.02 | +1.35 (1.04, 1.62) | 74% |
| > 60 | 131 (27) | −1.69 (−2.36, −1.41) | −3.31 (−4.24, −2.99) | −1.67 | +1.61 (0.64, 2.83) | 75% |
| All | 986 (143) | | | | +1.01 (0.67, 1.40); samples (0.71, 1.26) | 68% |

Split by level of the parent links instead (the earlier table, "passages
ahead" 1 / 2 / 3–4): KARST +0.47 / +0.27 / −0.81, alfak −0.16 / −1.10 /
−2.14, KARST − alfak +0.63 / +1.37 / +1.33. On the 354 level-1 rows
outside the p53 lines (all but 13 of them one passage, 15 days, ahead)
KARST scores +0.60 (0.37, 0.75) and alfak +0.10; the direction is the same
by levels and by days.

- The difference holds at every minobs (5: +1.53, 10: +1.02, 20: +0.26, all
  intervals above zero), whether or not the treatment changed between the
  parent and the later passage (+0.67 and +1.05), and for inputs of one
  interval (+1.36) and of more (+0.81).
- Within 15 days (one PDX passage), KARST's forecast scores above
  persistence (interval above zero); alfak's does not. From 16 to 60 days
  neither does significantly; alfak falls below it (intervals below zero),
  KARST does not. Beyond 60 days (mostly the p53 k.o. line) both fall below
  persistence. How much of KARST's margin over persistence is the ABM's
  missegregation spreading and not selection is not known (Design, Scores).
- In Wasserstein distance within 15 days, KARST's forecast is closer to the
  later passage than the parent in 59% of rows, alfak's in 21%. Both errors
  are close to the sampling noise of the later passage (median excess 0.17
  for KARST, 0.26 for alfak, 0.17 for the parent).
- KARST with `initial_state = "free"` is no better than the default
  (−0.06 per cell, interval −0.16 to 0.09; within 15 days −0.06, interval
  −0.08 to 0.02; it is worse in 56–68% of rows, by horizon): on these
  passages the balance prior does not cost accuracy, which supports keeping
  it as the default.
- With a mixture weight of 0.001 instead of 0.01 the ordering is the same.

**Reading.** alfak's forecasts move the population further than it moves
(its fitness differences are larger: `../cloud_level/README.md` and the
paper benchmark), and a proper score penalises probability put where the
cells do not go. Neither method forecasts the composition more than 15 days
ahead significantly better than assuming no change; KARST is closer to that
baseline. This is the comparison on real data that sampling noise does not
decide: the forecast's probability of the cells actually seen.

![Forecasts against later passages](results/report/figures/horizon.png)
(figure in `results/report/figures/`, not tracked)

## All time points (2026-11-16)

Requested by the user the same day: the same test for the fits of every
passage of a trajectory (the paper benchmark's "all" design, trajectories of
three or more passages): KARST anchored (the default), pairwise and
state_space, alfak, and KARST "free"; beside them the two-time-point fits of
the same trajectories. No new fits; 3075 tasks (array 21989042,
`make_tasks.R ... all`, outputs `results/rows_all/`), analysis `analyse_all.R`
(`results/report_all/`). 677 rows on 99 trajectories are shared by the four
main arms.

Log score per cell over persistence (lineage intervals), by days ahead;
recomputed 2026-11-18 with the current `analyse_all.R` from the local rows
(`results/e7e32f8/report_all/rows.Rds`, days from the two-time-point rows of
the same trajectory and later passage):

| Days ahead (rows) | Anchored | Pairwise | State space | alfak |
| --- | --- | --- | --- | --- |
| ≤ 15 (259) | +0.92 (0.57, 1.07) | +0.93 (0.56, 1.08) | +0.93 (0.56, 1.07) | +0.38 (0.31, 0.62) |
| 16–30 (184) | +0.96 (0.22, 1.44) | +0.94 (0.20, 1.43) | +0.91 (0.25, 1.33) | −0.21 (−0.57, 0.12) |
| 31–60 (168) | +0.61 (−0.10, 1.30) | +0.47 (−0.27, 1.14) | +0.18 (−0.25, 0.67) | −1.40 (−2.09, −0.68) |
| > 60 (66) | −1.62 (−2.01, −1.56) | −1.84 (−2.63, −1.71) | −1.98 (−3.03, −1.81) | −2.58 (−2.67, −2.57) |
| Minus alfak, all rows | +1.12 (0.55, 1.47) | +1.06 (0.47, 1.44) | +0.96 (0.40, 1.31) | |

Beyond 60 days the 66 rows are of two lineage types (57 p53 k.o., 9 SA532):
those intervals mean little. The earlier table by level of the parent links
(1 / 2 / 3–4) had anchored +0.78 / +0.74 / −0.01 and alfak +0.25 / −0.66 /
−2.18.

- Every KARST method scores above alfak in each days bin (minus alfak:
  within 15 days +0.53 to +0.55, 16–30 days +1.11 to +1.17, 31–60 days +1.58
  to +2.01, intervals above zero), except state space beyond 60 days (+0.60,
  interval −0.36 to 0.75); in 69% of the rows; the sample bootstrap agrees.
- Among KARST's methods the anchored default scores highest: state space
  −0.15 (−0.20, −0.05), pairwise −0.06 (−0.10, −0.03), independent initial
  abundances −0.15 (−0.24, −0.01).
- This reverses the paper's Fig. 2d with all time points, where alfak
  ranked the size of the change at 0.55 and the anchored method at 0.23
  (`../cloud_level/README.md`, step 2: that ranking follows sampling noise).
- **Fitting every passage against two time points: not established.** The
  first run reported, on 636 rows of the same trajectories, that the anchored
  fit of every passage scores 0.36 per cell above the two-time-point fit
  (0.14, 0.49; higher in 90% of rows, at every horizon) and alfak's changes
  by +0.01 (−0.71, 0.91). Those are each fit's gain over persistence on its
  own landscape: the "all" landscape is larger, so the uniform floor
  (0.01 / K) and the persistence baseline differ between the designs (by
  0.52 per cell on average in these rows), and part of the difference can be
  mechanical. Not recomputed with the current code: it needs
  `forecast_score.R` run again (forecasts are not kept) for the "ends" and
  "all" tasks of these trajectories, with the common support
  (`make_tasks.R` sets it); until then `analyse_all.R` leaves this
  comparison out.

![All time points](results/report_all/figures/horizon_all.png)
(figure in `results/report_all/figures/`, not tracked)
