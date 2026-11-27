# alfak2 348addd versus alfakR on the ALFA-K paper tests: comparison report (official results)

Date: 2026-11-17. alfak2 348addd (dev), alfakR dev2 f93b2d2. Results directory `results_348addd/`.

> Terms: the input columns are called "time points". A Salehi "trajectory" consists of several passages; SA609 and the like are "lines".

## 0. What this round tests

| commit | Change | Effect on this benchmark |
|---|---|---|
| 827c4b5 | With two time points, the default is `initial_state = "auto"`: first run a check fit on the interval; if the overall level (log level) of the non-major karyotypes is above 1, switch to `"free"`, otherwise keep `"balance"` | On one-passage windows the landscape accuracy drops clearly, see `docs/report_827c4b5.md`; not taken as official results |
| 348addd | Default changed back to `initial_state = "balance"`; `"auto"` kept as an option, with the documentation stating its limitation on one-passage intervals | All fits should be the same as in 05a8855 |

User decision: this round does not refit everything. First refit a verification subset and compare it fit by fit with 05a8855; if they are the same, link the remaining fits from `results_05a8855/`. Scoring, summaries and Fig. 4–6 are rerun in this round's directory.

## 1. Conclusions

1. **The results of 348addd are the same as those of 05a8855.**
   - The 104 verification tasks are identical to 05a8855 one by one: same status, landscape differences ≤ 1.1e-16, cross-validation and ABM forecasts exactly the same (Section 2). These 104 include the inputs that 827c4b5 switched to free.
   - The 114 rows of accuracy values in the summary table are the same as in 05a8855 row by row. The conclusions of Fig. 4–6 are also the same. The only difference is that this round has no auto arm.
   - Therefore the results of the paper tests are still those in Sections 3–6 of `docs/report_05a8855.md`, and this report does not repeat them.
2. **Salehi's Fig. 2d mainly reflects sampling noise.**
   - Splitting the cells of one passage at random into two halves, the Spearman between the distance between the two halves and the observed change is 0.73–0.77, higher than for any method's forecast.
   - After fixing this noise (partial Spearman), no alfak2 arm differs reliably from alfakR.
   - In the raw Fig. 2d alfakR leads (all time points: 0.55 vs anchored 0.23); this lead cannot be taken as evidence of a more accurate landscape (Section 4).
3. **When forecasting later passages, alfak2 beats alfakR in both designs.** This is the test from the multi-time-point session, reading the fits of 05a8855 (Section 5).
   - Two time points: the log score per cell is 1.01 higher than alfakR's (lineage interval 0.67 to 1.40), and higher at every forecast horizon.
   - All time points: anchored is 1.12 higher (0.55 to 1.47). Among the three alfak2 arms, anchored is best.
   - 1 passage ahead, alfak2 beats "no change" and alfakR does not.
4. **The lesson of 827c4b5**: with two time points, the check fit sets the history to the interval itself. On a one-passage window, the measured overall level mainly reflects that "the true history is much longer than the interval", not that the missegregation rate is set wrongly.
   - free was used even when the missegregation rate was set correctly, so the landscape accuracy dropped.
   - See `docs/report_827c4b5.md` for details.

## 2. Fits: verification subset plus links

Verification subset (`tasks/348addd/`):

| Group | Inputs | Notes |
|---|---|---|
| ABM w2_ends | 28 | inputs that 827c4b5 switched to free, 8 each for λ 0.4 / 0.8 / 1.6; 4 for λ 0.2 (not switched) |
| misseg w2_ends | 20 | inputs that 827c4b5 switched to free, 4 for each of the 5 true missegregation rates |
| Salehi ends | 20 | 16 switched to free by 827c4b5, 4 not switched |
| Controls | 28 | 4 each for ABM w4/w8 ends and all, 2 each for misseg, 4 for Salehi all |
| pairwise / state_space | 4 + 4 | ABM w4_all |

`scripts/compare_rounds.R` compares against `results_05a8855` (`analysis/compare_rounds/`):
- every fit has the same status, and for fits that succeeded on both sides the landscape differences are all ≤ 1.1e-16;
- cross-validation and ABM forecasts are exactly the same, and all fits use `"balance"`;
- fits that succeeded on neither side number 11 by minobs: 8 Salehi ends, 3 Salehi all. Their failure status agrees with 05a8855.

The links were made after all verification tasks had finished: 3011 for alfak2, 1251 for pairwise, 1251 for state_space, 3107 for alfakR (from `results/`). Scoring and analysis write only to this round's `scores/`, `tables/` and `analysis/`, not into the linked fit folders.

## 3. Paper tests

This round rescored and resummarised; the results are the same as in 05a8855 row by row, see `docs/report_05a8855.md`.
- Section 3: summary accuracy. The 114 rows of `analysis/summary/summary_accuracy.tsv`, values and pair counts all the same.
- Section 4: the individual paper tests.
- Sections 5–6: Fig. 4–6. Number of conclusions agreeing with the paper:
  - alfak2 two time points 11/19, anchored 10/19, pairwise 9/19, state_space 10/19;
  - alfakR two time points 14/19, all time points 15/19;
  - the paper's own fits (recomputed with this code) 19/19.

## 4. Fig. 2d: raw and after fixing sampling noise (`analysis/fig2d_noise/`)

Fig. 2d computes a Spearman correlation across trajectories between two distances of each forecast at the next passage: the observed change (parent to daughter) and the predicted change (parent to forecast), both in Wasserstein distance.

**The paper's rule**: fits selected by CV score (`analysis/salehi/fig2d_wasserstein.tsv`).

| Fit | Rows | Spearman |
|---|---|---|
| alfak2 two time points | 53 | 0.25 |
| alfak2 all time points (anchored) | 52 | 0.33 |
| alfak2 pairwise | 53 | 0.57 |
| alfak2 state space | 49 | 0.56 |
| alfakR two time points | 44 | 0.45 |
| alfakR all time points | 46 | 0.70 |
| Paper | 45 | 0.69 |

The subsets selected by the paper's rule differ, the row counts are small, and they cover only 6 lineage types, so the arms cannot be compared directly.

**Rows shared by all arms**, without CV selection (`scripts/fig2d_noise.R`):
- Noise = after splitting the parent and the daughter each at random into two halves, the mean distance between the halves (20 splits), computed on the karyotypes of the alfakR landscape.
- Interval = bootstrap 95% interval over the 6 lineage types.

| Design and subset | Rows | Noise itself | alfak2 raw | alfak2 noise fixed | alfakR raw | alfakR noise fixed |
|---|---|---|---|---|---|---|
| Two time points, all trajectories | 471 | 0.75 | 0.29 | 0.47 | 0.32 | 0.44 |
| Two time points, one interval (2 passages) | 147 | 0.77 | 0.45 | 0.49 | 0.60 | 0.55 |
| All time points, 3–10 passages: anchored | 346 | 0.73 | 0.23 | 0.40 | 0.55 | 0.52 |
| Same as above: pairwise | 346 | | 0.35 | 0.55 | | |
| Same as above: state_space | 346 | | 0.45 | 0.56 | | |

Difference from alfakR (alfak2 minus alfakR):

| Design and subset | Arm | Raw | Noise fixed |
|---|---|---|---|
| Two time points, all trajectories | alfak2 | −0.03 (−0.19 to 0.15) | +0.03 (−0.16 to 0.21) |
| Two time points, one interval | alfak2 | −0.15 (−0.19 to −0.02) | −0.07 (−0.25 to 0.04) |
| All time points | anchored | −0.32 (−0.53 to 0.09) | −0.12 (−0.41 to 0.13) |
| All time points | pairwise | −0.20 (−0.37 to 0.10) | +0.03 (−0.10 to 0.14) |
| All time points | state_space | −0.10 (−0.27 to 0.11) | +0.04 (−0.28 to 0.16) |

After fixing the noise, no arm differs reliably from alfakR. The differences in the raw Fig. 2d are mainly about how closely the size of the forecast follows the noise. This agrees with the result of the multi-time-point session (alfak2 `validation/sparse_joint/cloud_level`, step 2): one interval −0.07 (−0.25 to 0.04), anchored −0.13 (−0.46 to 0.12), state_space +0.03. This round's recomputation gives the same numbers.

## 5. Forecasting later passages (test from the multi-time-point session)

**Source**
- alfak2 b9a3b9e, `validation/sparse_joint/salehi_horizon/`. The tables are copied to `build/salehi_horizon_from_alfak2/` on HPC; for the meaning of the columns see `COLUMNS.md` there.
- Only the fits of 05a8855 are read, the same as in this round. Forecasts are redone with this benchmark's forecast engine: the ABM of `alfakR::predict_evo`, 5 replicates, pm 5e-5, nb 1e7.

**Design**
- Within 4 passages after the last passage of the input, every daughter passage serves as truth.
- Score: log score per cell, relative to "no change" (the parent's frequencies). A positive number means the forecast put more probability where the cells actually are.
- Computed on the karyotypes of the fitted landscape; the three methods have the same landscape karyotypes.
- Interval: bootstrap 95% interval over lineage types.

**Two time points**: the 986 rows shared across 143 trajectories.

| Passages ahead | alfak2 vs "no change" | alfakR vs "no change" | alfak2 minus alfakR | Rows where alfak2 is higher |
|---|---|---|---|---|
| 1 (471 rows) | +0.47 (0.17 to 0.67) | −0.16 (−0.76 to 0.22) | +0.63 (0.34 to 0.99) | 65% |
| 2 (280 rows) | +0.27 (−0.44 to 0.79) | −1.10 (−2.17 to −0.20) | +1.37 (0.90 to 2.04) | 70% |
| 3–4 (235 rows) | −0.81 (−1.96 to 0.43) | −2.14 (−3.53 to −0.73) | +1.33 (1.10 to 1.63) | 71% |
| All | | | +1.01 (0.67 to 1.40) | 68% |

- By the number of intervals in the input: one interval +1.36 (1.10 to 1.75), more than one interval +0.81 (0.38 to 1.33).
- 2 or more passages ahead, neither method beats "no change".

**All time points**: 99 trajectories with 3 or more passages, 677 rows shared by the four arms.

| Passages ahead | anchored | pairwise | state_space | alfakR all time points |
|---|---|---|---|---|
| 1 (346 rows), vs "no change" | +0.78 (0.35 to 1.01) | +0.78 | +0.77 | +0.25 (−0.15 to 0.61) |
| 2 (191 rows) | +0.74 | +0.69 | +0.63 | −0.66 (−1.35 to −0.25) |
| 3–4 (140 rows) | −0.01 | −0.22 | −0.55 | −2.18 (−2.79 to −1.80) |

| Paired difference (all horizons) | Difference | Interval | Rows where the first is higher |
|---|---|---|---|
| anchored minus alfakR | +1.12 | 0.55 to 1.47 | 69% |
| pairwise minus alfakR | +1.06 | 0.47 to 1.44 | 69% |
| state_space minus alfakR | +0.96 | 0.40 to 1.31 | 69% |
| state_space minus anchored | −0.15 | −0.20 to −0.05 | 35% |
| pairwise minus anchored | −0.06 | −0.10 to −0.03 | 46% |
| free minus anchored | −0.15 | −0.24 to −0.01 | 32% |

All time points versus two time points, same set of trajectories (636 rows):

| Method | Difference | Interval | Rows where the first is higher |
|---|---|---|---|
| alfak2 (anchored vs two time points) | +0.36 | 0.14 to 0.49 | 90% |
| alfakR | +0.01 | −0.71 to 0.91 | 54% |

The two designs choose landscape karyotypes in different ways, so across designs only each one's gain relative to "no change" is compared, not the raw log-likelihood. This test reverses the ranking of the raw Fig. 2d on all time points: in Fig. 2d it is alfakR 0.55 vs anchored 0.23, whereas for forecasts of later passages anchored is 1.12 higher than alfakR.

## 6. Runs

Job IDs are recorded in `build/348addd_jobs.tsv` on HPC:
- build 21986108;
- verification arrays 21986131 (alfak2, 96 tasks, 8G) and 21986133 (pairwise / state_space, 8 tasks, 32G);
- comparisons 21987401 / 403 / 405;
- scoring 21987448;
- collect_analyse 21987450;
- Fig. 4 21987451 / 453, Fig. 6 21987455 / 457, fig456_summary 21987459;
- Fig. 2d noise 21987525 / 527.

All jobs ran on qos xxlarge, 1 CPU per task, with an add-hours monitor for every array. No TIMEOUT, OOM or NODE_FAIL.

The 827c4b5 jobs still read `hpc/alfak2.lock` while queued, so this round specifies its own lock file via `ALFAK2_LOCK=hpc/locks/alfak2_348addd.lock` (`hpc/env.sh`). After 827c4b5 finished, `hpc/alfak2.lock` was changed to 348addd.

The alfak2 part of the resource table (`tables/resources.Rds`) is taken from the job records of 05a8855 (`RESOURCES_ARGS="alfakR=8327282 alfak2=05a8855"`), because the vast majority of the fits were run in that round.

## 7. Files

- `analysis/compare_rounds/results_05a8855*_summary.tsv`: comparison of the verification subset with 05a8855.
- `analysis/summary/`, `analysis/abm/`, `analysis/salehi/`, `analysis/fig456/`: the same tests as in 05a8855.
- `analysis/fig2d_noise/`: `fig2d_noise.tsv`, `fig2d_noise_differences.tsv`, `figures/fig2d_noise.pdf`.
- HPC `build/salehi_horizon_from_alfak2/`: the later-passage test of the multi-time-point session, including `all_time_points/`.
