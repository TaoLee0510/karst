# Missegregation-rate grid on Salehi: Fig. 2c / 2d for KARST and alfak (round 6ae2a24)

2026-11-24. KARST (this round's library, 1.1.3) and alfak (dev2@f93b2d2) fit and forecast with the same missegregation rate β, β = 5e-5 × {1, 2, 4, 8, 16, 32, 64}. ×1 is this round's own fit. Both designs were run: two time points (277 trajectories) and all time points (205).

## Conclusions

1. **At the default β (5e-5), KARST leads numerically on Fig. 2c but not significantly, and lags on Fig. 2d.**
   - Fig. 2c (no selection, paired): angle difference alfak − KARST is +1.4° (−0.1 to 4.9) for two time points and +1.1° (−0.3 to 4.3) for all time points.
   - Fig. 2d (paper rule): KARST 0.22, alfak 0.45.

2. **When both methods use the same β:**
   - **Fig. 2c, no selection, paired:** at every β the KARST angle is smaller. The interval lies entirely above 0 in only one place: ×64, all time points, +8.2° (0.3 to 23.1). The gap widens with β mainly because alfak gets worse (two time points 83.9° → 89.3°); KARST improves only slightly (79.7° → 73.9°).
   - **Fig. 2c, paper rule (the statistic used in the figure):** KARST is smaller only at ×1 (and at ×2 for two time points). From ×4 on, alfak is smaller at most β. Not significant at any β: there are only 5–6 lineage types, and the intervals are wide.
   - **Fig. 2d, paper rule:** KARST rises with β, from 0.22 to 0.74. From ×8 on for two time points and from ×32 on for all time points, KARST is above alfak, but all intervals include 0. The closest is ×64, all time points: +0.26 (−0.04 to 0.58).
   - **Fig. 2d, no selection (two time points, per trajectory):** on raw values KARST leads from ×16 on (0.46 vs 0.37 at ×32). After fixing the noise, KARST is slightly higher at every β (0.50 vs 0.36 at ×64), also not significant.

3. **Each method at its own best β** (picked on the same data, so an optimistic estimate):

   | test | design | KARST | alfak | difference (95% interval) |
   |---|---|---|---|---|
   | Fig. 2c, no selection | two time points | 73.9° (×64) | 82.6° (×16) | +2.3° (−1.1 to 8.0) |
   | Fig. 2c, no selection | all time points | 67.4° (×8) | 74.1° (×16) | +1.0° (−1.1 to 4.9) |
   | Fig. 2d, paper rule | two time points | 0.74 (×64) | 0.57 (×2) | +0.17 (−0.27 to 0.62) |
   | Fig. 2d, paper rule | all time points | 0.80 (×64) | 0.71 (×4) | +0.09 (−0.24 to 0.35) |

   With each at its best, KARST leads in point estimate on both tests, but neither is significant.

4. **Can KARST "lead":**
   - Fig. 2d and the no-selection Fig. 2c: yes in point estimate, as long as β is high enough (×16 to ×64).
   - A significant lead occurs in only one place: Fig. 2c no selection, all time points, ×64.
   - Fig. 2c drawn by the paper rule: raising β does not make KARST lead.

5. **"Data with high dispersion are more sensitive to β" is not supported across the Salehi cell lines.**
   - KARST's best β per cell line is unrelated to its dispersion (Spearman 0), and so is the improvement from switching to the best β (Spearman 0.10, p = 0.77, 10 cell lines).
   - The dispersed SA609, SA609R2 and SA532 improved by 11°–19°; but the low-dispersion SA1035U improved by 28°, and the most dispersed SA000 improved by 0.
   - The evidence for that statement comes from simulations (`docs/fig2d_misseg_6ae2a24.md`): when β is below the truth, the more dispersed the population, the larger KARST's loss.

6. **Suggested wording for reviewers** (Section 7): KARST is sensitive to the missegregation rate; on Salehi, Fig. 2d rises from 0.22 to 0.74 with β; give the default and the best values side by side, and note that the best value was picked on the same data.

7. Neither report Ver10 nor Ver11 has been changed; to be decided after you have looked at this.

## 1. Method

- **Fits:** this round's `run_task.R` and settings, with only `pm` changed to 5e-5 × multiplier, used for both fit and forecast. One KARST task covers minobs 5/10/20; alfak has one task per minobs. 1,928 tasks per β, results in one folder each.
- **×1:** links directly to this round's fits; the alfak fits themselves are also linked from the alfak2 era, with the same software version.
- **The later steps are the same as in this round:** `score_salehi.R`, `check_ksets.R salehi`, `collect_salehi.R`, `analyse_salehi.R`, `fig2d_noise.R`.
- **Comparisons (`scripts/beta_grid_compare.R`):**
  - **Fig. 2c, no selection:** pairs where both methods have a forecast on the same trajectory, the same minobs and the same evaluation passage. Difference = alfak angle − KARST angle, median taken; interval from a cluster bootstrap by Salehi cell line (2,000 times). Above 0 means KARST's predicted direction is closer to the actual one.
  - **Fig. 2c, paper rule:** each method picks minobs and trajectory length by its own CV score. The statistic is the same as in the figure: the median over lineage types of the median angle. Interval likewise from a bootstrap by cell line.
  - **Fig. 2d, paper rule:** Spearman of predicted and actual change (Wasserstein distance); difference bootstrapped by cell line.
  - **Fig. 2d, no selection:** as in `fig2d_noise.R`: rows that both methods have, raw values and partial Spearman after fixing the noise, bootstrapped by lineage type. For the all-time-points design there are only the KARST and alfak arms here, so the shared rows are more than in this round (346 rows), and the numbers of that group cannot be compared directly with this round's report.
  - **Best β:** three ways of reporting: the whole curve at the same β; each at its best; KARST's best β per cell line against that cell line's dispersion (`fig2d_misseg/dispersion.tsv`: fraction of karyotypes seen only once among 300 cells).

## 2. Checks

- **×1 reproduction:** rescoring ×1 with the current code is bit-identical to this round (spot check of 1,920 rows from trajectories with data, maximum difference 0). Fig. 2c, Fig. 2d (paper rule) and the noise-fixed Fig. 2d (two time points) are identical, word for word, to this round's report.
- **Out of memory:** at each β, 13–16 alfak tasks (all time points, minobs 5) exceeded 16 GB and were killed; at collection they were caught by the missing-output check. Rerun with 48 GB, all succeeded, peak 47 GB. Then the whole β was rescored, and the noise blocks and collection redone.
- **Timeout fallback for the Wasserstein distance:** the default networkflow of `transport::wasserstein` hangs on a few inputs (found in `fig2d_misseg`). The new `R/wasserstein.R` first uses networkflow and switches to revsimplex after 60 seconds; revsimplex is an exact solver of the same optimal transport problem. The number of tasks that triggered it during scoring grows with β: 1 at ×1, 93 at ×64. Calls that did not trigger it give unchanged results, so the ×1 reproduction is not affected.
- **Failure rates are the same at every β:**

  | | two time points | all time points |
  |---|---|---|
  | KARST | 14.4% | 11.2% |
  | alfak | 9.4% | 3.6% |

  KARST's failures are all determined by the data and do not depend on β: 126 domain errors, 63 "no frequent karyotypes at this minobs".

## 3. Fig. 2c

**No selection, paired** (median angle, smaller is better; difference = alfak − KARST):

| β | two time points KARST / alfak | difference (95% interval) | all time points KARST / alfak | difference (95% interval) |
|---|---|---|---|---|
| ×1 | 79.7 / 83.9 | +1.4 (−0.1 to 4.9) | 71.2 / 79.0 | +1.1 (−0.3 to 4.3) |
| ×2 | 78.5 / 84.3 | +1.8 (−0.3 to 6.5) | 69.1 / 77.3 | +1.7 (−0.3 to 4.7) |
| ×4 | 79.2 / 84.6 | +1.5 (−0.5 to 5.7) | 69.5 / 78.2 | +1.7 (−0.5 to 5.0) |
| ×8 | 79.5 / 84.8 | +0.8 (−1.0 to 7.4) | 67.4 / 76.5 | +1.4 (−0.4 to 5.3) |
| ×16 | 79.5 / 82.6 | +0.4 (−1.4 to 5.9) | 68.7 / 74.1 | +1.0 (−1.7 to 6.6) |
| ×32 | 74.6 / 85.3 | +3.0 (−1.1 to 7.4) | 70.6 / 82.4 | +4.9 (−1.5 to 9.3) |
| ×64 | 73.9 / 89.3 | +4.5 (−2.9 to 11.7) | 68.5 / 87.1 | **+8.2 (0.3 to 23.1)** |

Number of pairs: 394 for two time points, 409 for all time points, covering 11 cell lines. Pairs where the KARST angle is smaller make up 53%–63%.

**Paper rule** (median over lineage types of the median angle):

| β | two time points KARST / alfak | all time points KARST / alfak |
|---|---|---|
| ×1 | 72.1 / 83.2 | 66.0 / 79.9 |
| ×2 | 79.0 / 83.8 | 82.6 / 83.8 |
| ×4 | 82.3 / 48.0 | 78.1 / 64.7 |
| ×8 | 79.4 / 70.5 | 81.2 / 58.1 |
| ×16 | 82.0 / 55.3 | 86.2 / 57.0 |
| ×32 | 85.4 / 83.5 | 85.0 / 85.7 |
| ×64 | 94.6 / 71.7 | 91.8 / 87.1 |

- At every β the interval of the difference includes 0.
- This statistic is based on only 5–6 lineage types and also changes with the CV selection, so it is very unstable.
- The number of KARST fits selected falls from 50 to 31 as β rises.
- Figures: `analysis/figures/fig2c_grid.pdf` (no selection), `fig2c_rule_grid.pdf` (paper rule).

## 4. Fig. 2d

**Paper rule** (Spearman; difference = KARST − alfak):

| β | two time points KARST / alfak | difference (95% interval) | all time points KARST / alfak | difference (95% interval) |
|---|---|---|---|---|
| ×1 | 0.22 / 0.45 | −0.23 (−0.57 to 0.03) | 0.28 / 0.70 | −0.41 (−0.96 to −0.02) |
| ×2 | 0.24 / 0.57 | −0.33 (−0.77 to 0.06) | 0.23 / 0.64 | −0.41 (−0.95 to −0.08) |
| ×4 | 0.40 / 0.46 | −0.06 (−0.35 to 0.47) | 0.33 / 0.71 | −0.38 (−0.77 to −0.16) |
| ×8 | 0.50 / 0.42 | +0.08 (−0.29 to 0.55) | 0.38 / 0.57 | −0.19 (−0.84 to 0.16) |
| ×16 | 0.54 / 0.45 | +0.09 (−0.41 to 0.57) | 0.54 / 0.59 | −0.04 (−0.63 to 0.29) |
| ×32 | 0.56 / 0.27 | +0.29 (−0.34 to 0.76) | 0.56 / 0.46 | +0.10 (−0.49 to 0.41) |
| ×64 | 0.74 / 0.36 | +0.39 (−0.12 to 0.78) | 0.80 / 0.54 | +0.26 (−0.04 to 0.58) |

**No selection, two time points, per trajectory** (`fig2d_noise.R`, 471 rows):

| β | raw KARST / alfak | after fixing noise KARST / alfak |
|---|---|---|
| ×1 | 0.29 / 0.32 | 0.47 / 0.44 |
| ×2 | 0.28 / 0.32 | 0.46 / 0.42 |
| ×4 | 0.29 / 0.32 | 0.44 / 0.40 |
| ×8 | 0.33 / 0.34 | 0.42 / 0.40 |
| ×16 | 0.38 / 0.35 | 0.42 / 0.38 |
| ×32 | 0.46 / 0.37 | 0.43 / 0.30 |
| ×64 | 0.52 / 0.45 | 0.50 / 0.36 |

- Ranking by noise alone gives 0.75, still higher than either method at any β.
- The interval of KARST − alfak includes 0 at every β. At ×32, raw +0.09 (−0.08 to 0.26), after fixing noise +0.13 (−0.06 to 0.35).
- All time points (3–10 passages, 346 rows): on raw values KARST is above alfak only at ×64 (0.55 vs 0.50).
- Figure: `analysis/figures/fig2d_grid.pdf`.

## 5. Per cell line

Two time points, no selection, by KARST's median angle:

| cell line | pairs | dispersion (seen only once) | KARST best β | KARST default → best | alfak default → best (β) |
|---|---|---|---|---|---|
| SA039U | 9 | 0.03 | ×1 | 5.8 → 5.8 | 9.0 → 6.6 (×2) |
| SA1035U | 16 | 0.08 | ×64 | 135.1 → 107.6 | 127.4 → 87.4 (×64) |
| SA906b | 63 | 0.17 | ×2 | 75.7 → 74.5 | 83.9 → 73.8 (×64) |
| SA1035T | 18 | 0.19 | ×1 | 59.3 → 59.3 | 61.5 → 61.1 (×16) |
| SA906a | 45 | 0.19 | ×64 | 26.2 → 25.5 | 26.8 → 25.8 (×16) |
| SA532 | 62 | 0.26 | ×32 | 87.8 → 74.8 | 102.2 → 96.7 (×64) |
| SA609 | 92 | 0.58 | ×32 | 91.4 → 80.7 | 96.6 → 90.5 (×32) |
| SA535 | 47 | 0.63 | ×8 | 27.6 → 23.2 | 46.9 → 45.8 (×8) |
| SA609R2 | 18 | 0.78 | ×32 | 99.4 → 80.8 | 111.3 → 89.7 (×64) |
| SA000 | 11 | 0.81 | ×1 | 88.6 → 88.6 | 97.3 → 96.3 (×8) |
| SA609UnBU | 13 | no data | ×1 | 60.6 → 60.6 | 84.5 → 83.3 (×64) |

- The number of pairs per cell line is small, so the best β is largely a matter of chance.
- Figure: `analysis/figures/line_best_rate.pdf`.

## 6. Relation to earlier results

- **The 11-20 real-data β experiment** (E4, KARST fits only) also found that Salehi prefers a high forecast rate. After decomposition, the gain came mainly from a more dispersed forecast, not from the β used in the fit. Here β was changed in both fit and forecast, so KARST's improvement at high β may likewise come mainly from the dispersion of the forecast.
- **varyMisseg simulations** (`docs/fig2d_misseg_6ae2a24.md`): when β equals the truth, KARST leads alfak at every missegregation-rate level.
- Salehi has no truth, so the "best β" here only shows which β Fig. 2c/2d prefer, not which β is true.

## 7. Suggested wording for reviewers (for you to decide)

- "KARST uses a fixed missegregation rate of 5e-5. On Salehi, raising the missegregation rate of both KARST and alfak from 5e-5 to 3.2e-3 raises KARST's Fig. 2d from 0.22 to 0.74, with alfak between 0.27 and 0.57; the paired angle difference of Fig. 2c significantly favours KARST at ×64, all time points. Default and best values are reported side by side; the best value was picked on the same data."
- "In simulations, when the missegregation rate is below the truth, KARST's accuracy falls, and the more dispersed the population, the more it falls; with the true value, KARST always leads alfak."
- Saying "on Salehi, cell lines with high dispersion are more sensitive to the missegregation rate" is not recommended: the per-cell-line data do not support it.

## 8. Scale and cost

- **CPU:** 1,587 hours in total. Fits 1,431 hours (including the 48 GB reruns), scoring 141 hours, Fig. 2d noise 15 hours, collection 1 hour. More than the plan's estimate of 1,300–1,500 hours.
- **Wall clock:** 2026-11-23 21:06 to 11-24 04:17.
- **Memory:** fits default 16 GB; alfak all time points, minobs 5 needs up to 47 GB.
- **Jobs:** recorded in `results_6ae2a24/beta_grid/jobs.tsv`; every job had an add-hours monitor, cancelled after the job ended.

## 9. Files

**Code (not committed):**
- `R/wasserstein.R`: `wasserstein_safe()`, Wasserstein distance with a timeout fallback.
- `R/metrics_salehi.R`, `scripts/fig2d_check.R`: switched to `wasserstein_safe()`.
- `scripts/analyse_salehi.R`: additionally saves `sel_eval.Rds`.
- `scripts/beta_grid_compare.R`: comparisons and figures.
- `hpc/bg_job.sh`: job wrapper.

**Results:** `results_6ae2a24/beta_grid/`.
- On HPC: each β's fits and all analyses `x01`–`x64`, task tables `tasks/`, submit script `submit_bg.sh`, `jobs.tsv`, `logs/`.
- Locally: each β's analysis tables, and `analysis/`: `fig2c_paired.tsv`, `fig2c_rule.tsv`, `fig2d_rule.tsv`, `fig2d_noise_rho.tsv`, `fig2d_noise_differences.tsv`, `own_best.tsv`, `line_best.tsv`, `failures.tsv`, `figures/`.
