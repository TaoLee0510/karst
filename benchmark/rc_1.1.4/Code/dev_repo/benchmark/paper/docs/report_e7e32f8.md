# Comparison report of KARST e7e32f8 and alfak on the ALFA-K paper tests

Date: 2026-11-18. KARST e7e32f8 (1.0.0, dev, including the fixes from the 2026-11-17 code review). alfak dev2@f93b2d2 serves as the reference (fits linked from the alfak2 era); the paper's own fits are labelled alfak main@4ad82b0. Results directory on HPC: `/share/lab_crd/taoli/Project/karst/benchmark/paper/results_e7e32f8/`.

> Terms: the columns of the input are called "time points". A Salehi "trajectory" consists of several passages. alfak2 is KARST's name before 0.2.0; commit ids from the alfak2 era map to current commits via `doc/commit_map.tsv`.

## 0. What this round did

- After the rename to KARST, all KARST fits were redone: ABM, missegregation-rate simulations and Salehi, two designs (two time points; all time points, including pairwise and state space), 5617 tasks in total. In addition, 482 Salehi fits of KARST with `initial_state = "free"` were added, for the multi-horizon test.
- alfak is not refit: 3107 fit folders are linked from the alfak2-era round 8327282; these are the same directories read by the two rounds 05a8855 and 348addd.
- With the same image and the same KARST library, the multi-horizon prediction test and the development benchmark (including held-out, see Section 6) were also rerun.
- The previous round b9a275e (after the rename, before the review fixes) was stopped by the user's decision when the fits were mostly done but no tables had been produced yet, and is replaced by this round.

Among the review fixes, only one changes the results of this benchmark: when `as_alfak()` cross-validates a pairwise fit, the prior centre for unseen karyotypes in each pair of time points now uses the depth of that pair's starting point, consistent with the fit itself. The other fixes are not on this benchmark's computation path.

## 1. Conclusions

1. **The fits of KARST 1.0.0 are the same as those of alfak2 0.2.0.** One-by-one comparison with alfak2 05a8855 (`analysis/compare_rounds/`):
   - anchored (default): the status of all 8071 fits agrees; where both sides succeeded, the landscapes are identical within 1e-16, and the choice of history length, the cross-validation predictions and the ABM forecasts are exactly the same.
   - state space: the 3265 fits agree likewise.
   - pairwise: the landscapes, history lengths and forecasts of the 3265 fits are all the same; only the cross-validation changed (the fix above), affecting inputs with unequal depths.
   - Comparison with b9a275e: anchored and state space are exactly the same (b9a275e still lacked one Salehi fit when it was stopped); pairwise differs only in cross-validation.
2. **Tables that do not select by cross-validation are identical row by row to alfak2 348addd.** These include the paired landscape and forecast comparisons, the missegregation-rate simulations, the Salehi coverage and failure statistics, the Fig. 2d noise check, the Fig. 5 test, Fig. 6, and 109 of the 114 rows of the summary accuracy table (the other 5 rows are in Section 2). All rows of alfak and of the paper fits are the same.
3. **When fits are selected by CV score under the paper rule, ties are now handled by a fixed rule (Section 7).**
   - The paper rule keeps all minobs tied for the best CV score. KARST's minobs values are derived from the same fit, so their cross-validation predictions are often identical and the scores differ only by floating-point error; as a result one input is counted 1, 2 or 3 times, and the count changes with floating-point error. For KARST's 1921 inputs on ABM, the original rule kept 3282 rows; alfak never ties.
   - User decision on 2026-11-18: scores that differ by less than 1e-8 count as tied, the smallest minobs is taken, and only one is kept per input (`R/select.R`). For alfak and the paper's own fits, the selections are identical one by one to those of the original rule.
   - The tables before de-duplication are kept in `analysis_cv_ties_kept/` and `figures_cv_ties_kept/`. After de-duplication KARST has about 35–50% fewer samples (ABM); the pooled proportions change by about 0.03 or less (for example Fig. 1g at passage 10, window 2, the proportion beating no change 0.68 → 0.65), and finer subgroups change more (Section 7).
   - Compared with alfak2 348addd, the significance of none of the 19 conclusions of Fig. 4–6 changed; one state space conclusion (4f `treat_pairsame`, −0.0013 → +0.0009) changed sign near 0, and the number agreeing with the paper went from 10 to 11. The reason is that for one trajectory (SA609UnBU_X7_l_3) the CV scores of three minobs differ by only 1e-16, the two rounds each picked one of them, and these minobs list different karyotypes (minobs 5 lists 6831, 20 lists 2774). Now that the de-duplication rule always takes the smallest minobs, such a choice no longer changes with floating-point error.
4. **The main conclusions are the same as in the alfak2 era.**
   - On ABM and the missegregation-rate simulations, KARST's landscapes are more accurate than alfak's in every window, and better in 86–98% of paired inputs (Section 2).
   - Fig. 2d on Salehi mainly reflects sampling noise; with noise held fixed, no KARST arm differs reliably from alfak (Section 3).
   - When forecasting later passages, KARST beats alfak in both designs: two time points +1.01 per cell (lineage interval 0.67 to 1.40), all time points anchored +1.12 (0.55 to 1.47) (Section 4).

## 2. Paper tests: summary accuracy (`analysis/summary/summary_accuracy.tsv`)

Medians over paired inputs; the last column is the fraction of inputs in which KARST (two time points, or anchored for all time points) is better than alfak.

**ABM, landscape Pearson r**

| Design and window | Paired | KARST | pairwise | state space | alfak | KARST better |
|---|---|---|---|---|---|---|
| Two time points, window 2 | 1172 | 0.57 | | | 0.27 | 92% |
| Two time points, window 4 | 1173 | 0.68 | | | 0.42 | 97% |
| Two time points, window 8 | 1176 | 0.75 | | | 0.51 | 98% |
| All time points, window 4 | 1177 | 0.66 | 0.57 | 0.66 | 0.40 | 92% |
| All time points, window 8 | 1190 | 0.71 | 0.63 | 0.72 | 0.45 | 89% |

**ABM, landscape Spearman ρ**

| Design and window | KARST | pairwise | state space | alfak | KARST better |
|---|---|---|---|---|---|
| Two time points, window 2 | 0.57 | | | 0.28 | 94% |
| Two time points, window 4 | 0.68 | | | 0.42 | 95% |
| Two time points, window 8 | 0.74 | | | 0.50 | 98% |
| All time points, window 4 | 0.67 | 0.58 | 0.67 | 0.39 | 91% |
| All time points, window 8 | 0.71 | 0.63 | 0.72 | 0.43 | 88% |

**Missegregation-rate simulations, landscape Pearson r** (125 inputs, KARST uses the true missegregation rate)

| Design and window | KARST | pairwise | state space | alfak | KARST better |
|---|---|---|---|---|---|
| Two time points, window 2 | 0.67 | | | 0.39 | 95% |
| Two time points, window 4 | 0.78 | | | 0.61 | 90% |
| Two time points, window 8 | 0.87 | | | 0.79 | 90% |
| All time points, window 4 | 0.76 | 0.69 | 0.77 | 0.65 | 86% |
| All time points, window 8 | 0.81 | 0.74 | 0.83 | 0.80 | 66% |

**ABM, fraction in which the forecast beats no change at passage 10 (Euclidean distance)**

| Design and window | KARST | pairwise | state space | alfak |
|---|---|---|---|---|
| Two time points, window 2 | 0.70 | | | 0.36 |
| Two time points, window 4 | 0.85 | | | 0.45 |
| Two time points, window 8 | 0.86 | | | 0.40 |
| All time points, window 4 | 0.86 | 0.76 | 0.84 | 0.50 |
| All time points, window 8 | 0.91 | 0.78 | 0.89 | 0.53 |

**Salehi**

| Metric | Design | Paired | KARST | pairwise | state space | alfak |
|---|---|---|---|---|---|---|
| Forecast beats parent on the daughter (Wasserstein) | Two time points | 371 | 0.51 | | | 0.21 |
| | All time points (3–10 passages) | 283 | 0.54 | 0.52 | 0.49 | 0.27 |
| Fitness of novel karyotypes predicts their emergence | Two time points | 480 | 0.57 | | | 0.21 |
| | All time points | 351 | 0.52 | 0.66 | 0.55 | 0.24 |

These numbers are the same as in alfak2 348addd; only 5 rows have small differences, and none of them involve CV selection:
- In the Salehi novel-karyotype item, KARST two time points 0.565 → 0.567, all time points 0.521 → 0.516; in the same item, the paired proportions of pairwise and of state space against anchored differ by ≤ 0.005. The reason is that for one trajectory (SA039U_X51_l_3) the logistic regression on one daughter is completely separated, with coefficients reaching ±10¹⁴ and p = 0, so a 1e-16 difference in the landscape flips the sign (this is how the paper's 05 script does it, kept as is);
- ABM all time points window 4 Spearman: the fraction in which pairwise is better than anchored 0.322 → 0.323 (ties in the rank correlation).

**Failures and resources** (`analysis/summary/summary_performance.tsv`; alfak's resources are taken from the job records of round 8327282)

| Data | Design | Method | Inputs | Failure fraction | Median time per input (s) | Max (s) | Median memory (MB) | Max (MB) |
|---|---|---|---|---|---|---|---|---|
| ABM | Two time points | KARST | 1200 | 0 | 104 | 608 | 920 | 2025 |
| | | alfak | 1200 | 2.2% | 1228 | 5941 | 5436 | 23469 |
| | All time points | KARST | 800 | 0 | 281 | 2346 | 1815 | 6852 |
| | | pairwise | 800 | 0 | 412 | 4084 | 1928 | 6871 |
| | | state space | 800 | 0 | 837 | 20297 | 2404 | 8315 |
| | | alfak | 800 | 1.4% | 3525 | 25142 | 12577 | 88804 |
| Salehi | Two time points | KARST | 253 | 6.3% | 90 | 425 | 896 | 1725 |
| | | alfak | 253 | 6.7% | 529 | 2688 | 2359 | 10009 |
| | All time points | KARST | 187 | 2.7% | 237 | 1601 | 1585 | 4851 |
| | | state space | 187 | 2.7% | 715 | 31653 | 2351 | 13539 |
| | | alfak | 187 | 2.7% | 1586 | 12465 | 8003 | 41358 |

Time is the fit time needed to produce the outputs for all minobs, excluding forecasts; KARST fits each input only once. 24 Salehi trajectories were judged unanalysable before fitting and are excluded for all methods.

## 3. Fig. 2d: raw and with sampling noise held fixed

**The paper's rule**: fits selected by CV score, ties de-duplicated (`analysis/salehi/fig2d_wasserstein.tsv`). The selected subsets all differ and have few rows, so the arms cannot be compared directly.

| Fit | Rows | Spearman |
|---|---|---|
| KARST two time points | 50 | 0.23 |
| KARST all time points (anchored) | 49 | 0.30 |
| KARST pairwise | 50 | 0.56 |
| KARST state space | 48 | 0.54 |
| alfak dev2@f93b2d2 two time points | 44 | 0.45 |
| alfak dev2@f93b2d2 all time points | 46 | 0.70 |
| Paper (alfak main@4ad82b0) | 45 | 0.69 |

**Rows shared by all arms, no CV selection** (`analysis/fig2d_noise/`):
- "Noise" is the mean distance between the two halves after the parent and the daughter are each randomly split into two halves (20 splits), computed on the karyotypes of the alfak landscape.
- Intervals are bootstrap 95% intervals over the 6 lineage types.

| Design and subset | Rows | Noise itself | KARST raw | KARST noise fixed | alfak raw | alfak noise fixed |
|---|---|---|---|---|---|---|
| Two time points, all trajectories | 471 | 0.75 | 0.29 | 0.47 | 0.32 | 0.44 |
| Two time points, one interval (2 passages) | 147 | 0.77 | 0.45 | 0.49 | 0.60 | 0.55 |
| All time points, 3–10 passages: anchored | 346 | 0.73 | 0.23 | 0.40 | 0.55 | 0.52 |
| Same: pairwise | 346 | | 0.35 | 0.55 | | |
| Same: state space | 346 | | 0.45 | 0.56 | | |

Difference from alfak (KARST minus alfak):

| Design and subset | Arm | Raw | Noise fixed |
|---|---|---|---|
| Two time points, all trajectories | KARST | −0.03 (−0.19 to 0.15) | +0.03 (−0.16 to 0.21) |
| Two time points, one interval | KARST | −0.15 (−0.19 to −0.02) | −0.07 (−0.25 to 0.04) |
| All time points | anchored | −0.32 (−0.53 to 0.09) | −0.12 (−0.41 to 0.13) |
| All time points | pairwise | −0.20 (−0.37 to 0.10) | +0.03 (−0.10 to 0.14) |
| All time points | state space | −0.10 (−0.27 to 0.11) | +0.04 (−0.28 to 0.16) |

With noise held fixed, no arm differs reliably from alfak. The raw Fig. 2d differences come mainly from how close the size of the forecast is to the noise.

## 4. Forecasting later passages (`validation/sparse_joint/salehi_horizon`)

- This round's fits (each KARST arm and KARST free) forecast from the last passage of the input, using this benchmark's engine (the ABM of `alfakR::predict_evo`, 5 replicates, pm 5e-5, nb 1e7); every daughter passage within the next 4 passages is used as truth.
- alfak's forecast rows reuse the results of that alfak2-era run (same alfak fits, same engine); only the method name was changed.
- Score: log score per cell, relative to "no change" (the parent's frequencies), mixing weight 0.01, computed on the karyotypes of the fitted landscape.
- Intervals: bootstrap 95% intervals by lineage type.
- All tables of this round differ from that alfak2-era run by less than 1e-14.

**Two time points**: the 986 rows shared by 143 trajectories.

| Passages ahead | KARST vs "no change" | alfak vs "no change" | KARST minus alfak | Rows where KARST is higher |
|---|---|---|---|---|
| 1 (471 rows) | +0.47 (0.17 to 0.67) | −0.16 (−0.76 to 0.22) | +0.63 (0.34 to 0.99) | 65% |
| 2 (280 rows) | +0.27 (−0.44 to 0.79) | −1.10 (−2.17 to −0.20) | +1.37 (0.90 to 2.04) | 70% |
| 3–4 (235 rows) | −0.81 (−1.96 to 0.43) | −2.14 (−3.53 to −0.73) | +1.33 (1.10 to 1.63) | 71% |
| All | | | +1.01 (0.67 to 1.40) | 68% |

KARST free at 1 passage ahead is +0.44 (0.18 to 0.61), not higher than the default (balance).

**All time points**: 99 trajectories with 3 or more passages, the 677 rows shared by all four arms.

| Passages ahead | anchored | pairwise | state space | alfak |
|---|---|---|---|---|
| 1 (346 rows) | +0.78 (0.35 to 1.01) | +0.78 (0.31 to 1.02) | +0.77 (0.26 to 1.01) | +0.25 (−0.15 to 0.61) |
| 2 (191 rows) | +0.74 (−0.24 to 1.38) | +0.69 (−0.35 to 1.37) | +0.63 (−0.37 to 1.28) | −0.66 (−1.35 to −0.25) |
| 3–4 (140 rows) | −0.01 (−1.86 to 1.39) | −0.22 (−2.05 to 1.22) | −0.55 (−2.08 to 0.70) | −2.18 (−2.79 to −1.80) |

| Paired difference (all horizons) | Difference | Lineage interval | Rows where the former is higher |
|---|---|---|---|
| anchored minus alfak | +1.12 | 0.55 to 1.47 | 69% |
| pairwise minus alfak | +1.06 | 0.47 to 1.44 | 69% |
| state space minus alfak | +0.96 | 0.40 to 1.31 | 69% |
| state space minus anchored | −0.15 | −0.20 to −0.05 | 35% |
| pairwise minus anchored | −0.06 | −0.10 to −0.03 | 46% |
| free minus anchored | −0.15 | −0.24 to −0.01 | 32% |

All time points vs two time points, on the same trajectories (636 rows): anchored +0.36 (0.14 to 0.49, 90% of rows higher), alfak +0.01 (−0.71 to 0.91). The two designs choose the landscape karyotypes differently, so across designs only each design's gain relative to "no change" is compared.

## 5. One-by-one comparison with alfak2 (`analysis/compare_rounds/`)

| Method | Compared with | Fits (by minobs) | Same status | Same landscape (≤1e-16) | Same cross-validation predictions | Same ABM forecasts |
|---|---|---|---|---|---|---|
| anchored | alfak2 05a8855 | 8071 | all | all where both sides succeeded | all | all |
| state space | alfak2 05a8855 | 3265 | all | all | all | all |
| pairwise | alfak2 05a8855 | 3265 | all | all | ABM w4 222/1200, w8 144/1200, Salehi 18/556, missegregation rate 2/125 + 2/125 | all |
| anchored / state space / pairwise | KARST b9a275e | as above | all except the 1 Salehi fit missing from b9a275e | all | anchored, state space all; pairwise as above | all |

The pairwise fits whose cross-validation differs are all inputs whose depth differs between time points, which is exactly the case the fix targets. The largest differences in the ABM scores (including the cross-validation metrics and rank correlations): anchored 0.027, state space 0.0010, pairwise 0.095. Even when landscapes differ by only 1e-16, ties in the rank correlations can differ.

## 6. Development benchmark (`benchmark/sparse_two_timepoint`, including held-out)

- 1908 KARST fits were redone; the 811 fits of alfak (minobs 5, 20) and of the uninformative predictor were copied over. The distance 3–4 sweep did another 1920 KARST fits. One fit per task, with the same image and the same KARST library as this round.
- Compared with the development benchmark of round b9a275e, all metrics except run time are exactly the same, and the generated tex tables are byte-identical.
- On held-out, the 168 balance fits compared with the alfak2 dbc35fe rerun (2026-11-16, same model code as 0.2.0, run on a Mac): 146 have fitness differences < 1e-4, the largest is 0.021, and the coverage is exactly the same (0.977 / 0.989 / 0.982). The differences come from floating-point arithmetic on the Mac vs HPC.
- Compared with the earlier run of 2026-11-11, errors and rank correlations differ by ≤ 0.002; coverage at distances 1 and 2 is higher (held-out balance 0.96 / 0.94 → 0.99 / 0.98), because those fits predate the deficit term for "karyotypes not seen at any time point" added on 11-13 (kappa = 1.1).

| Distance | Uninformative | Two-layer pipeline | alfak (minobs 20) | sparse joint | balance (default) |
|---|---|---|---|---|---|
| Development 0 | 0.046 | 0.045 (0.30) | 1.073 (0.23) | 0.033 (0.44), 0.97 | 0.026 (0.47), 0.98 |
| Development 1 | 0.092 | 0.091 (0.24) | 1.186 (0.20) | 0.055 (0.43), 0.97 | 0.044 (0.45), 0.99 |
| Development 2 | 0.137 | 0.136 (0.15) | 1.297 (0.18) | 0.079 (0.36), 0.95 | 0.066 (0.38), 0.98 |
| held-out 0 | 0.044 | 0.043 (0.22) | 0.277 (0.16) | 0.032 (0.41), 0.97 | 0.025 (0.44), 0.98 |
| held-out 1 | 0.090 | 0.090 (0.16) | 0.371 (0.13) | 0.058 (0.36), 0.98 | 0.046 (0.37), 0.99 |
| held-out 2 | 0.136 | 0.135 (0.09) | 0.454 (0.12) | 0.083 (0.30), 0.96 | 0.068 (0.31), 0.98 |

RMSE, with Spearman in parentheses, followed by the coverage of the 95% interval; the eight core scenarios pooled, 96 datasets per cell. The numbers in the method document (`doc/karst_method.tex` Section 5 and the appendix, `doc/generated/`) are taken from this rerun.

## 7. Ties in CV scores (review item 8)

The user decided on 2026-11-18 to de-duplicate (see Section 1, item 3). What was rerun is the analyses that select by CV: `analyse_abm`, `analyse_salehi`, `summary_comparison`, Fig. 4 (80 tasks), Fig. 6 (8 tasks), `fig456_summary` (jobs 22041408–22041417), without refitting. Affected tables:

| Paper figure | Table | Change for KARST |
|---|---|---|
| Supp. Fig. 4c | `abm/supp4c_selected_by_cv.tsv` | medians of r, R², ρ of the selected fits change by ≤ 0.09 |
| Supp. Fig. 5a, 5b, Fig. 1g | `abm/supp5*`, `abm/fig1g_beat_baseline.tsv` | 35–50% fewer samples; the proportion beating no change changes by ≤ 0.03 (Supp. 5b), and by ≤ 0.07 in each cell of Fig. 1g split by λ and passage |
| Fig. 2b | `salehi/fig2b_*` | CV median 0.86 → 0.89 (two time points) |
| Fig. 2c | `salehi/fig2c_angle_test.tsv` | median angle for all time points anchored 55.3° → 66.7°, pairwise 51.7° → 64.5°; all still significant |
| Fig. 2d (paper rule) | `salehi/fig2d_wasserstein.tsv` | see Section 3, Spearman changes by ≤ 0.01 |
| Fig. 2e, 2f | `salehi/fig2e_*`, `salehi/fig2f_*` | 1–3 fewer trajectories per arm; pooled proportions change by ≤ 0.02, subgroups by lineage type by at most 0.25 |
| Fig. 3e–f (rows of the paper rule) | `salehi/fig3ef_novel_karyotypes.tsv` | 1–3 fewer regressions |
| Fig. 4, 6 | `fig456/` | unchanged within this round (these two were already de-duplicated by the smallest minobs; only ties were judged by strict equality) |

The two other items marked "possible" in the review (items 9 and 10), and the old-name mapping in the scripts that read alfak2-era results, are still unchanged.

## 8. Runs

Job ids are recorded in `build/e7e32f8_jobs.tsv` on HPC:
- build 22020580;
- fits 22020610–22020636 (11 groups plus KARST free and the two groups of the development benchmark), distance 3–4 sweep 22020668;
- scoring 22020638, collect_analyse 22020640, Fig. 4 22020641 / 643, Fig. 6 22020645 / 647, fig456_summary 22020649, Fig. 2d noise 22020650 / 652;
- comparisons 22020654–22020659;
- multi-horizon test 22020660 / 662, analysis 22040335 (the first attempt 22020664 failed because the node had no X11 to draw PNGs; rerun after switching to cairo) and 22020666;
- development benchmark scoring and report 22020670, document figures and tables 22020671.

Fits ran from 2026-11-17 21:09 to 2026-11-18 06:31, and all 9927 tasks completed; 5 NODE_FAILs were requeued automatically by Slurm. All jobs ran in qos xxlarge, 1 CPU per task, with an add-hours monitor for each array.

## 9. Files

- `analysis/compare_rounds/`: one-by-one comparisons with alfak2 05a8855 and KARST b9a275e.
- `analysis/summary/`, `analysis/abm/`, `analysis/salehi/`, `analysis/fig456/`: paper tests.
- `analysis/fig2d_noise/`: raw and noise-fixed versions of Fig. 2d.
- `validation/sparse_joint/salehi_horizon/results/` (inside the karst clone on HPC): rows and report of the multi-horizon test.
- `benchmark/results/sparse_two_timepoint/` (same): metrics, report and `generated/` of the development benchmark; `round_b9a275e/` holds the results of the stopped round.
