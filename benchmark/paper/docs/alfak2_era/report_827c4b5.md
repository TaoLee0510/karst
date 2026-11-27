# alfak2 827c4b5: automatic choice of the initial state with two time points (record of a regression)

Date: 2026-11-17. alfak2 827c4b5 (dev), alfakR dev2 f93b2d2. Results directory `results_827c4b5/`.

> **This round is not an official result.** It records a regression: the default setting of 827c4b5 lowered landscape accuracy on windows of one passage, and the user therefore changed the default back to `"balance"` in 348addd. For the official results see `docs/report_348addd.md` (the same as 05a8855).
>
> Terms: the columns of the input are called "time points". A Salehi "trajectory" consists of several passages.

## 0. What this round tests

From 827c4b5 on, with two time points `fit_sparse_joint()` defaults to `initial_state = "auto"`:
1. First a check fit is done on the interval: balance prior, the overall level of the non-major karyotypes kept as one free log factor (`balance_level = "free"`), and the history set to the interval itself.
2. When this log level is above `auto_level = 1`, `"free"` is used instead (the initial abundance of each karyotype estimated separately); otherwise `"balance"` is kept.
3. With more than two time points, `"auto"` is simply `"balance"`.

`run.Rds` records `initial_state` and `initial_state_check` (level, standard deviation, threshold, choice, whether the check fit converged).

**How this round was run**
- All of the alfak2 default arm refitted: 3107 tasks, two time points and all time points (anchored).
- pairwise and state_space each refitted on 44 check inputs. After they were found identical to 05a8855, the rest were linked from `results_05a8855/`.
- alfakR linked from `results/`.

## 1. Conclusions

1. **Switching to free happens almost only on windows of one passage (w2) and on Salehi PDX trajectories** (Section 2).
   - ABM w2: 54% (97% at λ 0.8). w4, w8: 0%.
   - misseg w2: 95%. These fits use the true missegregation rate.
   - Salehi: PDX 72%, in vitro 0%.
2. **Fits that switched to free have clearly lower landscape accuracy** (Section 4). Paired with the balance fit of the same input in 05a8855, both against the truth:
   - ABM w2: median landscape r 0.58 vs 0.75, 93% of fits worse; cross-validation held-out prediction against truth 0.47 vs 0.66.
   - misseg w2: worse at every true missegregation rate. At 5e-5 (missegregation rate set correctly, where it should not switch), 0.53 vs 0.72.
3. **Cause**: the check fit sets the history to the interval itself. When the window is only one passage, the "cloud" the population has accumulated over its whole history is much larger than what one passage can accumulate. So the level mainly measures "the true history is much longer than the interval", not a wrongly set missegregation rate.
   - The fraction switched is almost the same at every true missegregation rate (31%–33%, almost all in w2).
   - On the same simulations, the level is low in w4 and w8.
4. **Forecast scores rise slightly, but not enough to offset this.**
   - ABM w2 passage 10 beating the baseline: Euclidean 0.71 → 0.75, Wasserstein 0.70 → 0.74.
   - Salehi two time points: forecast beating the parent (Wasserstein) 0.51 → 0.55; Fig. 2d (paper rule) 0.25 → 0.58.
   - But Salehi "fitness predicts the appearance of new karyotypes" 0.57 → 0.38.
   - Fig. 4–6 conclusions agreeing with the paper: alfak2 two time points 11 → 12 (of 19). Fig. 4d PDX_v is still in the direction opposite to the paper (−1.15 → −0.80).
5. **Decision** (user):
   - the default is changed back to `"balance"` (348addd); `"auto"` is kept as an option, and the documentation states its limitation on intervals of one passage.
   - This round serves only as the record of this regression.

## 2. Choice (`analysis/initial_state/choice.tsv`, `tables/initial_state.tsv`)

Two-time-point inputs, grouped by source and design:

| Source | Design | Inputs | Switched to free | Median level (5% / 95%) | Check fit not converged |
|---|---|---|---|---|---|
| ABM | w2_ends | 400 | 54% | 1.07 (−0.33 / 1.61) | 47 |
| ABM | w4_ends | 400 | 0% | 0.09 (−0.64 / 0.57) | 38 |
| ABM | w8_ends | 400 | 0% | −0.82 (−1.66 / −0.04) | 84 |
| misseg | w2_ends | 125 | 95% | 1.38 (1.03 / 1.79) | 26 |
| misseg | w4_ends | 125 | 2 | 0.31 (−0.19 / 0.76) | 23 |
| misseg | w8_ends | 125 | 1 | −0.72 (−1.31 / −0.07) | 26 |
| Salehi | ends | 253 | 57% | 1.21 (−1.43 / 2.78) | 28 |

**ABM w2, by λ**
- λ 0.2: 0% (median level −0.04);
- λ 0.4: 31% (0.87);
- λ 0.8: 97% (1.32);
- λ 1.6: 86% (1.26).

**misseg w2, by true missegregation rate**: 5e-5 92%, 1e-4 100%, 2e-4 100%, 4e-4 100%, 8e-4 84%.

**Salehi**
- In vitro (the two lines SA039, SA906): 55 inputs, 0%; PDX: 198, 72%.
- PDX by number of passages: 2 passages 90%, 3 passages 72%, 4–5 passages 64%, 6 or more 64%.
- The interval used by the fit (days between the first and last passages):
  - in vitro median 100 days (10% / 90%: 25 / 175); PDX 45 days (15 / 90).
  - Within PDX, the median is 30 days for those switched to free and 60 days for those not switched.
  - The Spearman between level and interval is −0.63.
- So in vitro not switching and PDX switching is at least partly caused by the interval, and cannot be read as a biological difference between the two kinds of samples.

## 3. Comparison with 05a8855 (`analysis/compare_rounds/`)

**alfak2 default arm**
- All all-time-point fits are identical to 05a8855, landscape difference ≤ 2.2e-16, and cross-validation and ABM forecasts agree: ABM w4_all / w8_all 1200 each, misseg 125 each, Salehi all 546 (counted by minobs).
- The two-time-point fits that still use balance are also all identical: history, landscape, forecast all the same.
- Only the fits that switched to free changed: Salehi 429, ABM w2 642, misseg w2 119, misseg w4/w8 3.
- No fit changed its status (success or failure).

**pairwise / state_space**: the 44 + 44 check inputs are identical to 05a8855, landscape difference ≤ 3.3e-16; the rest linked.

**Differences in scoring**: for fits still using balance, scores differ by up to 0.03. 05a8855 against 8512694 had a difference of the same size (0.039); it is a scoring fluctuation that always exists between repeated runs, not a change of this round.

## 4. Against the truth: fits switched to free and their balance fits (`analysis/initial_state/vs_truth.tsv`)

`scripts/initial_state_vs_truth.R` pairs each fit switched to free (827c4b5) with the balance fit of the same input, same design and minobs in 05a8855, and scores both against the truth. The table gives medians.

**ABM w2**

| Group | Inputs (fits) | Landscape r | One-step neighbours rnn | Cross-validation held-out prediction against truth |
|---|---|---|---|---|
| All | 214 (642) | 0.58 vs 0.75 (93% worse) | 0.62 vs 0.77 | 0.47 vs 0.66 (83% worse) |
| λ 0.4 | 31 | 0.25 vs 0.42 | 0.28 vs 0.48 | 0.09 vs 0.30 |
| λ 0.8 | 97 | 0.62 vs 0.75 | 0.67 vs 0.78 | 0.56 vs 0.69 |
| λ 1.6 | 86 | 0.62 vs 0.80 | 0.65 vs 0.80 | 0.44 vs 0.67 |

Forecasts of ABM w2 (the scores of Fig. 1g):
- passage 10 beating the baseline: Euclidean 14% better, 6% worse; Wasserstein 13% vs 7%;
- average win over passages 1–10: 23% better, 12% worse;
- angle at passage 10: 37.6° vs 43.2°. Better at λ 1.6 (p = 8e-5); no difference at λ 0.4 and 0.8.

**misseg w2** (fits use the true missegregation rate)

| True missegregation rate | Inputs switched | Landscape r | Cross-validation held-out prediction against truth |
|---|---|---|---|
| 5e-5 | 23 / 25 | 0.53 vs 0.72 (91% worse) | 0.40 vs 0.62 |
| 1e-4 | 25 / 25 | 0.63 vs 0.74 | 0.58 vs 0.71 |
| 2e-4 | 25 / 25 | 0.61 vs 0.70 | 0.49 vs 0.72 |
| 4e-4 | 25 / 25 | 0.37 vs 0.60 | 0.26 vs 0.63 |
| 8e-4 | 21 / 25 | 0.34 vs 0.54 | 0.09 vs 0.45 |

p ≤ 1e-3 in every group. Inputs still using balance (ABM 186, misseg 6) have exactly the same scores on every measure.

**Convergence of the check fit**
- ABM w2: 17% of the 214 switched did not converge, 6% of the 186 not switched.
- misseg w2: 21% of the 119 switched.
- The level of switched inputs is unrelated to whether the check fit converged: the median is 1.29 for both in ABM, 1.37–1.40 in misseg. So switching is not caused by non-convergence.
- The final fits converge more often instead: in the 642 of ABM, non-converged went down from 126 to 48.

## 5. Changes in the paper tests (`analysis/summary/summary_accuracy.tsv`)

Compared with 05a8855, only 10 of the 114 accuracy rows changed, all alfak2 two time points (plus one row for new karyotypes of Salehi all time points, 0.521 → 0.519):

| Test | 05a8855 | 827c4b5 |
|---|---|---|
| ABM w2 landscape Pearson r (1172 pairs) | 0.57 | 0.42 |
| ABM w2 landscape Spearman rho | 0.57 | 0.44 |
| ABM w2 rescaled R² | 0.26 | 0.14 |
| misseg w2 landscape r (125 pairs) | 0.67 | 0.47 |
| ABM w2 passage 10 beating the baseline, Euclidean | 0.71 | 0.75 |
| ABM w2 passage 10 beating the baseline, Wasserstein | 0.70 | 0.74 |
| Salehi two time points, forecast beating the parent, Wasserstein (371 pairs) | 0.51 | 0.55 |
| Salehi two time points, forecast beating the parent, Euclidean | 0.44 | 0.44 |
| Salehi two time points, fitness predicting the appearance of new karyotypes (480 pairs) | 0.57 | 0.38 |

**Salehi inputs switched to free** (`analysis/initial_state/effect.tsv`): forecast beating the parent, Wasserstein 0.63 vs 0.54 (p = 0.004), Euclidean 0.43 vs 0.44.

**Fig. 2d**
- Paper rule (`analysis/salehi/fig2d_wasserstein.tsv`): alfak2 two time points 0.25 → 0.58 (55 fits).
- Shared rows (`analysis/fig2d_noise/`):
  - one interval: raw 0.45 → 0.59, level with alfakR's 0.60; after fixing the sampling noise, the difference from alfakR −0.07 → −0.04 (−0.20 to 0.11).
  - All trajectories: raw 0.29 → 0.37 (alfakR 0.32).
- Fig. 2d mainly reflects sampling noise (`docs/report_348addd.md` Section 4), so the gain here does not mean the landscape is more accurate. On simulations with truth, the landscape got worse instead.

**Fig. 4–6** (`analysis/fig456/conclusions.tsv`)
- alfak2 two time points, conclusions agreeing with the paper 12 / 19 (05a8855: 11). The extra one is `log1p(dk):pair_typeparallel` of Fig. 4f.
- Fig. 4d PDX_v: −1.15 → −0.80, still significant, direction opposite to the paper.
- All other arms are the same as 05a8855.
- In Fig. 6, the fraction of "landscapes whose dominant karyotype changes with the missegregation rate" changes slightly for every alfak2 arm (for example two time points 0.60 → 0.40), but the conclusions agree.

## 6. Runs

Job numbers are recorded in `build/827c4b5_jobs.tsv` on HPC:
- build 21975233;
- five alfak2 arrays 21975243 / 346 / 348 / 350 / 352;
- check arrays 21975354 / 356;
- scoring 21980388;
- collect_analyse 21992955, Fig. 4 21992958 / 962, Fig. 6 21992964 / 966, fig456_summary 21992968.

Problems during the run:
- **NODE_FAIL**: 28 tasks. Slurm requeues NODE_FAIL tasks automatically; I resubmitted 24 of them as well, and for 13 ABM tasks among them two copies wrote to the same fit folder at the same time. These 13 were deleted and each refitted once. From now on check `sacct -D` first, and do not resubmit by hand.
- **Links**: the links were made while 3 state_space check tasks were still running, so their locations were linked to 05a8855. These 3 links were found and deleted within two minutes; the files of 05a8855 were not changed that day. From now on, link only after all check tasks have finished.
- **collect_analyse failed the first time** (21980728): free fits have only one history and no `history_profile`, and `alfak2_chosen_history()` failed when the number of candidates was NA. It now records it as "dt" (`R/coverage.R`); the tables of earlier rounds are not affected.

## 7. Files

- `tables/initial_state.tsv`, `analysis/initial_state/choice.tsv`: choice and level.
- `analysis/initial_state/effect.tsv`: scores paired with 05a8855, grouped by choice.
- `analysis/initial_state/vs_truth.tsv`, `vs_truth_fits.tsv`: against the truth.
- `analysis/compare_rounds/`: comparison with 05a8855.
- `analysis/fig2d_noise/`: Fig. 2d raw and after fixing the noise.
- `analysis/summary/`, `analysis/abm/`, `analysis/salehi/`, `analysis/fig456/`: paper tests.
