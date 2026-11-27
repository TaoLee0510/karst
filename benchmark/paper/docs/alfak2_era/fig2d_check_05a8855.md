# Fig. 2d ranking check: why alfak2 ranks poorly with two time points (round 05a8855)

Date: 2026-11-16. Only the existing fits and forecasts of round 05a8855 are read; nothing is refitted and alfak2 is not changed.

> Terminology: the input columns are called "time points". A Salehi "trajectory" consists of several passages; "a trajectory of two passages" means there is only one interval between the first and the last.

## Question

Fig. 2d compares two distances for each forecast at the next passage:

- predicted change: the Wasserstein distance from the parental distribution to the forecast distribution;
- observed change: the Wasserstein distance from the parental distribution to the daughter distribution.

The Spearman correlation of these two distances across lineages, in round 05a8855, is:

| | n | Spearman | Pearson |
| --- | --- | --- | --- |
| alfak2, two time points | 53 | 0.251 | 0.477 |
| alfakR, two time points | 44 | 0.447 | 0.759 |

The question this check has to answer is why alfak2's Spearman is low.

## Conclusions

1. **0.25 vs 0.45 is not a stable difference between the two methods.**
   - On all 471 rows where both methods have forecasts (same trajectory, same minobs, same daughter passage), alfak2 is 0.29 and alfakR is 0.32. Paired difference −0.03, 95% interval −0.19 to 0.15.
   - The paper's rule picks fits by each method's own CV score, which leaves only 53 and 44 rows, and these come from only 6 lineages. The different picked subsets, together with the small number of rows, are what produced 0.25 vs 0.45.
2. **The real difference appears only on trajectories of two passages**, that is, with only one interval between the first and the last, mostly 15 days.
   - Here alfak2 is 0.45 and alfakR is 0.60, a difference of −0.15. The 95% interval clustered by lineage is −0.19 to −0.02, and by source line −0.29 to −0.01. On trajectories with 3 or more passages, the two are level, or alfak2 is slightly better.
   - The cause is the history length that alfak2 chooses. This history length determines the size of the predicted change: on the log scale the slope is −0.70, so a 10 times longer history gives a predicted change of about 1/5. But it is almost unrelated to the observed change (0.15).
   - After removing the effect of history length, alfak2's Spearman is 0.61, level with alfakR's 0.60 (alfakR after the same removal is 0.63).
   - In the subset that alfak2 picks by CV, such trajectories make up a larger share: 21 of 53 rows, against 10 of 44 rows for alfakR.
3. **None of the other three hypotheses holds.**
   - **H2 karyotypes used for scoring:** the landscape karyotypes of the two methods are exactly the same, and the observed axis agrees on all 471 rows.
   - **H3 amount of data or missegregation rate:** the predicted change of both methods comes almost entirely from selection. The predicted magnitude is unrelated to the amount of data.
   - **H4 properties of the model itself:** on ABM simulations with known truth, alfak2 ranks better than alfakR (0.82–0.87 vs 0.75–0.76; the upper bound using the true landscape is 0.89). Its predicted fitness variance is also closer to the truth.
4. **Both methods reach only about 0.3 on Salehi, because a large part of the observed change is unrelated to selection.**
   - The observed change correlates with parental diversity (0.45), and negatively with the fraction of parental cells inside the landscape (−0.54).
   - After removing diversity and karyotype count, both rise to about 0.45.

## H1: picking

Tables `s1_selection.tsv`, `s1_paired_difference.tsv`, `s1_by_length.tsv`.

The paired difference is alfak2 − alfakR. The intervals are clustered by two units: 6 lineages (the unit of the paper's Fig. 2c) and 12 source lines (such as SA609, SA906b).

| Row set | n | alfak2 | alfakR | Difference | 95% interval, by lineage | By source line |
| --- | --- | --- | --- | --- | --- | --- |
| Paper rule, each method picks its own | 53 / 44 | 0.251 | 0.446 | rows differ, cannot be paired | | |
| Common rows, all minobs | 471 | 0.294 | 0.323 | −0.03 | −0.19, 0.15 | −0.20, 0.13 |
| Common rows, minobs 5 / 10 / 20 | 177 / 163 / 131 | 0.29 / 0.29 / 0.25 | 0.36 / 0.18 / 0.24 | −0.07 / +0.10 / +0.02 | all contain 0 | |
| alfak2's picking | 53 | 0.251 | 0.397 | −0.15 | −0.38, 0.33 | −0.45, 0.20 |
| alfakR's picking | 44 | 0.269 | 0.446 | −0.18 | −0.48, 0.03 | −0.41, −0.01 |

On the picked subsets alfakR is higher by about 0.15–0.18, but only under one of the two clusterings does the interval narrowly exclude 0.

**By trajectory length** (common rows):

| Passages | n | alfak2 | alfakR | Difference | 95% interval, by lineage | By source line |
| --- | --- | --- | --- | --- | --- | --- |
| 2 (one interval) | 147 | 0.445 | 0.597 | −0.15 | −0.19, −0.02 | −0.29, −0.01 |
| 3 | 107 | 0.491 | 0.447 | +0.04 | −0.20, 0.44 | −0.20, 0.41 |
| 4 | 91 | 0.193 | 0.118 | +0.08 | −0.23, 0.43 | −0.21, 0.36 |
| 5 or more | 126 | 0.043 | −0.050 | +0.09 | −0.28, 0.50 | −0.22, 0.36 |

- **For two passages, the difference excludes 0 under both clusterings.**
- **The longer the trajectory, the worse both methods are.** The first and last passages are far apart, while the forecast only looks at the next passage.

## H2: which karyotypes are used on the axes

Tables `s2_axes.tsv`, `s2_by_coverage.tsv`, `s2_coverage_medians.tsv`.

- **The axes of the two methods are the same.** The landscape that alfak2 outputs has the same karyotype set as alfakR's, and the observed axes of the 471 rows all agree. The plan's original assumption, "the two sides have different axes", does not hold.
- **Coverage.** The median fraction of parental cells inside the landscape is 0.79, and inside the alfak2 model domain 0.69.
- **Grouping by coverage shows no consistent pattern.** Almost every group comes from one or two lineages.
- **Other axes:**

  | Axis | alfak2 | alfakR |
  | --- | --- | --- |
  | The paper's axis (karyotypes in the landscape) | 0.294 | 0.323 |
  | All observed karyotypes, cells outside the landscape kept unchanged in the forecast | −0.274 | −0.216 |
  | Euclidean distance between mean karyotypes | 0.359 | 0.361 |

  With all karyotypes, both become negatively correlated: no landscape can predict the change brought by cells outside the landscape.

## H3: what determines the size of the predicted change

Tables `s3_spearman_features.tsv`, `s3_partial.tsv`, `s3_scale.tsv`.

- **The predicted change comes almost entirely from selection.** The change computed from selection alone, i.e. parental frequency multiplied by exp(fitness × days), has a Spearman of 1.000 (alfak2) and 0.996 (alfakR) with the actual forecast, and a median ratio of 0.99. Missegregation diffusion and randomness in the forecast are negligible.

  The plan's step "re-forecast with the missegregation rate set to 0" is therefore replaced by this deterministic calculation; no new simulations are needed.
- **The predicted change follows the selection among parental cells.** Its Spearman with the Price vector (covariance of fitness and copy number) is 0.94 and 0.91, and with the fitness variance 0.79 and 0.81.
- **Different scales, similar ranking.**
  - alfak2's median fitness variance is 0.32 times alfakR's, and its predicted change 0.61 times.
  - The rank correlation between the two methods' predicted changes is 0.73.
  - Spearman is not affected by the overall scale.
- **Unrelated to the amount of data.**
  - The correlation of the variance ratio (alfak2/alfakR) with the number of cells is 0.05, and with the number of major karyotypes −0.15.
  - The correlation of alfak2's predicted change with the number of parental cells is 0.19, and with the number of karyotypes in the landscape 0.03.
- **The observed change itself (the same for both methods):**
  - correlation 0.45 with parental diversity (entropy);
  - correlation −0.54 with the fraction of parental cells inside the landscape: the fewer cells inside the landscape, the noisier the normalised distribution;
  - correlation −0.51 with the coverage of the alfak2 model domain.

  After removing diversity and karyotype count, the partial correlation is 0.46 for alfak2 and 0.44 for alfakR.

## Exploration outside the plan: the history length chosen by alfak2

Tables `s3_history_by_length.tsv`, `s3_history_groups_length2.tsv`, figure `figures/salehi_fig2d_length2_by_history.png`.

On trajectories of two passages, alfak2 chooses three kinds of history: equal to the interval, 2 intervals, or 3.3–12 intervals (mostly 30 days and 100 days).

| Passages | History/interval vs alfak2 predicted change | vs observed change | log slope | alfak2 | alfak2, history removed | alfakR | alfakR, history removed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | −0.50 | 0.15 | −0.70 | 0.445 | 0.608 | 0.597 | 0.632 |
| 3 | −0.18 | −0.14 | −0.28 | 0.491 | 0.479 | 0.447 | 0.436 |
| 4 | −0.07 | 0.06 | −0.34 | 0.193 | 0.198 | 0.118 | 0.118 |
| 5 or more | 0.17 | −0.13 | 0.49 | 0.043 | 0.067 | −0.050 | −0.027 |

Trajectories of two passages, grouped by alfak2's history:

| History | n | alfak2 | alfakR | Median predicted change, alfak2 / alfakR | Fitness variance ratio |
| --- | --- | --- | --- | --- | --- |
| Equal to the interval | 45 | 0.69 | 0.67 | 0.67 / 1.10 | 0.49 |
| 2 intervals | 54 | 0.64 | 0.63 | 0.42 / 1.11 | 0.14 |
| More than 2 intervals | 48 | 0.41 | 0.55 | 0.22 / 0.62 | 0.08 |

**Interpretation:**
- The longer the history, the more alfak2 treats the parental composition as the result of selection over a longer time, so it gives smaller fitness differences and a smaller predicted change.
- On data with only one interval, the history choice is unrelated to how much actually changed in that period, so it adds noise to the predicted magnitude. Within each history group, alfak2 and alfakR are level.
- This is the same problem as the "history" item (A) of the short-interval weakness in the 2026-11-12 handoff document.

**The ABM does not show this effect** (table `s4_abm_history.tsv`):
- w2_ends also has only one interval (about 7 days), and alfak2's history is also widely spread: equal to the interval in only 17%, the rest 30 to 3000 days.
- But the correlation between history and predicted change is only −0.05.
- Why Salehi differs was not investigated this time.

## H4: ABM simulations with known truth

Tables `s4_abm.tsv`, `s4_abm_pooled_minobs5.tsv`. From the forecast start to the next sampling time point (about 7 days), 400 simulations, minobs 5:

| Design | Method | Simulated version of Fig. 2d | Upper bound with the true landscape | Fitness variance, ranking against truth | Fitness variance / truth (median) |
| --- | --- | --- | --- | --- | --- |
| w2_ends | alfak2 | 0.821 | 0.894 | 0.774 | 0.80 |
| w2_ends | alfakR | 0.759 | 0.892 | 0.251 | 2.13 |
| w4_ends | alfak2 | 0.871 | 0.894 | 0.862 | 0.85 |
| w4_ends | alfakR | 0.750 | 0.892 | 0.442 | 1.77 |
| w8_ends | alfak2 | 0.842 | 0.894 | 0.814 | 0.78 |
| w8_ends | alfakR | 0.752 | 0.888 | 0.544 | 1.75 |

- **The same for every wavelength and every design.** alfak2's ranking is closer to the upper bound.
- **Scale of the fitness variance:** alfak2 is slightly smaller (about 0.8 times), alfakR too large (about 1.8–2.1 times).
- **Median angle between the Price vector and the truth:** alfak2 8–20°, alfakR 15–36°.

So ranking the magnitude of the predicted change is not in itself a weakness of alfak2.

## Implications for alfak2

- There is no need to change alfak2's defaults for Fig. 2d. The overall ranking is level with alfakR, and better on simulations.
- The only place that could be improved is the history choice with a single short interval. This needs a separate plan, for example one of the following, followed by validation on simulated data and on Salehi:
  - fix the history to the interval when there is only one interval;
  - average over the history candidates instead of taking only one;
  - choose a longer history only when the evidence is sufficient.
- This check did none of these.

## Methods and files

- **Scripts** (new files; existing scripts not changed):
  - `scripts/fig2d_check.R`: computes on HPC one fit directory at a time, `--stage=salehi|abm --index=N`;
  - `scripts/fig2d_check_analyse.R`: local analysis;
  - `hpc/fig2d_check_task.sh`: array task.
- **Checks:**
  - The two axes recomputed by this check are identical, row by row, to the values in the benchmark table: all 97 picked rows were found.
  - The Spearman values under the paper's rule are reproduced exactly: 0.251 / 0.477, 0.447 / 0.759.
- **Deviations from the plan:**
  - H2 was changed to checking coverage and alternative axes, because the karyotype sets of the two sides are the same;
  - H3 uses the deterministic "selection only" change instead of re-forecasting.
- **HPC** (qos xxlarge, 1 CPU, 8 GB, 4 hours per task):

  | Job | Content | Result |
  | --- | --- | --- |
  | 21969835 | Salehi test | done |
  | 21969836 | ABM test | failed: the simulated counts are a data frame; script fixed |
  | 21969982 | ABM test | done |
  | 21969984 | Salehi, 27 tasks | all done, about 1 minute each |
  | 21970013 | ABM, 99 tasks | all done, about 1.5 minutes each |

  - Add-hours monitors: 21969846 and 21969847 had started on small before the switch to xxlarge; 21969985, 21970012 and 21970014 ran on xxlarge. Each was cancelled after its job ended.
- **Results:**
  - On HPC in `results_05a8855/analysis/fig2d_check/`, with a local copy: the per-block rows are in `salehi/` and `abm/`, the tables are `s*_*.tsv`, and the figures are in `figures/`.
  - The tables and metadata used for the local analysis were copied from HPC: `results_05a8855/tables/salehi_scores.Rds`, `salehi_inputs.Rds`, and `results_05a8855/paper_meta/metadata.csv`.
