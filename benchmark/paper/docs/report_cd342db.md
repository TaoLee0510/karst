# KARST cd342db versus alfak on the ALFA-K paper tests: comparison report

Date: 2026-11-20. KARST cd342db (1.1.2, dev). alfak dev2@f93b2d2 serves as the reference, with its fits linked from the alfak2 era; the paper's own fits are labelled alfak main@4ad82b0. Results directory on HPC: `/share/lab_crd/taoli/Project/karst/benchmark/paper/results_cd342db/`.

> Terms: the input columns are called "time points". A Salehi "trajectory" consists of several passages. "Fits (by minobs)" means one row per minobs; the minobs values of one input share a single fit.

## 0. What this round did

The previous report was e7e32f8 (1.0.0). Since then KARST has gone through three versions, all of which affect this benchmark:

- **1.1.0** (2c60af1, fixes from the 2026-11-18 code review):
  - `correct_efflux` converts fitness into birth rate with alfak's formula: `(f − k) / (2 P_ii − 1)`. Here `k` makes the frequency-weighted mean over the frequent karyotypes at time 0 equal to `log(nb / n0) / dt`. Every alfak-format output of KARST changes as a result, and so do the forecasts, the scores, the CV selection and Fig. 2–6. `P_ii` differs between karyotypes, so the landscape Pearson r also changes slightly.
  - Two fixes to multi-time-point fits: how karyotypes seen only at intermediate time points are handled; the drift weights in pairwise and state space cross-validation.
  - When no posterior mode is found on the fit's scale, an error is raised (`karst_inner_failure`) instead of carrying on with the last point.
  - Changes to the benchmark itself:
    - Salehi fit failure is decided by run status. The number of two-time-point pairs therefore changes from 371 to 394.
    - Fig. 5b uses only one minobs per trajectory.
    - The comparison with the paper's conclusions now uses 16 per arm.
- **1.1.1** (50415cf): when the search for the posterior mode fails, search once more from the initial latent-variable values.
- **1.1.2** (cd342db): when comparing candidate histories, use the objective value at the posterior mode reported by each candidate's fit, no longer the nlminb objective value. Only candidates with a finite nlminb objective enter the comparison.

The work in this round:

- All KARST fits redone: ABM, missegregation-rate simulations and Salehi, both designs (two time points; all time points, including pairwise and state space), 5617 tasks in total. In addition, 482 Salehi fits with KARST `initial_state = "free"`.
- alfak is not refit: 3107 fit folders are linked from the alfak2-era round 8327282, the same set of directories that round e7e32f8 read.
- The multi-horizon prediction test is rerun for all methods, including alfak's forecasts (using the same set of alfak fits). This is the first time the two designs are scored on a common karyotype support set.
- The development benchmark is not in this round. The effect of 1.1.2 on it was already measured in `validation/sparse_joint/history_choice`: the history changed in 5 of 168 datasets, and the median RMSE ratio is 1.000.

The two previous rounds:
- 2c60af1 (1.1.0) ran to completion, but 85 fit tasks failed: for some candidate history that was not selected, no posterior mode was found, so the whole fit raised an error. This round got no separate report.
- Shortly after 50415cf (1.1.1) was submitted, the user decided to change the history choice rule first, so this round was stopped and `results_50415cf/` has only partial results.

## 1. Conclusions

1. **This round is complete.**
   - All 6099 fit tasks ran to completion.
   - 20 tasks hit a node failure (NODE_FAIL) and finished after Slurm requeued them automatically.
   - 3 tasks failed at container start and were rerun (Section 8).
   - No timeouts and no out-of-memory.
2. **Fit success and failure are exactly the same as in e7e32f8.**
   - Only 112 fits failed, all Salehi `karst_domain_coverage_error`: one time point has no cell inside the model's domain. By method: anchored 42, free 42, pairwise 14, state space 14.
   - All 84 fits that failed in 2c60af1 because no posterior mode was found are recovered. Counted by minobs: anchored 129, pairwise 54, state space 38 (Section 6).
3. **1.1.1 and 1.1.2 affect only very few fits, and accuracy shows no consistent change (Section 6).**
   - The 5547 successful fits have 30099 candidate histories in total. Of these, 72 candidates (in 71 fits) needed the second search to find the mode, and only 1 fit selected such a candidate.
   - 6 candidates failed both searches, all in Salehi pairwise fits. None of them was selected, and their fits were completed with other candidates.
   - In 52 fits (0.9%) the new rule chose differently from the choice by nlminb objective.
   - Of these, 44 inputs from ABM and the missegregation-rate simulations were also fit successfully in 2c60af1 and can be compared in pairs: the landscape Pearson r got better in 23 and worse in 21, median change 0.000 (p = 0.57).
4. **The main conclusions are the same as in e7e32f8.**
   - ABM and missegregation-rate simulations: KARST's landscapes are more accurate than alfak's in every window, better in 66–98% of paired inputs (Section 2).
   - Salehi: Fig. 2d mainly reflects sampling noise; after fixing the noise, no KARST arm differs reliably from alfak (Section 3).
   - When forecasting later passages, KARST beats alfak in both designs (Section 4):
     - two time points: +1.04 per cell, lineage interval 0.71 to 1.43;
     - all time points, anchored: +1.15, lineage interval 0.59 to 1.51.
5. **New result: fitting with all time points forecasts better than using only two time points (Section 4).**
   - Scored on the common support set, on the same set of trajectories anchored is 0.51 higher per cell, lineage interval 0.13 to 0.79, higher in 89% of rows; higher in every days bin.
   - alfak gives +0.14 (−0.45 to 0.85).
   - Previously this comparison had only been computed on each design's own landscape, so the multi-horizon test's README marked it as "not established".

## 2. Paper tests: summary accuracy (`analysis/summary/summary_accuracy.tsv`)

Each table takes medians over the paired inputs.

- Column "KARST better": among paired inputs, the fraction where KARST (two time points, or anchored with all time points) is better than alfak.
- Column "better than anchored": among paired inputs, the fraction where pairwise or state space is better than anchored.

**ABM, landscape Pearson r**

| Design and window | Paired | KARST | pairwise | state space | alfak | KARST better | Better than anchored (pairwise / state space) |
|---|---|---|---|---|---|---|---|
| Two time points, window 2 | 1172 | 0.57 | | | 0.27 | 92% | |
| Two time points, window 4 | 1173 | 0.68 | | | 0.42 | 97% | |
| Two time points, window 8 | 1176 | 0.75 | | | 0.51 | 98% | |
| All time points, window 4 | 1177 | 0.66 | 0.57 | 0.66 | 0.40 | 92% | 32% / 65% |
| All time points, window 8 | 1190 | 0.71 | 0.63 | 0.72 | 0.45 | 89% | 30% / 64% |

**ABM, landscape Spearman ρ**

| Design and window | KARST | pairwise | state space | alfak | KARST better |
|---|---|---|---|---|---|
| Two time points, window 2 | 0.57 | | | 0.28 | 94% |
| Two time points, window 4 | 0.68 | | | 0.42 | 96% |
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

**ABM, fraction where the forecast beats "no change" at passage 10**

| Design and window | Euclidean: KARST | pairwise | state space | alfak | Wasserstein: KARST | pairwise | state space | alfak |
|---|---|---|---|---|---|---|---|---|
| Two time points, window 2 | 0.70 | | | 0.36 | 0.69 | | | 0.32 |
| Two time points, window 4 | 0.85 | | | 0.45 | 0.87 | | | 0.44 |
| Two time points, window 8 | 0.86 | | | 0.40 | 0.88 | | | 0.35 |
| All time points, window 4 | 0.86 | 0.77 | 0.84 | 0.50 | 0.87 | 0.71 | 0.85 | 0.45 |
| All time points, window 8 | 0.91 | 0.78 | 0.89 | 0.53 | 0.91 | 0.80 | 0.91 | 0.49 |

**Salehi** (95% intervals in parentheses)

| Metric | Design | Paired | KARST | pairwise | state space | alfak |
|---|---|---|---|---|---|---|
| Forecast beats parent on daughters (Wasserstein) | Two time points | 394 | 0.51 (0.46–0.56) | | | 0.20 (0.17–0.24) |
| | All time points (3–10 passages) | 285 | 0.53 (0.48–0.59) | 0.52 | 0.49 | 0.27 (0.22–0.32) |
| Forecast beats parent on daughters (Euclidean) | Two time points | 394 | 0.44 (0.39–0.49) | | | 0.24 (0.20–0.28) |
| | All time points (3–10 passages) | 285 | 0.42 (0.37–0.48) | 0.44 | 0.43 | 0.26 (0.21–0.31) |
| Fitness of new karyotypes predicts their emergence | Two time points | 480 | 0.56 (0.52–0.61) | | | 0.21 (0.17–0.25) |
| | All time points | 351 | 0.51 (0.46–0.57) | 0.66 | 0.55 | 0.24 (0.20–0.28) |

**Failures and resources** (`analysis/summary/summary_performance.tsv`; alfak's resources are taken from the job records of round 8327282)

| Data | Design | Method | Inputs | Failure fraction | Median time per input (s) | Max (s) | Median memory (MB) | Max (MB) |
|---|---|---|---|---|---|---|---|---|
| ABM | Two time points | KARST | 1200 | 0 | 108 | 1099 | 1006 | 2442 |
| | | alfak | 1200 | 2.2% | 1228 | 5941 | 5436 | 23469 |
| | All time points | KARST | 800 | 0 | 292 | 2200 | 1982 | 7387 |
| | | pairwise | 800 | 0 | 405 | 3969 | 2061 | 8264 |
| | | state space | 800 | 0 | 771 | 17981 | 2474 | 8637 |
| | | alfak | 800 | 1.4% | 3525 | 25142 | 12577 | 88804 |
| Missegregation rate | Two time points | KARST | 375 | 0 | 100 | 627 | 1243 | 2697 |
| | | alfak | 375 | 0 | 949 | 3088 | 5245 | 16095 |
| | All time points | KARST | 250 | 0 | 299 | 2783 | 2891 | 8143 |
| | | pairwise | 250 | 0 | 583 | 2802 | 2944 | 9133 |
| | | state space | 250 | 0 | 1055 | 24245 | 3764 | 12110 |
| | | alfak | 250 | 0 | 4426 | 30045 | 25734 | 147156 |
| Salehi | Two time points | KARST | 253 | 6.3% | 88 | 483 | 957 | 1964 |
| | | alfak | 253 | 6.7% | 529 | 2688 | 2359 | 10009 |
| | All time points | KARST | 187 | 2.7% | 227 | 1550 | 1766 | 5808 |
| | | pairwise | 187 | 2.7% | 263 | 1870 | 1955 | 7213 |
| | | state space | 187 | 2.7% | 705 | 32303 | 2515 | 15073 |
| | | alfak | 187 | 2.7% | 1586 | 12465 | 8003 | 41358 |

- Time is the fit time needed to produce all minobs outputs, excluding forecasts; KARST fits each input only once.
- Memory is the peak of the task process.
- 24 Salehi trajectories were judged unanalysable before fitting and are excluded for all methods.

## 3. Fig. 2d: raw and after fixing sampling noise

**The paper's rule**: fits selected by CV score, ties deduplicated (`analysis/salehi/fig2d_wasserstein.tsv`). Each arm selects a different subset, and the row counts are small, so the arms cannot be compared directly.

| Fit | Rows | Spearman |
|---|---|---|
| KARST two time points | 50 | 0.218 |
| KARST all time points (anchored) | 49 | 0.282 |
| KARST pairwise | 50 | 0.565 |
| KARST state space | 48 | 0.558 |
| alfak dev2@f93b2d2 two time points | 44 | 0.447 |
| alfak dev2@f93b2d2 all time points | 46 | 0.695 |
| Paper (alfak main@4ad82b0) | 45 | 0.688 |

**Rows shared by all arms, without CV selection** (`analysis/fig2d_noise/`):
- "Noise": split the parent and the daughter each at random into two halves and compute the mean distance between the halves (20 splits), computed on the karyotypes of the alfak landscape.
- Intervals: 95% bootstrap intervals over the 6 lineage types.

| Design and subset | Rows | Noise itself | KARST raw | KARST noise fixed | alfak raw | alfak noise fixed |
|---|---|---|---|---|---|---|
| Two time points, all trajectories | 471 | 0.75 | 0.30 | 0.47 | 0.32 | 0.44 |
| Two time points, one interval (2 passages) | 147 | 0.77 | 0.45 | 0.49 | 0.60 | 0.55 |
| All time points, 3–10 passages: anchored | 346 | 0.73 | 0.23 | 0.40 | 0.55 | 0.52 |
| Same as above: pairwise | 346 | | 0.35 | 0.55 | | |
| Same as above: state space | 346 | | 0.45 | 0.56 | | |

Difference from alfak (KARST minus alfak):

| Design and subset | Arm | Raw | Noise fixed |
|---|---|---|---|
| Two time points, all trajectories | KARST | −0.03 (−0.19 to 0.15) | +0.03 (−0.16 to 0.21) |
| Two time points, one interval | KARST | −0.15 (−0.19 to −0.02) | −0.07 (−0.25 to 0.05) |
| All time points | anchored | −0.32 (−0.52 to 0.09) | −0.12 (−0.43 to 0.14) |
| All time points | pairwise | −0.20 (−0.37 to 0.10) | +0.03 (−0.11 to 0.14) |
| All time points | state space | −0.10 (−0.27 to 0.11) | +0.04 (−0.28 to 0.16) |

After fixing the noise, no arm differs reliably from alfak. The differences in the raw Fig. 2d come mainly from how close the forecast's magnitude of change is to the noise.

## 4. Forecasting later passages (`validation/sparse_joint/salehi_horizon`)

- **Forecast**: using this round's fits (each KARST arm, KARST free, alfak), forecast from the last passage of the input. The engine is the same as in this benchmark: the ABM of `alfakR::predict_evo`, 5 replicates, pm 5e-5, nb 1e7. All later daughter passages serve as truth.
- **Score**: log score per cell, relative to "no change" (i.e. the parent's frequencies), mixture weight 0.01, computed on the karyotypes of the fitted landscape.
- **Bins**: binned by days after the parent. Not binned by the parent level in the metadata, because adjacent levels of the p53 lines are 5 to 21 passages apart.
- **Intervals**: 95% bootstrap intervals over lineage types.

**Two time points**: 143 trajectories, 986 rows shared by the three methods.

| Days ahead (rows, trajectories) | KARST vs "no change" | alfak vs "no change" | KARST free | KARST minus alfak | Rows where KARST is higher |
|---|---|---|---|---|---|
| ≤ 15, one PDX passage (341, 99) | +0.61 (0.27 to 0.78) | +0.09 (−0.46 to 0.31) | +0.56 | +0.52 (0.29 to 0.94) | 62% |
| 16–30 (257, 84) | +0.53 (−0.01 to 0.89) | −0.55 (−1.64 to −0.01) | +0.41 | +1.08 (0.57 to 1.94) | 74% |
| 31–60 (257, 68) | +0.06 (−0.88 to 0.61) | −1.32 (−2.39 to −0.55) | +0.04 | +1.39 (1.10 to 1.65) | 75% |
| > 60 (131, 27) | −1.67 (−2.29 to −1.28) | −3.31 (−4.24 to −2.99) | −1.64 | +1.64 (0.71 to 2.96) | 74% |
| All (986, 143) | | | | +1.04 (0.71 to 1.43; by sample 0.75 to 1.29) | 70% |

- KARST free minus KARST: −0.05 (−0.15 to 0.11); KARST free is higher in 38% of rows.
- Rows beyond 60 days come from only 3 lineage types, so the interval in that bin says little.

**All time points**: 99 trajectories with 3 or more passages, 677 rows shared by the four arms.

| Days ahead (rows) | anchored | pairwise | state space | alfak |
|---|---|---|---|---|
| ≤ 15 (259) | +0.95 (0.58 to 1.10) | +0.96 (0.58 to 1.12) | +0.98 (0.58 to 1.11) | +0.38 (0.31 to 0.62) |
| 16–30 (184) | +0.99 (0.24 to 1.47) | +0.97 (0.21 to 1.46) | +0.94 (0.27 to 1.38) | −0.21 (−0.57 to 0.12) |
| 31–60 (168) | +0.66 (−0.09 to 1.37) | +0.50 (−0.26 to 1.18) | +0.24 (−0.22 to 0.74) | −1.40 (−2.09 to −0.68) |
| > 60 (66) | −1.57 (−1.83 to −1.53) | −1.79 (−2.42 to −1.70) | −1.95 (−2.97 to −1.79) | −2.58 (−2.67 to −2.57) |

| Paired difference (all 677 rows) | Difference | Lineage interval | Rows where the first is higher |
|---|---|---|---|
| anchored minus alfak | +1.15 | 0.59 to 1.51 | 71% |
| pairwise minus alfak | +1.09 | 0.49 to 1.48 | 71% |
| state space minus alfak | +1.01 | 0.44 to 1.37 | 71% |
| state space minus anchored | −0.14 | −0.20 to −0.05 | 40% |
| pairwise minus anchored | −0.06 | −0.10 to −0.04 | 45% |
| free minus anchored | −0.14 | −0.23 to −0.01 | 33% |

**All time points compared with two time points, on the common support set** (`report_all/paired_all_vs_ends.csv`)

- **Common support set**: the union of the landscape karyotypes of all fits of both designs for the same trajectory and the same minobs; both designs are scored on this set against the same "no change" baseline.
- **Rows**: 636 rows from 99 trajectories.

| Days ahead (rows) | anchored: all minus two time points | Rows higher | alfak: all minus two time points |
|---|---|---|---|
| ≤ 15 (237) | +0.55 (0.18 to 0.66) | 98% | +0.25 (−0.01 to 0.77) |
| 16–30 (173) | +0.60 (0.10 to 1.03) | 90% | +0.23 (−0.38 to 1.40) |
| 31–60 (160) | +0.49 (0.07 to 0.91) | 80% | −0.48 (−1.98 to 0.70) |
| > 60 (66) | +0.16 (0.15 to 0.23) | 73% | +1.03 (−0.46 to 1.26) |
| All (636) | +0.51 (0.13 to 0.79; by sample 0.17 to 0.84) | 89% | +0.14 (−0.45 to 0.85), 53% |

- **state space compared with anchored on two time points**: +0.37 (−0.03 to 0.68; by sample 0.05 to 0.67), higher in 76% of rows.
- **Gain of each relative to "no change" on this support set**:
  - anchored: all time points +0.64, two time points +0.13;
  - alfak: all time points −0.56, two time points −0.70.
- **Conclusion: fitting with all time points, KARST forecasts later passages better than with only two time points, while alfak shows no reliable improvement.** This time the two designs are scored on the same support set, so the difference is not caused by different landscape sizes.

## 5. Fig. 4–6 compared with the paper's conclusions (`analysis/fig456/conclusions.tsv`)

1.1.0 changed the comparison rule; there are now 16 conclusions per arm, which cannot be compared directly with the 19 of e7e32f8. The reference is the paper's own fits (alfak main@4ad82b0).

| Arm | Agree |
|---|---|
| KARST two time points | 10 / 16 |
| KARST all time points (anchored) | 7 / 16 |
| KARST pairwise | 8 / 16 |
| KARST state space | 8 / 16 |
| alfak dev2@f93b2d2 two time points | 12 / 16 |
| alfak dev2@f93b2d2 all time points | 13 / 16 |

Compared with 2c60af1 (same comparison rule), the significance of 3 conclusions flipped, and all 3 intervals lie right against 0:
- 4d `cisplatin_v`, pairwise: 0.264, interval −0.010 to 0.479, flipped to not significant;
- 4d `cisplatin_v`, state space: 0.281, interval 0.013 to 0.613, flipped to significant;
- 4f `log1p(dk):pair_typeparallel`, anchored: −0.080, interval −0.180 to −0.010, flipped to significant.

This report looks only at accuracy and does not interpret the biological meaning of these conclusions.

## 6. Posterior mode search and history choice

`history_profile` is read fit by fit (script `build/restart_counts_cd342db.R`, result `build/restart_counts_cd342db.tsv`, both only on HPC). The table counts only successful fits; "candidates" means candidate histories.

| Data, design and method | Fits | Candidates | Needing the second search (fits / candidates) | Failed both (candidates) | New rule chooses differently from nlminb |
|---|---|---|---|---|---|
| ABM two time points | 1200 | 6759 | 19 / 19 | 0 | 14 |
| ABM all time points anchored | 800 | 4359 | 13 / 14 | 0 | 18 |
| ABM pairwise | 800 | 4359 | 11 / 11 | 0 | 3 |
| ABM state space | 800 | 4359 | 8 / 8 | 0 | 10 |
| Missegregation rate two time points | 375 | 2123 | 6 / 6 | 0 | 4 |
| Missegregation rate anchored | 250 | 1373 | 3 / 3 | 0 | 1 |
| Missegregation rate pairwise | 250 | 1373 | 4 / 4 | 0 | 0 |
| Missegregation rate state space | 250 | 1373 | 2 / 2 | 0 | 0 |
| Salehi two time points | 253 | 1284 | 1 / 1 | 0 | 0 |
| Salehi anchored | 187 | 899 | 2 / 2 | 0 | 1 |
| Salehi pairwise | 191 | 919 | 0 / 0 | 6 | 0 |
| Salehi state space | 191 | 919 | 2 / 2 | 0 | 1 |
| Total | 5547 | 30099 | 71 / 72 | 6 | 52 |

- **Second search**: 72 candidates needed it, 0.24%, spread over 1.3% of the fits. Only 1 fit (Salehi, two time points) selected such a candidate.
- **Failed both**: 6 candidates, all in Salehi pairwise fits; none was selected, and their fits were completed with other candidates.
- **New rule chooses differently from nlminb**: 52 fits (0.9%), none of which selected a candidate that needed the second search.
- **KARST free fits** do not compare candidate histories, and no second search was used in fitting.
- **The 84 fits that failed in 2c60af1 and were recovered in this round**:
  - 77 are fits in the table above that had a second search, or a candidate that failed both. In 1.1.0, if the search of any one candidate failed, the whole fit raised an error.
  - The other 7 fits were themselves successful in 2c60af1; the failure occurred in the search when `as_alfak()` subsequently rebuilt the posterior, so no minobs had any output. 1.1.1 likewise searches once more during the rebuild.
- **Unrecorded retries**: the retries when `as_alfak()` rebuilds the posterior and the retries in the cross-validation folds are not written to `run.Rds`; the table above does not include them.

**Effect of the history changes on accuracy**:

- **How compared**: against round 2c60af1. It has the same correct_efflux formula as this round; for fits whose history did not change, the landscapes of the two rounds are exactly the same.
- **Subjects**: the 44 inputs from ABM and the missegregation-rate simulations whose history changed and which were fit successfully in both rounds, each input counted once.

| Method | Inputs | Median change in Pearson r | Better / worse | Spearman ρ better / worse |
|---|---|---|---|---|
| ABM two time points | 11 | 0.000 | 7 / 4 | 9 / 2 |
| ABM anchored | 16 | +0.001 | 10 / 6 | 11 / 5 |
| ABM pairwise | 3 | +0.036 | 2 / 1 | 2 / 1 |
| ABM state space | 9 | −0.014 | 2 / 7 | 2 / 7 |
| Missegregation rate two time points | 4 | 0.000 | 1 / 3 | 2 / 2 |
| Missegregation rate anchored | 1 | +0.109 | 1 / 0 | 1 / 0 |
| Total | 44 | 0.000 | 23 / 21 (p = 0.57) | 27 / 17 (p = 0.29) |

- The 10% and 90% quantiles of the change are −0.017 and +0.042; the largest single change is +0.109.
- state space got worse in 7 / 9 inputs, but there are only 9 inputs, so on its own this says little.
- Overall consistent with the `validation/sparse_joint/history_choice` study: changes in history choice are rare, and accuracy shows no consistent direction.

## 7. Fit-by-fit comparison with previous rounds (`analysis/compare_rounds/`)

| Method | Compared with | Fits (by minobs) | Status | Both succeeded, history unchanged | Both succeeded, history changed |
|---|---|---|---|---|---|
| anchored | 2c60af1 | 8071 | failure → success 129, success → failure 0, both failed 189 | 7664: landscapes and ABM forecasts exactly the same; 6 differ in cross-validation | 89: landscapes differ |
| pairwise | 2c60af1 | 3265 | failure → success 54, success → failure 0, both failed 59 | 3143: exactly the same | 9 |
| state space | 2c60af1 | 3265 | failure → success 38, success → failure 0, both failed 59 | 3138: landscapes and forecasts the same; 15 differ in cross-validation | 30 |
| anchored | e7e32f8 | 8071 | all the same (both failed 189) | 7778 | 104 |
| pairwise | e7e32f8 | 3265 | all the same (both failed 59) | 3197 | 9 |
| state space | e7e32f8 | 3265 | all the same (both failed 59) | 3173 | 33 |

- **Compared with 2c60af1**:
  - every fit whose landscape differs had its history changed;
  - the 21 fits whose history did not change and which differ only in cross-validation should come from the 1.1.1 retry: when the mode search of a cross-validation fold failed the first time, 1.1.0 dropped that fold, while 1.1.1 searches once more. Between the two rounds' code, this is the only change that alters cross-validation without altering the landscape.
- **Compared with e7e32f8**: all landscapes differ, which is caused by the correct_efflux of 1.1.0; fit success and failure are exactly the same.
- **"Both failed"**: besides the 112 failed fits of Section 1, this also includes cases where the fit succeeded but `as_alfak()` produced no output at some minobs. The ones that failed in both rounds are the same set.

## 8. Runs

Job IDs are recorded in `build/cd342db_jobs.tsv` on HPC:
- build 22111287;
- fits 22111330–22111350 (11 tiers) and KARST free 22111352; reruns 22124841, 22124843;
- scoring 22111354, collect_analyse 22111356, Fig. 4 22111357 / 59, Fig. 6 22111461 / 63, fig456_summary 22111465, Fig. 2d noise 22111466 / 68;
- comparisons 22111470–22111475;
- multi-horizon test 22111476 / 78, analysis 22111480 / 82;
- the counts for this report 22141113, not part of this round.

**Timeline**:
- Fits ran from 2026-11-19 14:22 to 2026-11-20 00:02. The last job finished at 01:03.
- The longest task was Salehi SA609_X10_l_10 (all time points, state space, 32G tier), about 9 hours.
- Node failures (NODE_FAIL) hit 20 tasks in total, all of which completed after Slurm requeued them automatically: 12 in the ABM 4G tier, 6 in the 8G tier, 2 in the multi-horizon test.

**The 3 rerun tasks**: they failed at container start ("unknown userid", nodes hpctpa3pc0028, 0021, 0024):
- ABM 16G tier task 708;
- ABM 32G tier tasks 137, 488.

They were rerun at 17:16 with `build/resubmit_cd342db.sh`, and collect_analyse and the comparison jobs were made to wait for the reruns to finish before starting.

All jobs ran on qos xxlarge, 1 CPU per task, and every array had an add-hours monitor.

## 9. Files

- `analysis/summary/`, `analysis/abm/`, `analysis/salehi/`, `analysis/fig456/`: paper tests.
- `analysis/fig2d_noise/`: the raw and noise-fixed versions of Fig. 2d.
- `analysis/compare_rounds/`: fit-by-fit comparison with 2c60af1 and e7e32f8.
- `validation/sparse_joint/salehi_horizon/results/report/` and `report_all/` (inside the karst clone on HPC): tables and figures of the multi-horizon test. The results of rounds e7e32f8 and 2c60af1 are in `results/e7e32f8/`, `results/2c60af1/`.
- `build/restart_counts_cd342db.tsv` (HPC): the fit-by-fit records of Section 6.
