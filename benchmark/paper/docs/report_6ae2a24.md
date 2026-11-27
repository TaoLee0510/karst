# Comparison report of KARST 6ae2a24 and alfak on the ALFA-K paper tests

Date: 2026-11-20. KARST 6ae2a24 (1.1.3, dev). alfak dev2@f93b2d2 is the reference, with its fits linked from the alfak2 era; the paper's own fits are denoted alfak main@4ad82b0. The results directory is on HPC: `/share/lab_crd/taoli/Project/karst/benchmark/paper/results_6ae2a24/`.

> Terminology: the input columns are called "time points". A Salehi "trajectory" consists of several passages. "Fits (by minobs)" means one row per minobs; the minobs values of the same input share one fit. "Draws" means the `bootstrap_res` and `landscape_posterior_samples` that `as_alfak()` draws from the posterior.

## 0. What this round did

The previous report was cd342db (1.1.2). While re-examining the 7 fits that failed in `as_alfak()` in round 2c60af1, a defect was found in `sparse_joint_rebuild()` of 1.1.1 and 1.1.2.

**Defect**: when `as_alfak()` rebuilds the posterior for the draws, it runs a new search starting from the stored mode. All 5987 fits of cd342db were checked (on HPC only, `build/check_rebuild_cd342db.R`, `build/rebuild_scan_cd342db.R`):
- **7 fits fail the search at the stored mode**:
  - Cause: the gradient there is about 1e-9, but the smallest eigenvalue of the inner Hessian is about 1e-8, so TMB judges it "not positive definite" and returns NaN.
  - Consequence: the 1.1.1 retry started from the initial values and landed on a worse mode. The two most severe:
    - Salehi SA039U_X60_l_4, all time points: objective +14.1, g off by 0.13;
    - ABM w_0p8_rep_72, w8_all: +10.3, g off by 0.033.
- **In another 3 fits, the stored mode itself is not a stationary point**: the fit's own search stopped early, for example ABM w_0p8_rep_52, w8_all, state space, gradient 6.25. The search during the rebuild carried on from there and moved g by up to 0.40.
- **Scope of the effect**: only the draws are affected, including their correlations and the dispersion of the initial abundances. The landscape mean, median and sd are taken directly from the fit and are not affected.

**Fix (1.1.3, 6ae2a24)**:
- `sparse_joint_rebuild()` takes the fit's own posterior directly: it evaluates at the stored mode (`objective$env$par`), without any replacement and without searching again.
- Each fit records the maximum absolute value of the latent-variable gradient at that point, written in `rebuild$gradient` of `run.Rds`.
- The cross-validation folds still search for their own modes (searching once more if the first attempt fails). They are counted in three classes, written in `rebuild$xval_folds`:
  - folds where the mode was found;
  - of these, those where only the second search found it;
  - folds where neither search found it; these folds are not counted in the cross-validation score.

**Work in this round**:
- All KARST fits were redone: ABM, missegregation-rate simulations and Salehi, in two designs (two time points; all time points, including pairwise and state space), 5617 tasks in total. In addition, 482 Salehi fits of KARST `initial_state = "free"`.
- alfak was not refitted: 3107 fit folders were linked from round 8327282 of the alfak2 era.
- The multi-horizon prediction test was rerun for all methods.
- The fits were compared one by one with rounds cd342db and e7e32f8.
- The development benchmark is not part of this round.

## 1. Conclusions

1. **This round is complete.**
   - All 6099 fit tasks ran to completion.
   - 6 tasks hit node failures (NODE_FAIL) and completed after Slurm requeued them automatically.
   - 2 tasks failed at container start-up and were rerun (section 6).
   - There were no timeouts and no out-of-memory failures.
2. **1.1.3 did not change any fit (section 2).**
   - Compared with cd342db, for the fits that succeeded on both sides under the three methods, the following are all identical:
     - success/failure status and the chosen history;
     - landscape (maximum difference of the mean 4.4e-16, of the sd 3.6e-15);
     - cross-validation;
     - forecasts (ABM and Salehi).
   - The per-fit records of candidate histories (number of candidates, second searches, choice result) are also identical.
3. **The results that differ from cd342db are not caused by any change of fit (section 2).** The causes fall into three classes:
   - **Draws**: the draws are new random draws in every fit, so every result that uses them changes with them. This concerns the example projections of Supp. Fig. 2d–e and the dispersion of the initial abundances.
   - **Ties in rank correlations**: landscape differences of order 1e-16 broke the ordering of tied values and changed Spearman ρ by up to 0.027; Pearson r is unchanged.
   - **Logistic regressions with complete separation**: in the Salehi new-karyotype test, 38 regressions have coefficients between several thousand and 1e18, with no definite solution. They give "the fitness of a new karyotype predicts its appearance" one more significant regression in each arm, for example from 0.565 to 0.567 for two time points.
   - **Completely unchanged**: the multi-horizon test, Fig. 2d (original and fixed-noise) and the comparison with the paper's conclusions have output that is byte-identical to cd342db. The ABM and Salehi forecast scores are also identical; only the cross-validation scores in the tables have rounding differences of no more than 3e-10.
4. **The stored modes are almost all stationary points (section 3).**
   - Of the 5985 fits with rebuild records, 4 have a gradient above 1e-3 and 124 above 1e-6.
   - These 4 were all already flagged convergence = 1 by the convergence check of step 7. The largest is ABM w_0p8_rep_52 (w8_all, state space), gradient 6.25.
   - Cross-validation: 174731 folds found the mode, 23 of them only with the second search; 90 folds did not find it, 84 of them in ABM pairwise fits (78 in window 8).
5. **The main conclusions are the same as for cd342db (sections 4, 5).**
   - ABM and missegregation-rate simulations: KARST's landscape is more accurate than alfak's in every window.
   - Salehi: Fig. 2d mainly reflects sampling noise; after the noise is fixed, no KARST arm differs reliably from alfak.
   - At the passages after the forecast, KARST beats alfak in both designs. Fitting with all time points forecasts better than using only two time points.

## 2. Fit-by-fit comparison with cd342db (`analysis/compare_rounds/`)

| Method | Compared with | Fits (by minobs) | Status | Succeeded on both sides | History | Landscape: maximum difference of mean / sd | Cross-validation | Forecasts (ABM and Salehi) |
|---|---|---|---|---|---|---|---|---|
| anchored | cd342db | 8071 | all identical | 7882 | all identical | 2.2e-16 / 8.9e-16 | all identical | all 6681 identical |
| pairwise | cd342db | 3265 | all identical | 3206 | all identical | 4.4e-16 / 3.6e-15 | all identical | all 2687 identical |
| state space | cd342db | 3265 | all identical | 3206 | all identical | 2.2e-16 / 8.9e-16 | all identical | all 2687 identical |

- The landscape differences of order 1e-16 arise because the fits ran on different nodes, which round slightly differently.
- The comparison with e7e32f8 gives the same result as round cd342db's comparison with e7e32f8: all landscapes differ (caused by the correct_efflux of 1.1.0), and the success or failure of the fits is identical.

**Outputs that changed, and the reasons**:

| Output | Change | Reason |
|---|---|---|
| Spearman ρ in the ABM scores (`rho`, `rhofq`, `rhonn`, `rhod2`) | up to 0.027 (`rhonn`); the medians in the summary tables are unchanged, only the last digit of a few p values changed | ties: landscape differences of 1e-16 broke the ordering of ties. The subset with few karyotypes (`rhonn`) changes the most |
| Pearson r and R² in the ABM scores, cross-validation scores in the ABM and Salehi tables | ≤ 3e-10 | rounding |
| ABM paired win fractions (`paired_landscape.tsv`) | a few windows differ by one input | when the scores of the two methods are equal, the win flips with rounding |
| Mean absolute error of the example projections in Supp. Fig. 2d–e | KARST two time points 0.0094 → 0.0098; one pairwise row also changed | the projections start from each draw |
| Dispersion of the initial abundances (`note32_initial_frequencies.tsv`) | last digit | draws |
| Salehi new karyotypes (Fig. 3e–f) | one more significant regression in each arm: two time points 0.565 → 0.567, all time points 0.513 → 0.516, pairwise 0.661 → 0.664, state space 0.547 → 0.550 | 38 logistic regressions with complete separation, with coefficients between 6e3 and 1e18 that vary arbitrarily with input differences of 1e-16; another 4 coefficients are NA because of collinearity |
| Fig. 4 estimates, Fig. 5 models | ≤ 7e-9; no significance flipped | rounding |
| Run time and memory | see section 4 | different node load |

**The draws are indeed new**: the `bootstrap_res` of the two rounds were compared in 8 ABM fits, including the Supp. Fig. 2 example rep_74 (script `build/draws_diff_6ae2a24.R`, HPC). The maximum difference of the fitness draws is between 0.11 and 0.43, which is 5 to 13 times the draw sd (median 0.016 to 0.034). `compare_rounds.R` does not compare draws in the first place.

## 3. Records of the posterior rebuild

The `rebuild` of `run.Rds` is read fit by fit (script `build/restart_counts_6ae2a24.R`, results `build/restart_counts_6ae2a24.tsv`, both only on HPC; local copy in `results_6ae2a24/report_extra/`).

**Fits with records**: 5985 of the 5987 successful fits. The other 2 are Salehi SA609UnBU_X7_l_3_d1_0_d2_0 (two time points, KARST and KARST free): every minobs failed with "no frequent karyotypes", so there is no rebuild.

| Data, design and method | Fits | Gradient > 1e-6 | Gradient > 1e-3 | Maximum gradient | Cross-validation: folds with the mode found | Of which by the second search | Not found |
|---|---|---|---|---|---|---|---|
| ABM two time points (window 2 / 4 / 8) | 1200 | 16 | 0 | 5.4e-5 | 33460 | 3 | 0 |
| ABM all time points anchored | 800 | 12 | 0 | 6.9e-4 | 25760 | 0 | 0 |
| ABM pairwise | 800 | 2 | 0 | 5.6e-6 | 28871 | 0 | 84 |
| ABM state space | 800 | 19 | 2 | 6.25 | 28952 | 16 | 3 |
| Missegregation rate two time points | 375 | 13 | 0 | 2.8e-5 | 6298 | 0 | 0 |
| Missegregation rate anchored | 250 | 11 | 0 | 2.7e-4 | 4955 | 0 | 0 |
| Missegregation rate pairwise | 250 | 5 | 0 | 7.2e-6 | 5959 | 0 | 3 |
| Missegregation rate state space | 250 | 14 | 0 | 4.0e-5 | 5962 | 1 | 0 |
| Salehi two time points | 252 | 5 | 0 | 2.9e-5 | 5043 | 0 | 0 |
| Salehi two time points free | 252 | 2 | 0 | 2.1e-5 | 5043 | 0 | 0 |
| Salehi anchored | 187 | 6 | 0 | 3.5e-4 | 5422 | 3 | 0 |
| Salehi anchored free | 187 | 3 | 0 | 2.2e-6 | 5422 | 0 | 0 |
| Salehi pairwise | 191 | 5 | 2 | 0.010 | 6792 | 0 | 0 |
| Salehi state space | 191 | 11 | 0 | 6.4e-5 | 6792 | 0 | 0 |
| Total | 5985 | 124 | 4 | 6.25 | 174731 | 23 | 90 |

- **The median gradient** is 5e-11 to 3e-9 across the rows: the stored modes are almost all stationary points.
- **The 4 fits with a gradient above 1e-3** were all already flagged convergence = 1 by the convergence check of step 7:
  - ABM w_0p8_m_5e-05_rep_52, w8_all, state space: 6.25;
  - ABM w_0p2_m_5e-05_rep_19, w8_all, state space: 0.026;
  - Salehi SA535_CISPLATIN_CombinedH_X9_l_5, all time points, pairwise: 0.010;
  - Salehi SA609_X6_l_3, all time points, pairwise: 0.0011.
  The draws of these fits are centred on the point that the fit reports, consistent with the reported landscape. That point itself is not an exact mode, which the convergence flag already states.
- **The 90 folds where no mode was found**: 84 in ABM pairwise (78 in window 8, 6 in window 4), 3 in ABM state space, 3 in missegregation-rate pairwise. These folds are not counted in the cross-validation score. The folds were handled in the same way in 1.1.2, so the cross-validation agrees completely with cd342db.
- **Candidate histories**: the columns outside this table (number of candidates, candidates that needed a second search, candidates that failed both searches, the chosen history, whether it differs from the nlminb choice) are identical, one by one, to the records of cd342db. See section 6 of `docs/report_cd342db.md`.

## 4. Paper tests: summary accuracy (`analysis/summary/summary_accuracy.tsv`)

All tables take medians over the paired inputs.

- "KARST better" column: among the paired inputs, the fraction where KARST (two time points, or anchored with all time points) is better than alfak.
- "Better than anchored" column: among the paired inputs, the fraction where pairwise or state space is better than anchored.

**ABM, landscape Pearson r**

| Design and window | Pairs | KARST | pairwise | state space | alfak | KARST better | Better than anchored (pairwise / state space) |
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

**ABM, fraction of forecasts that beat "no change" at the 10th passage**

| Design and window | Euclidean: KARST | pairwise | state space | alfak | Wasserstein: KARST | pairwise | state space | alfak |
|---|---|---|---|---|---|---|---|---|
| Two time points, window 2 | 0.70 | | | 0.36 | 0.69 | | | 0.32 |
| Two time points, window 4 | 0.85 | | | 0.45 | 0.87 | | | 0.44 |
| Two time points, window 8 | 0.86 | | | 0.40 | 0.88 | | | 0.35 |
| All time points, window 4 | 0.86 | 0.77 | 0.84 | 0.50 | 0.87 | 0.71 | 0.85 | 0.45 |
| All time points, window 8 | 0.91 | 0.78 | 0.89 | 0.53 | 0.91 | 0.80 | 0.91 | 0.49 |

**Salehi** (95% interval in parentheses)

| Metric | Design | Pairs | KARST | pairwise | state space | alfak |
|---|---|---|---|---|---|---|
| Forecast beats the parent on the daughter (Wasserstein) | Two time points | 394 | 0.51 (0.46–0.56) | | | 0.20 (0.17–0.24) |
| | All time points (3–10 passages) | 285 | 0.53 (0.48–0.59) | 0.52 | 0.49 | 0.27 (0.22–0.32) |
| Forecast beats the parent on the daughter (Euclidean) | Two time points | 394 | 0.44 (0.39–0.49) | | | 0.24 (0.20–0.28) |
| | All time points (3–10 passages) | 285 | 0.42 (0.37–0.48) | 0.44 | 0.43 | 0.26 (0.21–0.31) |
| Fitness of new karyotypes predicts their appearance | Two time points | 480 | 0.57 (0.52–0.61) | | | 0.21 (0.17–0.25) |
| | All time points | 351 | 0.52 (0.46–0.57) | 0.66 | 0.55 | 0.24 (0.20–0.28) |

**Failures and resources** (`analysis/summary/summary_performance.tsv`; alfak's resources are taken from the job records of round 8327282)

| Data | Design | Method | Inputs | Failure fraction | Median time per input (s) | Maximum (s) | Median memory (MB) | Maximum (MB) |
|---|---|---|---|---|---|---|---|---|
| ABM | Two time points | KARST | 1200 | 0 | 117 | 968 | 976 | 2238 |
| | | alfak | 1200 | 2.2% | 1228 | 5941 | 5436 | 23469 |
| | All time points | KARST | 800 | 0 | 354 | 4200 | 1866 | 6431 |
| | | pairwise | 800 | 0 | 415 | 3403 | 1907 | 7372 |
| | | state space | 800 | 0 | 729 | 17818 | 2363 | 8030 |
| | | alfak | 800 | 1.4% | 3525 | 25142 | 12577 | 88804 |
| Missegregation rate | Two time points | KARST | 375 | 0 | 100 | 682 | 1179 | 2431 |
| | | alfak | 375 | 0 | 949 | 3088 | 5245 | 16095 |
| | All time points | KARST | 250 | 0 | 311 | 2758 | 2743 | 7638 |
| | | pairwise | 250 | 0 | 561 | 2648 | 2830 | 7790 |
| | | state space | 250 | 0 | 989 | 21302 | 3625 | 11236 |
| | | alfak | 250 | 0 | 4426 | 30045 | 25734 | 147156 |
| Salehi | Two time points | KARST | 253 | 6.3% | 92 | 543 | 872 | 2161 |
| | | alfak | 253 | 6.7% | 529 | 2688 | 2359 | 10009 |
| | All time points | KARST | 187 | 2.7% | 237 | 1771 | 1671 | 5012 |
| | | pairwise | 187 | 2.7% | 268 | 1954 | 1773 | 6407 |
| | | state space | 187 | 2.7% | 688 | 31821 | 2421 | 13415 |
| | | alfak | 187 | 2.7% | 1586 | 12465 | 8003 | 41358 |

- Time is the fit time needed to produce the output for all minobs, excluding forecasts; KARST fits each input only once.
- Memory is the peak of the task process.
- Compared with cd342db:
  - the median times differ by −8% to +21%, with no consistent direction: anchored became slower and state space faster. The node load differed, so from this round alone the change cannot be attributed to the code.
  - the median memory of the KARST rows is 3% to 9% lower, consistent with 1.1.3 no longer searching again for the draws, but this was not measured separately.
- 24 Salehi trajectories were judged unanalysable before fitting and are excluded for all methods.

## 5. Fig. 2d, the later passages and Fig. 4–6

The inputs of these three parts (Salehi landscapes, forecasts) are the same as in cd342db, and the analysis outputs are also byte-identical to cd342db. So only the key points are given below; the full tables are in sections 3–5 of `docs/report_cd342db.md`, and those numbers are also the numbers of this round.

**Fig. 2d** (`analysis/salehi/fig2d_wasserstein.tsv`, `analysis/fig2d_noise/`)

- Under the paper's rule (picking by CV score), the Spearman values are:
  - KARST: two time points 0.218, anchored 0.282, pairwise 0.565, state space 0.558;
  - alfak: two time points 0.447, all time points 0.695.
- After fixing the sampling noise on the rows that all arms share, no arm differs reliably from alfak. For example, two time points, all trajectories: KARST minus alfak is +0.03 (−0.16 to 0.21).

**The later passages after the forecast** (`validation/sparse_joint/salehi_horizon/results/6ae2a24/`)

| Comparison | Rows | Difference in per-cell log score | Lineage interval | Rows where the former is higher |
|---|---|---|---|---|
| Two time points: KARST minus alfak | 986 | +1.04 | 0.71 to 1.43 | 70% |
| Two time points: KARST free minus KARST | 986 | −0.05 | −0.15 to 0.11 | 38% |
| All time points: anchored minus alfak | 677 | +1.15 | 0.59 to 1.51 | 71% |
| All time points: pairwise minus alfak | 677 | +1.09 | 0.49 to 1.47 | 71% |
| All time points: state space minus alfak | 677 | +1.01 | 0.44 to 1.37 | 71% |
| Common support set: anchored, all minus two time points | 636 | +0.51 | 0.13 to 0.79 | 89% |
| Common support set: alfak, all minus two time points | 636 | +0.14 | −0.45 to 0.85 | 53% |

**Comparison of Fig. 4–6 with the paper's conclusions** (`analysis/fig456/conclusions.tsv`, 16 per arm, with the paper's own fits as reference)

| Arm | Agree |
|---|---|
| KARST two time points | 10 / 16 |
| KARST all time points (anchored) | 7 / 16 |
| KARST pairwise | 8 / 16 |
| KARST state space | 8 / 16 |
| alfak dev2@f93b2d2 two time points | 12 / 16 |
| alfak dev2@f93b2d2 all time points | 13 / 16 |

This report looks only at accuracy and does not interpret the biological meaning of these conclusions.

## 6. Runs

The job numbers are recorded in `build/6ae2a24_jobs.tsv` on HPC:
- build 22143738;
- fits 22143765–22143785 (11 tiers) and KARST free 22143787; rerun 22154608;
- scoring 22143789, collect_analyse 22143791, Fig. 4 22143792 / 94, Fig. 6 22144288 / 90, fig456_summary 22144292, Fig. 2d noise 22144293 / 95;
- comparisons 22144297–22144302;
- multi-horizon test 22144303 / 05, analysis 22144307 / 09;
- counting of the rebuild records 22154613.

**Timeline**:
- 2026-11-20 02:16 build, 02:19 submission.
- The fits ended at 11:49. The last job, fig456_summary, ended at 12:54.
- The longest task was Salehi SA609_X10_l_10 (all time points, state space, 32G tier), 8 hours 50 minutes.

**Failures**:
- **Node failures (NODE_FAIL)**: 6 tasks in total, all completed after Slurm requeued them automatically:
  - ABM 4G tier: 1;
  - ABM 8G tier: 4;
  - ABM 16G tier: 1.
- **Container start-up failures**: 2 tasks ("unknown userid", nodes hpctpa3pc0023, 0016):
  - ABM 32G tier task 412 (w_0p8_m_5e-05_rep_06, w8_all, state space);
  - ABM 32G tier task 541 (w_0p8_m_5e-05_rep_70, w4_all, state space).

  Rerun at 07:25 with `build/resubmit_6ae2a24.sh`, excluding these two nodes. The collect_analyse and comparison jobs were changed to start only after the rerun had completed.

All jobs ran on qos xxlarge, with 1 CPU per task, and every array had an add-hours monitor.

## 7. Files

- `analysis/summary/`, `analysis/abm/`, `analysis/salehi/`, `analysis/fig456/`: paper tests.
- `analysis/fig2d_noise/`: the original and fixed-noise versions of Fig. 2d.
- `analysis/compare_rounds/`: fit-by-fit comparisons with cd342db and e7e32f8.
- `validation/sparse_joint/salehi_horizon/results/6ae2a24/` (inside the karst clone on HPC): tables and figures of the multi-horizon test. `report/` and `report_all/` were also synced to the same local path.
- The following files are only in `build/` on HPC:
  - `restart_counts_6ae2a24.tsv`: the per-fit records of section 3;
  - `score_diff_6ae2a24.R`: the source of the Spearman differences in section 2;
  - `draws_diff_6ae2a24.R`: the comparison of the draws in section 2;
  - `nonstat_6ae2a24.R`: the convergence flags of the 4 fits in section 3.
