# dev autorun: preregistered criteria and decision rules

Written: 2026-11-19 02:36 EST, before any diagnosis or scoring. User approval: 2026-11-19 ~02:30 EST (phases 0–3 all run automatically; no commits during the run; no HPC limit, use xxlarge; location is `dev/` at the repository root).

Source: `doc/literature/karst_limitations_review.html` (candidates R1–R10, diagnostics D1–D6).

**Rules**
- The criteria in this file are fixed before the corresponding results are seen, and are not modified afterwards. If a criterion itself turns out to be unreasonable, say so plainly in the report, do not change the criterion, treat it as "not passed", and stop and report.
- "Not evaluable" (missing data, missing baseline, cannot be computed) is treated as "not passed".
- The package itself (`R/`, `src/`, `tests/`, defaults) is not changed; model changes are kept as patches in `dev/patches/`, applied on top of baseline b319cfd, and built into a separate library.
- The strongest possible conclusion is "recommend adoption (as an optional argument, default unchanged)"; nothing is ported into the package.
- Clarification (fixed before scoring): in the autorun plan, D4 originally said "compute the dispersion after fitting at the true rates". Because each non-major karyotype has two latent variables for its two counts, the fitted values are almost saturated and the residual dispersion would be underestimated; D4 therefore uses instead the true expectation obtained by rerunning the ABM over the same window (see below). The threshold values are unchanged.

## General definitions

- **Baseline**: karst b319cfd (package code identical to 2c60af1), without patches, labelled `base`.
- **e sets**: e1 passaged culture (288 datasets), e2 clonal replacement, e3 neutral culture, e4 true missegregation rate 4×/16×. Fits use the first and last time points (`ends`). Inputs are read only from the HPC backup directory `alfak2_results_backup_2026-11-17`.
- **Grid cells**: as defined in `validation/sparse_joint/beta_choice/analyse_series.R` (24 grid cells for e1: bottleneck × interval × number of time points × cell count).
- **Distance**: L1 distance to the major karyotypes, 0, 1, 2 (the far-distance study also has 3, 4).
- **Scoring**: `stt_score_predictions()`; the truth is the net growth rate at the true rates (`sm_true_net_rate()`).
- **Coverage criterion (L4 general, hereafter "the general coverage criterion")**: nominal 95% intervals. Categories: observed with cells <5; observed with cells ≥5; unobserved distance 1; unobserved distance 2. Each non-stress grid cell (category × scenario × wavelength or family × design) has coverage ≥0.90, the pooled coverage of each category is within 0.93–0.97, the median width is ≤ 1.25 times the current one, and at least 90% of datasets have their own coverage ≥0.85. Stress grid cells are the scenarios after 1080 days and e4.

## Phase 0: diagnostics

**D6+D4 process-noise study** (a copy of the e-set ABM with passage records added, `dev/abm/`)
- Neutral part: constant fitness 0.9, β = 0, 2e6 starting cells split evenly among 10 different karyotypes; bottleneck 200 / 2000 / 20000; runs of 60 days and 360 days; 40 replicates per combination. For each karyotype F = Var(p_T)/[p₀(1−p₀)] (across replicates), n is the number of passages, and the per-passage constant c_B = B·[1 − (1−F)^(1/n)]. c = the median of c_B over the three bottlenecks; once fixed, it is used by R1.
  - **Rule G0**: c not in [0.5, 4] → stop and report.
- Passage records: for each e1–e4 dataset, replay the window with the same seed and record the number of survivors at each passage between the first and last time points; check that the time-point populations from the replay agree with the dataset's true frequencies (disagreement → stop and report).
- Landscape part: 12 landscapes of e1 × 3 bottlenecks × 2 intervals, 30 replicates per window, starting from the same t₀ population, taking the 4th and 8th time points. Per-karyotype process F_k = Var_rep(p)/[m(1−m)], with m the replicate mean; with n₁ cells the count dispersion is D_k = 1 + (n₁−1)F_k.
  - **Rule G1**: at bottleneck 200 and n₁ = 1000, if the median D_k of rare karyotypes (m < 1% and n₁·m ≥ 0.5) is ≥1.5 → R1 enters phase 1; otherwise R1 is skipped.
  - **Rule G2**: if the median of F_k / F_pred (computed from c and the passage records) over rare karyotypes is >2 → flag "jackpot excess", and allow R5 as a fallback for R1.

**D3** The e1 datasets with 200 cells, fit with base at ×1 and ×4. Data term = −Σ log Poisson(y | expected) (two time points); the rest = objective − data term. Report by grid cell the shares of the two parts in gain4. Report only, no action.

**D2** Upper bound on ranking at distance 3–4 for the development GRF landscapes (dev split): at each wavelength, estimate the covariance function of the GRF family (by L2 distance) by Monte Carlo, condition on the true log fitness of the karyotypes inside the base fit's domain, and take the Spearman of the GP posterior mean for the distance 3–4 karyotypes.
  - **Rule**: if KARST's distance 3–4 Spearman at both λ = 0.4 and 0.8 is ≥ 0.8 of the upper bound → R3 is not judged on ranking, only on bias and RMSE; otherwise R3's ranking criterion is a pooled distance 3–4 Spearman gain ≥0.05 with a landscape bootstrap 95% interval excluding 0.

**D5** First 10 datasets of the development set: the ratio of the latent-variable sd of non-major karyotypes between tmbstan (NUTS, hyperparameters fixed at ψ̂) and Laplace. Report only; if the median ratio is <0.8, note that Laplace is one of the causes of under-coverage.

**D1** Implement in dev two classes of landscape, RMF and integer-lattice NK (K = 2) (hash-determined independent noise, 22 chromosomes), 6 landscapes each × the 8 scenarios of the development benchmark, and do base fits.
  - **Rule**: for each family, the median n_major, the median fraction of cells inside the domain, and the median number of distinct karyotypes all fall within the 10–90% quantile range of the GRF development set; otherwise adjust the parameters once, and if they still do not match → stop and report.

## Phase 1: observation model

**R1** (NB1, second time point, φ_i = (N₁−1)·F·(1−p̂_i), F = 1 − Π_k(1 − c/B_k), where B_k is the number of survivors at each passage between the first and last time points and p̂_i is the observed frequency at the second time point; without passage records it is bit-for-bit identical to base, confirmed by unit tests)
- Fits: e1–e3 ends at ×1 and ×4; e4 at ×1 and ×4; the corresponding base fits under the same build.
- **To pass, all of the following must hold**:
  1. e3: the fraction of karyotypes with |ρ|/s_model > 4.5 does not rise in any of the 24 grid cells, and is at least halved for the bottleneck-200 grid cells pooled.
  2. e1: the R1/base RMSE ratio is ≤1.03 in ≥18 of the 24 grid cells at each distance, with median ≤1.02; e2 median ≤1.03.
  3. The test power on e1 at effect 0.05 (each method using its own threshold for 5% false positives on e3) drops by no more than 0.02.
  4. In grid cells with median φ ≥ 1, the fraction of datasets with gain4 > 0 is at least halved, and no grid cell rises. (Prediction recorded in advance: grid cells with φ < 0.2 do not change.)
  5. After refitting the calibration constants on the development GRF, the 95% coverage of e1–e3 is ≥0.95 at every distance.
- **Decision**: pass → recommend adoption as an optional argument. Not passed and G2 flagged → run R5 (same criteria). Not passed and not flagged → the R1 line ends.

**R2**
- (i) Hyperparameter uncertainty: J·V(ψ̂)·J′, where J is obtained by finite differences of the inner mode with respect to ψ, and V is the inverse of the Hessian of the Laplace objective with respect to ψ.
- (iii) Missegregation-rate mixture: β ∈ 5e-5 × {1, 4, 16}, equal weights on the log scale; the point estimate takes the lowest rate (consistent with `fit_sparse_joint_rates()`), and the variance is the mixture variance.
- Data: development and held-out GRF, e1–e4 ends, ALFA-K missegregation simulations (w2/w4/w8 ends), the RMF of D1.
- **Retention rule**: each item, added on its own and with the calibration constants refit on the development GRF, satisfies the general coverage criterion, and either its median width is ≥5% lower than that of base (with constants refit likewise), or base does not satisfy the criterion while the item does.

**R4** (miscall observation kernel, ε known)
- Data: the cloud_level miscall simulations (ε at the simulated value), clean development simulations (ε = 0.5%), ALFA-K w2 ABM (ε = 0.5%).
- **To pass, all of the following must hold**: on the miscall simulations, the error ratio relative to base (balance) is ≤0.85 at distances 0–2, and the single-copy direction accuracy at distance 1 is ≥0.55; on the clean development simulations the error ratio is ≤1.02; the median landscape r on the ALFA-K w2 ABM is ≥0.73.

## Phase 2: far-distance prediction

**R3** (shell-aggregated likelihood) and **R9** (mean-reverting decay + link function, used only together with R3)
- Tune first on the development GRF; the held-out GRF and RMF are evaluated only once.
- **To pass, all of the following must hold**: at λ = 0.4 and 0.8, the mean error (bias) at distances 3 and 4 drops by ≥40% (the baseline at distance 4 is 0.123, 0.071); in no (λ, distance) grid cell at distances 0–4 is the RMSE ratio worse than base by more than 0.02; coverage ≥0.90 in each grid cell; the same holds on RMF; the ranking criterion follows the D2 rule.
- R9 additionally needs: a pooled distance 3–4 Spearman gain ≥0.05 with an interval excluding 0, and no worse RMSE on GRF and on the additive family.

## Phase 3

**R5** (random rare-karyotype cloud): if it was not run as a fallback, apply R1's five criteria relative to R1 (relative to base if R1 did not pass), and recommend it only if it is better on criterion 1 or 4 and fails none of the others. The young-cloud part additionally needs: the median of step 0's log2(selected/true) within ±0.5, the fraction within 4-fold ≥0.93, and unchanged e1–e3 null results.

**R6**
- β prior (used with R1; with base if R1 did not pass): on e1–e3, a prior centred at ×1 with log sd = log 4 passes the step 2 criteria (≥18 of the 24 grid cells within 1.03, median ≤1.02, e2/e3 median ≤1.03); on e4, with a prior that contains the true value, the median RMSE ratio at distance 1 is <1.
- Dosage prior: on the held-out GRF, the RMSE ratio changes by ≤ +0.01 and Spearman by ≥ −0.02 in each grid cell, and the κ interval contains 0 in ≥90% of datasets; on a score landscape built from another source, the distance 3–4 Spearman gain is ≥0.10, the bootstrap interval excludes 0, and coverage is ≥0.90. If the required published score table cannot be obtained → not evaluable.

**R7** (history anchored by passage records): on the ALFA-K ABM one-passage datasets, median landscape r ≥0.75 and cross-validation ≥0.66; on patient_regrowth s0020/s0200 with N_b given, the skill at distances 1–2 is ≥0.05 higher than base, with a loss on observed major karyotypes ≤0.02; error ratio on the development scenarios ≤1.02.

**R8**
- Per-dataset parametric bootstrap (50 refits per dataset, 2 fixed-point rounds): satisfies the general coverage criterion on 48 datasets each of the development GRF and RMF, with median width ≤ 1.25 times that of base.
- Cross-family conformal (GRF, RMF, NK, leave one family out): satisfies the general coverage criterion on the held-out family.

**R10**
- WGD: on an ABM with WGD, the error on doubled karyotypes and their neighbours ≤ 0.85 of base; on an ABM without WGD ≤1.01; the per-cell horizon log score on the SA906 WGD trajectories is not below base.
- Piecewise fitness: in ≥2/3 of the Salehi PDX drug-withdrawal branches, the per-cell log score is better than with an unchanged landscape; error ratio on an ABM with a switching landscape ≤0.9. If the required treatment schedule cannot be obtained → not evaluable.

## When to stop and ask the user

A build fails for an unknown reason; more than 5% of the tasks of an array fail (NODE_FAIL is requeued by Slurm itself and does not count); G0 triggers; a criterion is unreasonable; a situation the rules do not cover.

## Clarification 2 (fixed at 2026-11-19 02:55 EST, before any R1 result was available)

Computation details, fixed under the reading closest to the original text:
- **Grid cells and ratios (criterion 2)**: same as step 2 (`beta_choice/analyse_series.R`). e1's grid cells = bottleneck × interval × cell count × number of time points (24 grid cells); each dataset's RMSE ratio (R1/base) at d0, d1, d2; the median within each grid cell; a grid cell "within 1.03" means that the largest of its grid-cell medians over the three distances is ≤1.03; "median ≤1.02" means that at each distance the median of the 24 grid-cell medians is ≤1.02. e2 and e3: at each distance, the median of the dataset ratios is ≤1.03.
- **Criterion 1**: in each grid cell, pool all scored karyotypes of all datasets in that grid cell and compute the fraction with |ρ|/s_model > 4.5; only R1 ≤ base counts as "does not rise" (no tolerance). For the 8 bottleneck-200 grid cells pooled, R1 ≤ 0.5 × base.
- **Criterion 3**: each method's threshold z* = the smallest z on the grid `sm_z_grid` at which the false-positive fraction of the pooled e3 null karyotypes (stratum all, kind null) is ≤0.05; power = the fraction of pooled e1 (stratum all, kind effect_0.05) that exceed the threshold at z* with the correct direction. R1 ≥ base − 0.02.
- **Criterion 4**: gain4 = objective(×1) − objective(×4) (the truth for e1 is ×1, so gain4 > 0 means favouring a too-high β). Each dataset's φ = (N₁−1)·F; a grid cell's φ = the median φ of the datasets in it.
- **Criterion 5**: each method refits (floor, deficit) with the Gaussian likelihood of `calibration_fit.R`, using its own fits on the development GRF (dev split, scoreboard datasets); the deficit for unobserved karyotypes = deficit × 1.1/0.6; R1's fits on the development GRF use the passage records of the development benchmark (obtained by replay). The refit constants are applied to the ×1 fits of e1, e2, e3: the coverage of each dataset at each distance is averaged, and it passes only if e1, e2 and e3 are each ≥0.95 at d0, d1 and d2.
- **R2's "median width"**: the median width of the 95% intervals of the scored karyotypes, pooled over datasets.

## Clarification 3 (fixed at 2026-11-19 03:09 EST, before any R4 or R2 result was available)

- **R4 miscall simulations**: all `call1`/`call2` datasets in the three cloud_level groups `diverse`, `diverse2`, `paper3` (108 in total); R4 uses each dataset's true `call_error`. "Error ratio ≤0.85" means that at each distance (d0, d1, d2) the median over datasets of the R4/base RMSE ratio is ≤0.85; "single-copy direction accuracy" is the dataset median of the scoring metric `edge_sign_accuracy` at d1.
- **R4 clean development simulations**: the core datasets of the dev split of the development GRF, ε = 0.005, median R4/base RMSE ratio ≤1.02 at each distance.
- **The general coverage criterion for R2 and R8** is computed with `dev/R/calibration.R` (each dataset is first averaged within a category, then averaged over datasets); "stress grid cells" are the stress scenarios of the development GRF and e4.

## Clarification 4 (fixed at 2026-11-19 03:36 EST, before any R7 fit)

- **Implementation of R7**: the history of the balance prior at the first time point is fixed to the population age in the records (days from establishment to the first time point, capped at 3000 days), instead of being chosen among candidates by the marginal likelihood. Reason: an unbiased bottleneck does not change the mean fraction of the cloud (Simons 2014), so on average the cloud has been accumulating since the population was established. Ages per set: 1440 days for the development GRF, e1, e3, e4 and patient_regrowth; the window start for e2 (240 or 480 days); for the ALFA-K ABM, the day of the first time point (the interval, if it is less than one interval).
- **Verdict for R7**: on the ALFA-K ABM one-passage set (w2_ends), the median landscape r over all minobs rows pooled is ≥0.75, and the median cross-validation r against the truth (`r_heldout_true`) is ≥0.66; in the s0020 and s0200 scenarios of patient_regrowth, the skill (dataset median) at distances 1 and 2 is each ≥0.05 higher than base, and the skill drop for karyotypes with ≥5 cells at distance 0 is ≤0.02; the median RMSE ratio is ≤1.02 at each distance in the core scenarios of the development GRF.

## Clarification 5 (fixed at 2026-11-19 03:37 EST, before any R8 computation)

- **R8(b) cross-family conformal**: the families are GRF (dev + heldout scoreboard of the development benchmark), RMF and NK, all fit with base. Score s = |prediction − truth| / model_sd (scored karyotypes, distances 0–2). Categories: observed <5, observed ≥5, unobserved d1, unobserved d2. When one family is held out, take the 95% quantile of s per category on each of the other two families (each dataset weighted equally within each category), and take the larger of the two families as q_c; the interval is prediction ± q_c·model_sd. The three families are held out in turn, and each time the general coverage criterion must be satisfied on the held-out family (grid cell = category × scenario, and for GRF also × wavelength), with a median width ≤ 1.25 times that of the family with the current constants (0.009, 0.6, 1.1). It passes only if all three are satisfied.
- **R8(a) per-dataset parametric bootstrap**: the simulator is the dev ABM, and the landscape is the fit result (posterior-mean birth rates inside the domain, the trend outside); start from the fitted abundances at the first time point (population size 2e6 as in the development benchmark, with the recorded passages), and sample with the recorded n0, n1, dt; 50 simulations and refits per dataset give the empirical distribution of error / model_sd, and the 95% quantile per category is taken as that dataset's q_c; the second round uses the sd calibrated in the first round to generate the dispersion of the truth (fixed point). Evaluated on the 48 core datasets of the development GRF (4 scenarios per landscape) and on 48 RMF datasets.

## Time correction (2026-11-19 08:26 EST)

The writing times originally stated for this file and for Clarifications 2–5 were estimates and not accurate. They have been changed to the actual writing times according to the session log: this file 02:45 → 02:36, Clarification 2 04:05 → 02:55, Clarification 3 05:05 → 03:09, Clarification 4 06:30 → 03:36, Clarification 5 06:45 → 03:37. Only the times were changed; the content of the criteria was not. Checked against the HPC job records (`results/jobs.tsv`), the ordering stated in each clarification still holds:
- Clarification 2 was written before any R1 result was available: R1's e-set fits were submitted at 02:51, and R1's development GRF fits at 03:00.
- Clarification 3 was written before any R4 or R2 result was available: R4 was submitted at 03:20; R2(i) was submitted at 03:05 and had no results yet at that time.
- Clarification 4 was written before any R7 fit: R7 was submitted at 03:37.
- Clarification 5 was written before any R8 computation: R8(a) was submitted at 03:45; R8(b)'s scoring script was run only after it was written at 03:38.

