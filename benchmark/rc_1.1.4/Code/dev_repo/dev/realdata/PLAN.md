# Real-data plan: can KARST's limitations be overcome, judged from the features of real data

Confirmed by the user 2026-11-19 14:25 EST ("Confirmed, run A–B automatically, stop after B"). Phases A and B run automatically; after B ends, stop and the user chooses the direction; C–E are confirmed separately.

## The user's rules

- **Add the diploid back.** For the inputs of the four NUMBAT cancer types (GBM, HNSC, Ovarian, Melanoma), add back the count of the fully diploid karyotype before handing them to KARST. PANcanKFLs's `ALFAK_input_preparation.R` removed it; the count is taken from `karyo_in` of `karyotype_count/<sample>.RObj` of the same cancer type.
- **Do not change any file of PANcanKFLs.** PANcanKFLs still feeds alfak by default. When its processing logic is needed, call it read-only or rewrite it in dev/.
- **The GBM karyotypes** are scDNA data processed with NUMBAT; there is no problem with the source.
- **Data roles** as in the plan:
  - GBM and Salehi are target data, used only for the final evaluation, not for tuning;
  - Ovarian is used for development;
  - HNSC and Melanoma are used for validation.
- **Use the current latest KARST** (cd342db, 1.1.2). Build it separately into HPC `karst_dev/lib/karst_cd342db`, and do not touch any operation of the "Karst bug fix" session, including the HPC karst clone and its jobs.
- **Everything goes in `dev/realdata/`**; no commits during the process.
- **HPC**: use xxlarge, no task limit, an add-hours monitor with every array, one fit per task.
- **Patient data** are read only in place on HPC and locally, not copied into the repository; reports show only summary numbers (summary statistics such as the number of cells and of karyotypes per dataset may appear; karyotype-level counts do not).

## Input

- **NUMBAT cancer types ("added-back" input).** Take `x` of `processed/ALFA_K/<cancer type>/ALFA-K_inputs/<sample>.Rds`, then add the diploid row.
  - Each column is matched to the column of `karyo_in` whose non-diploid counts are exactly the same. If none or more than one is found, the dataset is marked "cannot be added back" and does not enter B.
  - Times are taken from the column names (days). For those with more than two time points (Melanoma), take the first and the last.
  - Ovarian EOC153 has three versions, treated as the same patient: in cross-patient tests, when one of them is held out, the other two also stay out of the training set.
- **NUMBAT cancer types ("as-is" input).** That is, the input PANcanKFLs gives to alfak, without the diploid. Used only for the diploid sensitivity in B1.
- **Salehi.** The input of the paper benchmark (`processed/salehi/alfak_inputs`, read-only); use `salehi_input(ends = TRUE)` of `benchmark/paper` to take the first and last passages, with the same settings as the paper benchmark (beta 5e-5, `min_cn = 0`, n0 1e5, nb 1e7).
- **Salehi development/confirmation split**: written into this file after A3 is done and before B starts. Split only by the size from A3 (number of trajectories per cell line, treatment status), without looking at any fit result.
- **Salehi split (fixed 11-19 14:39, B has no result yet).** A3: 277 trajectories, SA609 135 (untreated 93, treated 42), SA906 49, SA535 39 (24/15), SA532 28, SA1035 20 (10/10), SA039 6.
  - **Development**: SA609 and SA039, 141 in total. It has treated and untreated, PDX, and also a p53 wild-type cell line.
  - **Confirmation**: SA906, SA532, SA535, SA1035, 136 in total. It includes a p53-knockout cell line, an untreated PDX, and two treated PDXs.
  - B is a diagnosis; all cell lines are run, and the report lists the two groups separately.

## Fit settings

- **NUMBAT cancer types**, as in `validation/sparse_joint/real_cohort.R`:
  - the division rate is unknown; take three rates in 0.02–1 /day (centre 0.141), with `fit_sparse_joint_rates`;
  - beta 5e-5, `min_cn` = min(1, smallest observed value), `max_cn` = largest observed value + 1;
  - otherwise the package defaults (balance prior, history chosen by marginal likelihood).
- **Salehi**: the paper benchmark's `fit_karst()`, writing the alfak format (minobs 5/10/20).
- Build: cd342db, no patches.

## Phase A: data profile (read-only)

- **A1 Size and structure.** For each input (after adding back), report:
  - the number of cells and the diploid fraction at each time point;
  - the number of karyotypes, and the fraction of cells whose karyotype is seen in a single cell;
  - the number of chromosomes that vary, mean ploidy;
  - the fraction of cells covered by the KARST model domain (`diagnostics$domain_coverage` of cd342db, taken from the B1 fits, added in the B report);
  - metadata only as summaries: treatment, sampling site, time interval.
- **A2 NUMBAT measurement process.** Take each patient's `joint_post` from the round in which NUMBAT converged (the same as PANcanKFLs's `getCNVmatrix`), and report:
  - the number of consensus segments, the number of chromosomes with consensus segments (all other chromosomes are recorded as 2 and cannot be seen); the number of chromosomes with several segments (arm-level changes are merged).
  - The distribution of p_cnv over cell × segment: how much is ≥ 0.9 (called a change), 0.5–0.9 (uncertain, recorded as unchanged), < 0.5.
  - **Miss rate given by NUMBAT's own posterior** = Σ_{p<0.9} p_cnv / Σ p_cnv (the fraction of expected change events recorded as "unchanged"); **false-call rate** = Σ_{p≥0.9}(1 − p_cnv) / number called as changed.
  - **Reproduction check**: in dev/, recompute each cell's karyotype with the rules of `getCNVmatrix` + `NumbatPostProcess`, and compare with `karyo_in`. When they agree completely, also report the fraction of cells whose karyotype changes when the threshold is set to 0.5; when they do not agree, report only the segment-level numbers.
  - The number of diploid cells removed at each time point (now added back).
- **A3 Salehi.** The same size measures, plus:
  - the number of trajectories and the treatment status per cell line;
  - passage intervals;
  - the fraction of suspected WGD cells (cells with mean copy number ≥ 2.8).
- **Output**: `A/report.md`.

## Phase B: diagnosing KARST on real data (no truth)

**B1 Current state**, all inputs:
- fit status (with error class) and run time;
- Pearson dispersion at the second time point (definition of `real_cohort.R`, 1 = Poisson);
- **zero-count neighbour check**: for non-major karyotypes in the model domain, the ratio of the number observed at the second time point to the number expected under Poisson(expected count), Σ(1 − e^{−expected});
- rank agreement between different division rates (NUMBAT cancer types).

**B1 Sensitivity**, at the centre rate (Salehi with its default settings). Each pair of fits is compared on the observed karyotypes that both have (t0 + t1 count > 0), reporting:
- the Spearman of relative_fitness;
- the fraction with |Δ|/√(sd₁² + sd₂²) > 1 (Δ first has the mean difference of the two on karyotypes with ≥10 cells subtracted, because the two fits can have different references);
- the same two numbers for majors only.

Pairs of fits compared:
- beta ×1/4, ×4 against ×1;
- `initial_state = "free"` against the default balance;
- diploid "as-is" (removed) against "added back", NUMBAT cancer types only.

**B2 Salehi:**
- **(a) Log score on later passages.** Follow the procedure of `validation/sparse_joint/salehi_horizon/forecast_score.R`, called as is, only replacing the fits with this round's cd342db fits:
  - the per-cell log score gain of the prediction over "no change" (persistence), stratified by days (≤15, 16–30, 31–60, >60);
  - the fraction of cells of later passages that fall in the landscape (coverage).
- **(b) Coverage of prediction intervals.** Only minobs 5 and later passages ≤60 days:
  - draw 20 independent sets of fitness from the landscape's mean ± sd, predict with the ODE of `alfakR::predict_evo`, then draw multinomially with that passage's number of cells, to get a 95% prediction interval for each karyotype frequency;
  - for karyotypes in the landscape that are observed at that passage, or have predicted frequency ≥ 1/number of cells, count the fraction whose observed frequency falls in the interval.
  - This item ignores the correlation between fitness values and is for reference only.
- **(c) Agreement of parallel lineages.** Compare pairwise the trajectories of the same cell line with the same treatment status and non-overlapping input passages. On karyotypes observed in both summaries, compute z = (Δ − mean Δ)/√(sd₁² + sd₂²), and report the fraction with |z| > 1.96.

**B3 NUMBAT cancer types:**
- **(a) Cross-patient prediction (leave-one-out).** Within a cancer type (GBM, Ovarian, HNSC; Melanoma has only one eligible patient, not done):
  - fit a cohort with the added-back inputs of the other patients (`fit_sparse_joint_cohort`, sharing = "partial", other defaults, centre rate);
  - score the held-out patient's karyotypes with the shared landscape (`predict(..., sample = NULL)`);
  - on karyotypes with total count ≥ 3 in the held-out patient, compute the Spearman with the observed log ratio r = log((y₁+0.5)/N₁) − log((y₀+0.5)/N₀);
  - reference: the Spearman of the relative_fitness from the same patient's own fit with r (in-sample, upper bound);
  - the mean across patients and a bootstrap 95% interval (resampling patients).
- **(b) Split-half agreement.** For datasets with ≥200 cells at both time points (NUMBAT and Salehi), as in `real_cohort.R`, split the cells randomly into two halves, fit each, and report:
  - the Spearman;
  - the fraction of differences that fall in the reported 95% interval ("within_reported");
  - in three groups: 10+, 2–9, 1 cell.

## How to read B (fixed in advance)

B is only a diagnosis and does not judge pass or fail. Each item below, if met, is marked as "active on real data":

| ID | Limitation | Marking condition |
|---|---|---|
| L1 | Fit failure | Failure rate > 5% in any data group |
| L2 | Observation model (NUMBAT) | Any of: median diploid-sensitivity Spearman < 0.9; median fraction with |Δ|/sd > 1 above 0.2; median zero-count neighbour ratio outside [0.5, 2]; median A2 miss rate > 0.1 |
| L3 | Intervals insufficient on real data | Any of: median within_reported of the 10+ group in split-half agreement < 0.90; fraction of Salehi parallel lineages with |z| > 1.96 above 0.10; B2(b) coverage < 0.85 |
| L4 | Sensitive to β or history | Median Spearman < 0.9 between β ×1/4 and ×4, or between free and balance |
| L5 | Distant prediction | In the Salehi strata >15 days, gain over persistence not greater than 0 (interval contains 0 or is negative); or median coverage of later passages < 0.7 |
| L6 | Overdispersion (drift or founder effects) | Median Pearson dispersion at the second time point > 2 (in either the Salehi or a NUMBAT group). This is the check for "evidence of drift": only if this is marked is a drift term reconsidered |
| L7 | Cross-patient information | Mean leave-one-out Spearman > 0 and interval excludes 0: supports D2; otherwise D2 has no basis |
| L8 | WGD | Fraction of suspected WGD cells > 0.1 in any Salehi cell line |
| L9 | Too little data | Report only: datasets with a time point of < 30 cells |

The B report (`B/report.md`) lists the evidence and the mark for each item, ordered by impact, **and then stops**.

## Stopping rules

- Stop and report when the build fails for unknown reasons, or when > 5% of the tasks of an array fail and this is not a cost of the method itself (for example memory, which can be rerun).
- Stop and report when adding back the diploid fails for more than 2 datasets.
- Once this file is fixed, the reading is not changed; corrections are written in the "Corrections" section at the end, with time and reason.

# Phases C–E (confirmed by the user after 2026-11-20 00:29 EST: "proceed as you suggest")

This section was written at 2026-11-20 00:35 EST, before any fit and scoring of C–E. It is not changed afterwards; corrections are written in "Corrections" at the end.

## Scope

- **Do:** D4 (intervals), O (overdispersed likelihood), D3 (WGD).
- **Not done:**
  - D1: decided by the user 2026-11-20: KARST by default takes the input karyotypes as correct, and does not model for a particular copy-number software.
  - D2: L7 was not marked.
- **D5 not done for now, left for the user to decide.** Before writing the plan, the B fits were checked:
  - the 24 failures of SA609 have no karyotype shared between the two passages, and are a median of 11–13 steps from the nearest major;
  - in successfully fitted datasets, cells outside the domain are a median of 3–5 steps from the nearest major (the 10%–90% quantiles for SA609 are 3–5).
  - Extending the neighbours by only one layer (distance 2) basically cannot recover these cells, so there is no local change of the model domain to try.
- **Version:** patches are based on cd342db, placed in `dev/patches/`, and built into a separate library under HPC `karst_dev/lib/`. The package itself is not changed.
- **Other rules as before:**
  - no commits;
  - HPC with xxlarge, no limit;
  - add-hours monitor with the small QOS, 3 days;
  - patient data only as summaries;
  - while submitted tasks are running, do not change the scripts they use.

## Data roles (carried over from A–B)

- **For development:** Salehi development group (SA609, SA039), Ovarian.
- **Only for the final verdict:** Salehi confirmation group (SA906, SA532, SA535, SA1035), GBM.
- **Validation:** HNSC, Melanoma PatientB.
- **Simulations:** the dev split of the development GRF is used for development; held-out GRF, RMF, NK, e1–e3 only for verdicts.

## C: simulations matched to Salehi (C2)

- **Purpose:** provide truth in a Salehi-like setting, for O, D3, D4.
- **Simulator:** dev ABM (`dev/abm/abm_dev.cpp`); `abm_wgd.cpp` for the WGD scenario. GRF landscapes, wavelengths 0.8 and 1.2, 6 landscapes each.
- **Design:**
  1. start from the diploid, and obtain a diverse population through a period of pre-evolution;
  2. then passage every 15 days, diluting to bottleneck B;
  3. sample n₀ and n₁ cells at two passages, with intervals of 15, 30, 45, 60 days;
  4. the true frequencies of later passages are kept, for distant tests.
- **Tunable:** the ABM missegregation rate, bottleneck B, the length of pre-evolution, and the gain g in the WGD scenario. At most two grid rounds per scenario.
- **Acceptance for the three scenarios** (the median over simulated datasets must fall within the interval):

  | Scenario | Corresponds to | Singleton-karyotype fraction at the last passage | Model-domain coverage (worst time point, cd342db base fit) | Pearson dispersion at the second time point |
  |---|---|---|---|---|
  | H | SA609 | [0.34, 1.0] | [0.10, 0.40] | [0.86, 2.5] |
  | M | Confirmation group | [0.10, 0.50] | [0.55, 0.92] | [2.0, 6.0] |

  - **W (WGD, corresponding to SA906a):** the median fraction of ploidy-4 cells at the first time point is in [0.02, 0.5], and in ≥ 25% of datasets the ploidy-4 fraction rises by ≥ 0.1 between the two time points.
  - **The tuning order for W is fixed in advance** (the lesson of R10): first g = 0.005, then g = 0; ω = 1e-4 unchanged.
- **Still not acceptable after two rounds:** stop and report. Criteria that depend on that scenario are not automatically judged as failed.

## D and E: candidates and verdicts

### O: overdispersed likelihood (patch `O_nb1_est.patch`, two time points only)

**Method:** two steps.
1. First do the base fit. On model-domain states with expected count ≥ 5, compute the Pearson dispersion D₅ at the second time point.
2. Take φ̂ = max(0, D₅ − 1); when fewer than 3 states are usable, φ̂ = 0.
3. Refit the counts of the second time point with NB1 (variance μ(1 + φ̂)), and multiply the posterior information by 1/(1 + φ̂), as in R1.

With φ̂ = 0 it is bit-identical to base. There is no tunable parameter.

**Difference from R1:** φ is estimated from the data, not computed from passage records. So on data with a large bottleneck, close to Poisson, it falls back to base automatically.

**A pass requires all of:**
- **O-1 No worse (e group):**
  - in e1, ≥ 18 of 24 grid cells have an RMSE ratio (O/base) ≤ 1.03 at d0–d2;
  - the median over grid cells at each distance ≤ 1.02;
  - for e2, e3, the median at each distance ≤ 1.03;
  - for the held-out part of the development GRF and for RMF, the median at each distance ≤ 1.02.
- **O-2 False positives and power (e3, e1, algorithm as in Clarification 2):**
  - with bottleneck 200 pooled, the fraction with |ρ|/s > 4.5 is not higher than base;
  - power at the 5% false-positive threshold ≥ base − 0.02.
- **O-3 C2 truth:** in scenarios H and M, the median RMSE ratio at d0–d1 ≤ 1.00.
- **O-4 Real data, parallel lineages (Salehi confirmation group):**
  - the gain in pooled Spearman ≥ 0.03, with a bootstrap over lineage pairs whose interval excludes 0;
  - the fraction with |z| > 1.96 using the model sd drops by ≥ 0.05.
- **O-5 Real data, later passages (Salehi confirmation group, ≤ 60 days):** the paired difference of the per-cell gain against base ≥ −0.05.
- **O-6 GBM split halves:** with the model sd, the agreement fraction of the 10+ cell group is not lower than base by more than 0.02.

### D3: WGD (patch `D3_wgd.patch`, i.e. R10_wgd moved onto cd342db, content unchanged)

A pass requires all of the following, with W1–W3 defined as in `dev/phase3/r10_spec.md`; both ω ∈ {1e-4, 1e-3} must meet them:
- **W1:** on C2 scenario W, for scored karyotypes with mean copy number ≥ 3, the median RMSE ratio ≤ 0.85.
- **W2:** on the development GRF without WGD (dev split), the median RMSE ratio at d0–d2 ≤ 1.01; ≤ 1.01 on C2 scenarios H and M.
- **W3:** on the eligible ends trajectories of SA906a, the paired mean of the per-cell log score (`ll_fc_e2`) ≥ base. base uses cd342db.

### D4: intervals (computed only after the fit, likelihood unchanged)

**Problem:** in the current form deficit 0.6 > 1/1.96, so the 95% interval by construction always contains 0, and therefore carries no information.

**Development** (only the dev split of the development GRF and the Salehi development group):
1. First describe the error structure: how the sign accuracy varies with |estimate| / model sd, and the distribution of the ratio of truth to estimate.
2. Choose one of the three forms below:
   - F1: the current form, with the constants refitted;
   - F2: a sign-preserving multiplicative form, i.e. the |rel| term widens only in the direction away from 0, with the lower and upper sides calibrated separately;
   - F3: conformal by class, i.e. q_c × √(model sd² + floor²(1 + shell)²), with classes as in Clarification 5.
3. **Selection rule:** among the forms that meet "pooled coverage per class on the development GRF ≥ 0.95, and fraction of Salehi development-group parallel lineages with |z| > 1.96 ≤ 0.07", take the one with the highest "fraction of karyotypes with ≥ 10 cells whose interval excludes 0" in the Salehi development group.

**Verdict** (evaluated only once):
- **D4-1 Simulation coverage:** on held-out GRF and RMF, pooled coverage per observed class ≥ 0.93, and unobserved d1, d2 ≥ 0.90. The same holds on e1–e3. NK is reported only, and is required not to be lower than the current constants by more than 0.02.
- **D4-2 Parallel lineages:** the fraction with |z| > 1.96 in the Salehi confirmation group ≤ 0.10.
- **D4-3 Split halves:** in the 10+ cell group of GBM and Ovarian, the fraction falling inside the interval ≥ 0.90.
- **D4-4 Informative:**
  - in the Salehi confirmation group and GBM, the fraction of karyotypes with ≥ 10 cells whose interval excludes 0 has a median ≥ 0.20;
  - on held-out simulations, among karyotypes whose interval excludes 0, the fraction with the correct sign ≥ 0.95.
- **Report only:** frequency prediction intervals (propagating the sampling error of the starting passage as well). They are alfakR's prediction, not KARST's output.

**Combination:** if O passes, the form chosen for D4 has its constants refitted on the O fits by the same steps, and is judged once more with the same criteria.

### Stopping rules

- The build fails for unknown reasons;
- > 5% of the tasks of an array fail, and not because of resources;
- C2 acceptance fails in two rounds;
- a case the rules do not cover.

In any of these cases, stop and report. After all verdicts, write `CE/report.md`, and stop for the user to decide whether it goes into the package.

## Clarification C1: concrete definitions of the three D4 forms

Written 2026-11-20 00:43 EST, before any D4 computation on held-out data (held-out GRF, RMF, NK, e1–e3, Salehi confirmation group, GBM). Before this, the error structure was looked at only on the dev split of the development GRF: measured with the model sd, the sign accuracy is 0.93–0.97 when |estimate| / model sd is 2–4, and 0.98–0.99 when it is 4–8; on karyotypes with |z| > 2, the median ratio of truth to estimate is 1.05–1.27.

**Notation**
- z = estimate / model sd.
- u = √(model sd² + (0.009 · division rate · (1 + shell))²).
- **Class c:** observed, total count < 5; observed, ≥ 5; unobserved, distance 1; unobserved, distance 2. Checks on real data use only observed karyotypes.
- **|z| bin b:** [0, 1], (1, 2], (2, 4], (4, 8], > 8. A bin with fewer than 200 karyotypes in the development data is merged into the next bin.
- All quantiles are computed with weights within a class (for F2, within class × bin), with the same weight for each dataset.

**The three forms**
- **F1:** the current form. Refit (floor, deficit) on development karyotypes with `dev_calib_fit` (`dev/R/calibration.R`), with the unobserved deficit multiplied by 1.1/0.6. The interval is ±1.96 sd.
- **F2:** by class × |z| bin, calibrated separately on the side toward 0 and on the side away from 0.
  - Half-width on the side toward 0 = q_t · u, q_t = the 97.5% quantile of max(0, −error · sign(estimate)) / u;
  - half-width on the side away from 0 = q_a · u, q_a = the 97.5% quantile of max(0, error · sign(estimate)) / u.
- **F3:** symmetric conformal by class, half-width = q_c · u, q_c = the 95% quantile of |error| / u.

**z for parallel lineages:** when F2 is asymmetric, each estimate takes the half-width on the side facing the other estimate, divides it by 1.96 to use as the sd, and puts it into the z formula of L3.

**Selection:** by the rule in PLAN. The threshold "pooled coverage per class on the development GRF ≥ 0.95" is the coverage computed on the development data. The quantiles of F2 and F3 may be raised from 0.975 and 0.95 to at most 0.98 and 0.96 to reach this threshold, taking the smallest value that reaches it.

**O is paused** (see `CE/o_stop.md`), so the "Combination" item does not apply for now.

## User decision 1 (written 2026-11-20 01:15 EST, after the D4 selection result and before any held-out evaluation)

The user confirmed ("proceed as you suggest"), as follows:
- **D4 selection threshold relaxed:** per-class coverage on the development GRF changed from ≥ 0.95 to ≥ 0.945. The quantile search takes the smallest value that reaches the same threshold. Parallel lineages ≤ 0.07 unchanged; the held-out verdict criteria D4-1 to D4-4 unchanged, evaluated only once.
  - Under the new threshold, F1 (d2 0.939) and F3 (lineages 0.081) still fail, so only F2 is evaluated.
  - This is a change made after seeing the development results, decided by the user; the report states it as it is.
- **O′ not decided yet:** first refit C2's H (5e-4) and M (1.5e-4) with the simulated missegregation rate β_sim, to see whether base being "worse than predicting 0" on H is caused by a wrongly set β, then decide whether to do O′.

## User decision 2 (written 2026-11-20 01:25 EST)

The user: "Karyotype fitness landscapes are certainly not like RMF. So there is no need to do RMF."
- All later evaluations no longer use RMF. Any criterion that says "held-out GRF and RMF" looks only at GRF. NK is not affected.
- The existing verdicts have been rechecked on GRF only: R3 and R9 of phase 2, R8(a) and R8(b) of phase 3, and F2 of D4 each have a criterion not met on GRF (or on GRF and the e group), so the verdicts are unchanged with RMF removed.
