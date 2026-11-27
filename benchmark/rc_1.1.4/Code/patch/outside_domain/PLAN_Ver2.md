# Underestimated out-of-domain karyotypes: plan and pre-registered criteria

English version of `PLAN.md` (Ver2, 2026-11-25). The content is unchanged; the
Chinese original is kept.

Confirmed by the user on 2026-11-25 at 01:20 EST (D1a, D2a, D3a, D4b), with
these conditions:
- all code goes in `karst/patch/`;
- that folder stays out of git (`.git/info/exclude`);
- existing KARST is not modified.

HPC results go outside the KARST clone:
`/share/lab_crd/taoli/Project/karst_patch/outside_domain/`.

## Problem

KARST's fitting domain = the major karyotypes (share ≥ 1% at some time point
and ≥ 3 cells in total) + their one-step neighbours.
- Observed karyotypes outside the domain do not enter the likelihood.
  `predict()` gives them the average level of the observed non-major
  karyotypes inside the domain.
- Never-observed karyotypes are extrapolated along the trend, lower the
  farther they are from the majors.

On codex's 44 Salehi daughters (retrospective):
- the out-of-domain share actually rose in 30/44, while KARST predicted a fall
  in 40/44;
- the "near out-of-domain" share (2–4 steps from a major) actually rose in
  32/44, while KARST predicted a rise in only 3/44;
- the change of the groups between the two training time points did not
  carry over to the next passage.

## Definitions

- **Fitting domain:** the default settings of `fit_sparse_joint()`. Salehi
  uses `min_cn = 0`, and `max_cn` = max(8, the largest copy number in the
  data).
- **Groups:**
  - inside the domain;
  - near out-of-domain ("halo": out of the domain, ≤ 4 steps from the nearest
    major);
  - far out-of-domain ("island": > 4 steps).

  3 and 6 steps are used as sensitivity analyses.
- **List E:** the karyotypes seen at either end of the training data, plus
  their one- and two-step neighbours. Newly generated neighbours never have 0
  copies.
- **Forecast:** deterministic propagation of expected values (`R/engine.R`),
  with the same rules as the benchmark ABM (alfakR dev2 f93b2d2).
  - It starts from the complete parent, with a step of 0.1 days.
  - Offspring outside E are dropped.

## Step 0: engine check

On codex's 35 Salehi windows, use the same list E, the same initial
allocation and the same rate maps (KARST native, ALFA-K), and compare against
the mean of codex's 5 ABM runs.

- **Pass:** the median absolute difference between the deterministic and ABM
  forecasts must be:
  - ≤ 3° for each daughter's Fig. 2c angle;
  - ≤ 2 percentage points for the three group shares.

  The two rate maps are judged separately.
- **Fail:** stop and ask the user.

## Step 1: diagnosis on simulations with known truth (no refitting)

**Data:**
- fig2d_truth step 2: 400 main simulations × 3 designs, at Salehi cell
  counts, KARST 1.1.3, β = 5e-5 (the truth);
- fig2d_misseg: the β-fixed and β-true sets.

**A1, fitness bias.** The objects are the out-of-domain karyotypes seen in
training, and their unseen one-step neighbours that are not in the domain.
- Bias = KARST's relative fitness (`relative_fitness` from `predict()`) minus
  the true relative fitness.
- The truth is zeroed on the same reference as KARST: the count-weighted
  mean of the majors.
- Results are stratified by distance (2, 3–4, > 4 steps), reporting the
  median bias and the Spearman correlation.

**A2, group-share forecasts.**
- Start from the parent sample used for fitting and forecast to the second
  sampling time (about 15 days), with the same engine and the same list E.
- Three fitness inputs:
  - the true landscape (at the fitted β);
  - KARST;
  - the true landscape at the true β of the missegregation simulations.
- The observed value is the 1000-cell sample stored by the simulation.

**Branch (development half only):**

| Situation | Verdict | Next |
|---|---|---|
| Among inputs whose near out-of-domain share actually rose, KARST predicts a fall in fewer than half | Not reproduced | Stop |
| Reproduced, and the true landscape (fitted β) gets the direction right in ≥ 70% of those inputs | The problem is fitness | Go to step 2 |
| Reproduced, but the true landscape gets the direction right in < 70% | The problem is flux / the list | Stop; the forecast rate is the user's decision (the β line) |

### Refinements before step 1 was run (2026-11-25 01:50, before any step-1 result had been seen)

1. **KARST's division rate:** native division rate = `fitness_mean / (2P_ii − 1)`
   from `predict()`, with P_ii at the fitted β, as codex does.
2. **Aligned KARST:** the native level is set by the prior birth rate
   (1/day) and differs from the truth's level. A different level also gives
   a different missegregation flux into new karyotypes.

   So A2 gets one more arm, "aligned KARST": a constant is added to all of
   KARST's division rates, so that the parent's count-weighted mean division
   rate equals that of the true landscape. It then differs from the true
   landscape only in relative fitness.
   - "KARST predicts a fall" in the branch uses the native version, as in the
     Salehi case; the aligned version is only reported.
   - C2 in step 2 uses the aligned version; the native version is only
     reported.
3. **Computed in one pass:** to save job overhead, the step-2 candidate arms
   (B1, B2, B1+B2; native and aligned) are computed in the same jobs as
   step 1.
   - The step-1 analysis script reads only the step-1 arms.
   - The step-2 analysis runs only if the branch is "the problem is fitness".
4. **True relative fitness:** true net growth = GRF value × (2P_ii − 1), with
   P_ii at the fitted β. The reference is the majors' mean weighted by their
   total counts over the two time points, the same as KARST's `reference`.

## Step 2: candidate rules (prediction-time assignment only)

- **B1, neutral level:** out-of-domain karyotypes seen in training get
  relative fitness = 0 (log division rate equal to the reference level).
  - sd² = the old observed level's sd² + the variance of the log division
    rates of the observed in-domain karyotypes.
- **B2, inherit from the nearest:** a karyotype k that is not in the domain
  and was not seen in training, if it is closer to the nearest observed
  out-of-domain karyotype o than to the nearest major, gets
  log division rate = level(o) + [trend(k) − trend(o)].
  - level(o) is KARST's current observed level when B2 is used alone, and the
    neutral level when combined with B1.
  - The trend is KARST's fitted trend (dosage, curvature, shell distance).
  - All other karyotypes are unchanged.
- **B1 + B2.**

**Data split:**
- simulations are split by replicate parity into a development half and a
  confirmation half, and judged on the confirmation half;
- the own benchmark uses its held-out set.

| | Criterion (confirmation half) |
|---|---|
| C1 | Median bias of the observed out-of-domain karyotypes: absolute value ≤ half the current one; Spearman with the truth falls by ≤ 0.02 |
| C2 | Absolute error of the near out-of-domain share change: paired median ratio ≤ 0.8; W1 of the full distribution: paired median ratio ≤ 1.00 |
| C3 | Own benchmark held-out set: RMSE at distances 0, 1 and 2 each worse by ≤ 2%; 95% coverage at distance 0 falls by ≤ 0.02 |

- A candidate passes only if it meets all three.
- **Salehi, 44 daughters:** only W1, overlap, angle and near out-of-domain
  direction are reported, with no verdict.
- After a pass, stop and show the user. Merging into the package needs the
  user's separate approval, and the default stays unchanged.

## Rules

- The package is not modified, `dev/` is read only, nothing is committed.
- HPC:
  - one input per task;
  - qos xxlarge;
  - the excluded-node list (including 0003);
  - an add-hours monitor for every job.
- At most 3 heavy local processes at a time.
- Budget ≤ 1,000 CPU-hours; ask before going over it.
