# Underestimation of out-of-domain karyotypes: plan and criteria fixed in advance

Confirmed by the user 2026-11-25 01:20 EST (D1a, D2a, D3a, D4b), with the requirements: all code goes in `karst/patch/`, not into git (`.git/info/exclude`), and existing KARST is not modified. Results on HPC go outside the KARST clone: `/share/lab_crd/taoli/Project/karst_patch/outside_domain/`.

## Problem

KARST's fit domain = major karyotypes (fraction ≥1% at some time point and ≥3 cells in total) + their one-step neighbours.
- Observed karyotypes outside the domain do not enter the likelihood; `predict()` gives them the average level of "observed non-major karyotypes inside the domain".
- Karyotypes never seen are extrapolated by the trend, lower the farther they are from the major karyotypes.

On codex's 44 Salehi daughters (retrospective):
- the out-of-domain fraction actually rose in 30/44, while KARST predicted a fall in 40/44;
- the "near out-of-domain" 2–4 steps from the major karyotypes actually rose in 32/44, while KARST predicted a rise in only 3/44;
- the change of group between the two training time points does not continue into the next passage.

## Definitions

- **Fit domain:** the default settings of `fit_sparse_joint()`. Salehi uses `min_cn = 0`, and `max_cn` takes max(8, largest copy number in the data).
- **Groups:** inside the domain; near out-of-domain (outside the domain, ≤4 steps from the nearest major karyotype); far out-of-domain (>4 steps). 3 steps and 6 steps as sensitivity analyses.
- **Karyotype table E:** karyotypes seen at the two training ends, plus their one-step and two-step neighbours (newly generated neighbours exclude 0 copies).
- **Forecast:** deterministic expectation propagation (`R/engine.R`), with rules consistent with the benchmark ABM (alfakR dev2 f93b2d2). Start from the full parent, step 0.1 days, and discard descendants outside E.

## Step 0: engine check

codex's 35 Salehi windows: the same E table, the same initial allocation, the same rate maps (KARST native, ALFA-K), compared with codex's average of 5 ABM runs.

- **Pass:** for each daughter, the median absolute difference of the Fig. 2c angle between deterministic and ABM ≤ 3°; the median absolute difference of the three group fractions ≤ 2 percentage points. The two rate maps are judged separately.
- **Fail:** stop, ask the user.

## Step 1: diagnosis on simulations with truth (no refit)

**Data:**
- fig2d_truth step 2: 400 main simulations × 3 designs, Salehi cell numbers, KARST 1.1.3, β = 5e-5 (i.e. the true value);
- fig2d_misseg: two groups, β fixed and β at the true value.

**A1 Fitness bias:** out-of-domain karyotypes seen in training, and their unseen one-step neighbours (not inside the domain). Bias = KARST relative fitness (`relative_fitness` of `predict()`) minus true relative fitness. The true values use the same reference as KARST as zero: the count-weighted average of the major karyotypes. Stratified by distance (2, 3–4, >4 steps), report the median bias and the Spearman correlation.

**A2 Forecast of group fractions:**
- Start from the parent sample used for the fit, and forecast to the second sampling (about 15 days) with the same engine and the same E table.
- Three kinds of fitness: the true landscape (with the fitted β), KARST, the true landscape (with the true β of the missegregation-rate simulation).
- The actual values use the samples of 1000 cells stored by the simulation.

**Branches (development half only):**

| Case | Verdict | Next step |
|---|---|---|
| Among inputs whose near out-of-domain fraction actually rose, KARST predicts a fall in fewer than half | Not reproduced | Stop |
| Reproduced, and the true landscape (with the fitted β) has the right direction in ≥ 70% of these inputs | Problem is in the fitness | Go to step 2 |
| Reproduced, but the true landscape has the right direction in < 70% | Problem is in the flux / karyotype table | Stop; the forecast rate is for the user to decide (β line) |

### Refinements before step 1 starts (2026-11-25 01:50, no step 1 result seen yet)

1. **Where KARST's division rate comes from:** native division rate = `fitness_mean / (2P_ii − 1)` from `predict()`, with P_ii from the fitted β, as codex does.
2. **Aligned KARST:** the native level is set by the prior birth rate (1/day), which differs from the level of the truth; with a different level, the flux of new karyotypes produced by missegregation also differs. So A2 adds an arm "aligned KARST": add one constant to all of KARST's division rates so that the parent's count-weighted mean division rate equals the value of the true landscape. Then it differs from the true landscape only in relative fitness.
   - "KARST predicts a fall" in the branch verdict uses the native version, as in the Salehi case; the aligned version is reported only.
   - C2 of step 2 uses the aligned version; the native version is reported only.
3. **Computed in one go:** to save job overhead, the candidate arms of step 2 (B1, B2, B1+B2; native and aligned) are computed in the same batch of jobs as step 1. The step 1 analysis script reads only the step 1 arms; the step 2 analysis is run only if the branch is judged "problem is in the fitness".
4. **True relative fitness:** true net growth = GRF value × (2P_ii − 1), with P_ii from the fitted β. Reference point = the average of the major karyotypes weighted by total counts over the two time points, the same as KARST's `reference`.

## Step 2: candidate changes (only the assignment at prediction time changes)

- **B1 Neutral level:** out-of-domain karyotypes seen in training get relative fitness = 0 (log division rate equal to the reference level). sd² = sd² of the original observed level + the variance of the log division rates of observed karyotypes inside the domain.
- **B2 Nearest inheritance:** for karyotypes neither inside the domain nor seen in training, if the nearest observed out-of-domain karyotype o is closer than the nearest major karyotype, log division rate = level(o) + [trend(k) − trend(o)].
  - level(o) takes KARST's existing observed level when B2 is used alone, and the neutral level when combined with B1.
  - The trend is the trend fitted by KARST (dose, curvature, shell distance).
  - All other karyotypes are unchanged.
- **B1 + B2.**

**Data split:** simulations are split into a development half and a confirmation half by odd/even rep, and judged on the confirmation half; the own benchmark uses the held-out set.

| | Criterion (confirmation half) |
|---|---|
| C1 | Median bias of observed out-of-domain karyotypes, in absolute value ≤ half the current one; drop in Spearman with the truth ≤ 0.02 |
| C2 | Absolute error of the change in near out-of-domain fraction, paired median ratio ≤ 0.8; W1 of the full distribution, paired median ratio ≤ 1.00 |
| C3 | Own benchmark held-out set: RMSE at distances 0, 1, 2 each worse by ≤ 2%; drop in 95% coverage at distance 0 ≤ 0.02 |

- A candidate passes only if it passes all three.
- **Salehi 44 daughters:** report only W1, overlap rate, angle, near out-of-domain direction, with no verdict.
- After a pass, stop and show the user. Merging into the package needs separate consent from the user, and the default stays unchanged.

## Rules

- The package is not changed, `dev/` is read-only, no commits.
- HPC: one input per task, xxlarge, exclusion list (including 0003), an add-hours monitor with every job.
- Locally at most 3 heavy processes at the same time.
- Budget ≤ 1,000 CPU hours; ask before exceeding it.
