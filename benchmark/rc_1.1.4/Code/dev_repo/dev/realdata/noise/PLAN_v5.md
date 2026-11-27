# Pre-fit preprocessing v5: look only at frequent karyotypes, and use the growth of their neighbours to set a lower bound on the missegregation rate

**Fixed:** 2026-11-20, before any v5 computation.

**User confirmation:** 2026-11-20 around 19:55. Quote: "A frequency threshold of 10 cells is fine. Only a lower-bound correction: only raise, capped at 10-fold, no point estimate. Agreed. Run steps 1 and 2 automatically; if step 2 passes, go on to 3 and 4; if it does not pass, stop and report. Fine. Confirmed"

**Premises (user):**
- Sampling and drift are both relatively small perturbations.
- Changes that even the upper missegregation-rate limit of 1e-2 cannot explain are taken as high-confidence real changes (selection), and are no longer treated as noise.
- Sudden large changes between the two time points correspond to selection; the stable frequent karyotypes, which make up the majority, indicate stable fitness.
- Correct the missegregation rate with the observed counts of the one-step neighbours of frequent karyotypes. The threshold is borrowed from alfak's minobs, but plays a different role.

**Not done:**
- Do not change the KARST package; no commit.
- Do not generate new simulations; do not touch another session's tasks.
- Do not use the sealed held-out landscapes (seed 98000).

## Data

- **Simulations:** the v3 simulations, taking the no-injection version from the v4 screening records (`variant = none`), 1632 datasets in total.

  | Set | True missegregation rate | Datasets |
  |---|---|---|
  | e1, e2, e3 | 5e-5 (×1) | 288, 192, 288 |
  | e4 | 2e-4 (×4), 8e-4 (×16) | 144 each |
  | C2 | 1.5e-4 (×3), 5e-4 (×10), 1.5e-3 (×30), each split into bottleneck 300, 3000 | 96 each |

  - 300 or 1000 cells per time point.
  - Only the first and the last time point are used.
- **Noise-injected versions** (dm20, reg10): report only, not part of the criteria.
- **Halves:** reuse v4's `v4_half()`; the development half chooses the version, the confirmation half decides.
- **Real data:**
  - Salehi: 155 trajectories, first and last passage, ρ = 1 division/day.
  - NUMBAT: 41, ρ = 0.141, only the summary is reported.

## Definitions

**Frequent karyotype:** ≥ 10 cells at at least one time point.

**Significant change:**
- Reuse the envelope test `env_analyse()`. Not explained by any missegregation rate in [5e-5, 1e-2], and |z| > 3.
- Direction = the sign of z.

**Reference karyotypes:** the frequent karyotypes used to set the lower bound.
- Only frequent karyotypes that are not declining count.
- "Declining" uses a two-proportion test that ignores missegregation: z = (p₁ − p₀) / √(p̄(1 − p̄)(1/n₀ + 1/n₁)); z < −1.96 means declining.
- The envelope's decline test is not used. Over long times, the upper-limit missegregation rate expects almost all cells to leave, so that test has no power.

**Neighbour cloud:**
- Karyotypes at L1 distance 1 from some reference karyotype that themselves have < 10 cells at both time points.
- Karyotypes at distance 1 from any non-reference frequent karyotype (the declining ones) do not count, so that inflow from those does not mix in.
- Two versions, chosen on the development half:
  - **a:** remove neighbours called significantly rising (possibly selected, fitter than the parent karyotype);
  - **b:** on top of a, also remove the 2 neighbours with the largest total cell count (sum over both time points).

**Number of inflow routes:**
- The expected number of one-step-neighbour cells produced per division of reference karyotype P is about β·G_P.
- G_P = Σ_i c_i · (1 + 1[c_i ≥ 2]): gaining one copy is always possible; losing one requires c_i ≥ 2.
- No deduction is made for the copy-number cap, so G_P is on the high side and the lower bound is conservative.
- Ḡ is the weighted mean of G_P over the reference karyotypes, with weights x₀/n₀ + x₁/n₁.

**Lower bound:**
- r_t = neighbour-cloud cells ÷ reference-karyotype cells, at time point t.
- Point value: β_pt = (r₁ − r₀) / (ρ · Ḡ · Δt).
- Derivation: when the neighbours are on average no fitter than the reference karyotypes, dr/dt ≤ ρβḠ. ρ is the upper limit of the division rate, with values as in the v3 and v4 screening (1 for simulations and Salehi, 0.141 for NUMBAT).
- Lower confidence bound L: p_t ~ Beta(neighbour cloud + 0.5, reference + 0.5), r_t = p_t/(1 − p_t), 4000 draws, take the 5% quantile of (r₁ − r₀)/(ρḠΔt). The seed is determined by the dataset name.
- No reference karyotype, or a neighbour cloud of 0 at both time points: L = NA, no trigger.

**Trigger and the missegregation rate used:**
- Trigger: L > 5e-5.
- Missegregation rate used β_used = min(5e-4, max(5e-5, L)): only raise, capped at 10-fold.

## Step 1: test the premises (report only, counted on HPC)

1. **Neutral e3:** fraction of frequent karyotypes called significantly changed, split by bottleneck (200, 2000, 20000) and cell count.
   - ≤ 5% counts as "noise is a small perturbation" holding.
   - Also report the fraction of datasets with at least one.
2. **Simulations with selection (e1, e2, e4, C2):** among frequent karyotypes called significantly changed, the fraction whose direction agrees with the true fitness.
   - Agrees: rising with a true net growth rate above the sample mean (weighted by the cells pooled over both time points); or declining with one below the sample mean.
   - Truth: `sm_true_net_rate()` for the e sets, `c2_true_rate()` for C2.
   - Expected ≥ 90%.
3. **Salehi parallel lineages:**
   - Trajectory pairs: same cell line, same treatment at both ends, non-overlapping samples, as in v4.
   - For the same frequent karyotype called significantly changed in both trajectories, the fraction with the same direction. ≥ 75% supports, close to 50% does not support.
   - Also report: the fraction called significant in one trajectory and frequent but not significant in the other.
4. **Report only:** the fraction of reference karyotypes in the simulations whose true division rate is above ρ = 1. Above ρ, the lower bound will be too high.

## Step 2: lower bound on the missegregation rate (counted on HPC, evaluated locally)

**Criteria (simulations without injection):**

| Criterion | Requirement |
|---|---|
| Valid | datasets with L above the true missegregation rate ≤ 5% (overall), and ≤ 10% within each set |
| Useful | among datasets whose true missegregation rate is ≥ 4 times the default (both e4 levels, C2 5e-4 and 1.5e-3), ≥ 20% trigger |

- Each set: e1, e2, e3, e4×4, e4×16, and the 6 settings of C2, 11 in total.
- **Development half:** evaluate both a and b. If both pass, take the one with the higher "useful" fraction; if equal, take b. If neither passes, step 2 fails; stop and report.
- **Confirmation half:** the chosen version must pass the same criteria before going on to steps 3 and 4.
- **Report only:**
  - the noise-injected versions;
  - trigger rates by interval, bottleneck and cell count;
  - trigger counts and the distribution of L for Salehi and NUMBAT (NUMBAT summary only).

## Step 3: fits (HPC, triggered datasets without injection in the confirmation half)

**Fits:**
- T0: β = 5e-5.
- T3: β = β_used.
- Library karst_v4, without `noise_phi`, bit-for-bit identical to 1.1.3. Settings the same as v4's `fit_task.R`.
- Existing v4 T0 fits are used directly; missing ones are run.

**Error ratio:** the mean of RMSE_T3 / RMSE_T0 over distances 0, 1, 2.

**Criteria:**
- (i) Triggered datasets whose true missegregation rate is above the default: median error ratio ≤ 0.95, and the 95% upper limit of a landscape-cluster bootstrap (2000 replicates) < 1.
- (ii) Triggered datasets whose true missegregation rate equals the default: median error ratio ≤ 1.02. With fewer than 10, report only, and use step 2's validity criterion instead.
- Triggered datasets in the development half are also fit, report only.

## Step 4: Salehi (HPC, triggered trajectories)

**Fits and predictions:**
- T0 and T3 fits, settings the same as v4's `salehi_fit.R`.
- Prediction of later passages (≤ 60 days). The missegregation rate used for prediction is fixed at 5e-5 for both, so only the β used in fitting is compared, avoiding the known prediction dispersion effect (`beta/decomp`).
- T3J: T3 fit, predicted with β_used, report only.

**Criteria:**
- Confirmation group (SA906, SA532, SA535, SA1035): per-cell log score gain ≥ 0, and the 95% lower limit of a trajectory-cluster bootstrap > −0.1.
- If fewer than 10 trajectories trigger in the confirmation group, report only, no verdict.
- Development group: report only.

## Outputs

- Code: `dev/realdata/noise/R/anchor.R`, `dev/realdata/noise/v5/`.
- Results: HPC `results/realdata/noise/v5/`, local copy in `dev/results/realdata/noise/v5/`.
- Report: `dev/realdata/noise/report_v5.md`.
