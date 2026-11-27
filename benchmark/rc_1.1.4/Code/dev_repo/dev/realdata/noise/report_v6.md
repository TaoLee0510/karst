# Pre-fit preprocessing v6: correcting the missegregation rate with the level of the neighbour cloud

**Basis:**
- Plan `PLAN_v6.md`, fixed 2026-11-20 before any v6 computation. The user confirmed at about 20:30: "go ahead with all of your suggestions".
- Decision numbers come from `v6/calib_eval.R`; the output is in `dev/results/realdata/noise/v6/calib_eval.log`.
- Section 4 is a post-hoc description, from `v6_posthoc.R` → `posthoc.log`, not a decision.
- Computation: job 22236501, all 17 tasks COMPLETED. 2579 records in total: v3 simulations 1632, paper simulations 750, Salehi 155, NUMBAT 42, all succeeded.

## Conclusion: step 2 passed, step 3 failed on both test sets; stopped according to the plan

KARST keeps the default missegregation rate 5e-5. Step 4 (fitting) and step 5 (Salehi) were not done.

**Step 2 (development half, 816 datasets): passed**

| | A (neighbour cloud of frequent karyotypes) | B (adjacent cell pairs ÷ identical cell pairs) |
|---|---|---|
| AUC (×1 vs ≥ ×4) | 0.525 | **0.838** |

- B chosen. A cannot be computed on 136 datasets at ≥ ×4 that have no frequent karyotypes; it can be computed on all ×1 datasets.
- Calibration: log₂ S_B = −7.085 + 0.434 · log₂(true multiple).
  - 95% interval of the slope [0.385, 0.483], satisfying d > 0.
  - Residual standard deviation σ = 1.31, R² = 0.27.

**Step 3: failed**

| Test | False raise (×1, β_used ≥ 2 times; required ≤ 10%) | Useful (≥ ×4, β_used ≥ 2 times; required ≥ 50%) | No overshoot (β_used > 2 × true value; required ≤ 5%) |
|---|---|---|---|
| Confirmation half (816) | **36.2%** (of 384) | 76.8% | **19.6%** |
| Paper simulations (750) | **24.7%** (of 150) | 95.3% | **16.4%** |

- Both test sets passed "useful" and failed "false raise" and "no overshoot".
- That is: when the true missegregation rate is high, most datasets are adjusted upwards, but data at the default missegregation rate are also often adjusted upwards.

## 1. Discrimination is very good on the paper simulations, moderate on v3 (report only)

| Test | AUC (×1 vs ≥ ×4) |
|---|---|
| Development half | 0.838 |
| Confirmation half | 0.820 |
| Paper simulations | **0.968** |

Paper simulations by window and true multiple, fraction with β_used ≥ 2 times:

| True multiple | Window 2 (about 11 days) | Window 4 (about 32 days) | Window 8 (about 74 days) |
|---|---|---|---|
| ×1 | 0.32 | 0.30 | 0.12 |
| ×2 | 0.80 | 0.78 | 0.44 |
| ×4 | 0.84 | 0.86 | 0.88 |
| ×8, ×16 | 1.00 | 1.00 | 1.00 |

## 2. Fraction adjusted upwards by set (report only, halves separate)

| Set | True multiple | Development: β_used ≥ 2 times | Confirmation: β_used ≥ 2 times | Confirmation: median β_used |
|---|---|---|---|---|
| e1 / e2 / e3 | ×1 | 0.18 / 0.35 / 0.26 | 0.31 / 0.43 / 0.37 | ×1.2–×1.7 |
| C2 1.5e-4, bottleneck 300 / 3000 | ×3 | 0.63 / 0.67 | 0.40 / 0.65 | ×1.7 / ×2.9 |
| e4 ×4 | ×4 | 0.54 | 0.44 | ×1.5 |
| C2 5e-4, bottleneck 300 / 3000 | ×10 | 1.00 / 1.00 | 0.98 / 1.00 | ×5.5 / ×10 |
| e4 ×16 | ×16 | **0.10** | 0.58 | ×5.3 |
| C2 1.5e-3, bottleneck 300 / 3000 | ×30 | 1.00 / 0.88 | 1.00 / 0.85 | ×9.8 / ×6.6 |

- e4 ×16 differs a lot between the two halves (0.10 vs 0.58). This comes only from the different landscapes; the design of the two halves is the same.
- The difference between 300 and 1000 cells per time point is small: e sets 0.35 vs 0.33.

## 3. Real data (report only)

| Data | Trajectories | Adjusted upwards | ≥ 2 times | Median β_used |
|---|---|---|---|---|
| SA906 | 36 | 1 | 0 | ×1 |
| SA039 | 3 | 0 | 0 | ×1 |
| SA532 | 21 | 10 | 1 | ×1 |
| SA1035 | 12 | 7 | 3 | ×1.1 |
| SA535 | 16 | 9 | 4 | ×1.3 |
| SA609 | 67 | 67 | 67 | ×8.6 |

- NUMBAT (summary): GBM 14 of 20 ≥ 2 times, ovarian cancer 12 of 13, HNSC 4 / 4, melanoma 1 / 5.
- Because step 3 failed, these numbers cannot be taken as estimates of the missegregation rate.

## 4. Post-hoc description: why it failed (not a decision)

**Median log₂ S_B, by simulation source and true multiple:**

| Source | ×1 | ×2 | ×3 | ×4 | ×8 | ×10 | ×16 | ×30 |
|---|---|---|---|---|---|---|---|---|
| e sets | −6.20 | | | −5.78 | | | **−5.97** | |
| C2 | | | −5.68 | | | −5.09 | | −4.91 |
| Paper simulations | −6.38 | −5.51 | | −5.06 | −4.57 | | −4.22 | |

**Slopes fitted separately per source (report only):**

| Source | Slope d | Residual σ | R² |
|---|---|---|---|
| e sets | 0.29 | 1.43 | 0.08 |
| C2 | 0.24 | 0.53 | 0.27 |
| Paper simulations | 0.54 | 0.51 | 0.69 |

**Three observations:**

1. **In the e sets, S_B hardly changes with the missegregation rate.**
   - The median of e4 ×16 (−5.97) is even lower than that of ×4 (−5.78), and the spread at the same multiple is large.
   - Fitted separately, the slope of the e sets is only 0.29, R² 0.08.
2. **Different sources have different baselines.**
   - At the same multiple, the level of S_B differs by 0.7–1 log₂ units, e.g. at ×4, e sets −5.78 and paper simulations −5.06. With the slope of the paper simulations, this corresponds to a 2–4 times difference in missegregation rate.
   - The slope of the pooled calibration, 0.434, is larger than the separate slopes of both sources used in the calibration (development half: e sets 0.29, C2 0.20). This shows that part of it comes from the difference between sources: C2 has both a high missegregation rate and a high S_B. This is the risk mentioned in the plan.
3. **Drift changes the baseline.**
   - In the e sets at ×1, the fraction with β_used ≥ 2 times rises as the bottleneck gets larger: 17% at bottleneck 200, 32% at 2000, 39% at 20000.
   - Strong drift clears out rare neighbours at passage, so S_B becomes lower; when drift is weak, ×1 data also appear to have "many neighbours".

**Conjecture from this (not verified):**
- The level of the neighbour cloud is mainly set jointly by "missegregation rate ÷ neighbour disadvantage" and passage drift.
- In the paper simulations, whose design is uniform, it discriminates the missegregation rate very well (AUC 0.97).
- But the calibrated intercept cannot be moved from one kind of simulation to another, let alone to real data: neither the neighbour disadvantage nor the passage bottleneck of real cells is known.

**Median log₂ S_B on real data:**

| Cell line | SA039 | SA906 | SA532 | SA1035 | SA535 | SA609 |
|---|---|---|---|---|---|---|
| log₂ S_B | −8.30 | −7.37 | −6.84 | −6.14 | −6.07 | −4.86 |

- SA906 and SA039 are below the ×1 median of all simulations. SA609 is comparable to ×10–×16 in the simulations.
- The order agrees with the direction of the β study: SA906 low, the others high. But as shown above, the absolute value of the level cannot be converted into a missegregation rate.

## Next steps (waiting for your decision; I have not done these)

1. **Accept the current state:** the missegregation rate fixed at 5e-5. S_B (or the neighbour-cloud ratio) as a descriptive quantity in the data-quality report, giving only relative position, not converted into a missegregation rate.
2. **If this is to continue, external information that can pin down the baseline is needed.** For example:
   - a control of the same cell line with a known missegregation rate (an experimentally measured missegregation rate, or a parental line known to have a low missegregation rate). Use the S_B of the control to set the intercept, and use only the slope to convert multiples.
   - This needs your judgement on whether such data exist. This round has already shown that setting the intercept from simulations alone does not work.
3. **Not recommended:** changing the statistic or the calibration method again on simulations. The problem in this round is that the intercept changes with the source, not the choice of statistic.

## Files

- Code:
  - `dev/realdata/noise/R/cloud_level.R`;
  - `dev/realdata/noise/v6/stats_task.R`, `calib_eval.R`;
  - task table `v6/tasks/stats.tsv`.
- Post-hoc description: `dev/results/realdata/noise/v6/v6_posthoc.R` → `posthoc.log`, in the results directory, not put in the repository.
- Results:
  - HPC `results/realdata/noise/v6/stats/`;
  - local `dev/results/realdata/noise/v6/`: `stats/`, `calibration.rds`, `calib_eval.log`, `posthoc.log`.
- Job: 22236501 (17 tasks), library karst_v4, no fitting.
