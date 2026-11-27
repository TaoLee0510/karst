# Decomposing E4: the β used in fitting and the missegregation rate used in prediction

On 2026-11-20 (after 05:09, before 07:26) the user replied "try 2", i.e. item 2 of Section 8 of `../report.md`. This file was written at 2026-11-20 07:26 EST, before any new fits and predictions. It is not changed afterwards; anything that needs correcting goes into the "Corrections" section at the end, with the time and the reason.

## Question

E4 passed on Salehi per the plan: choosing β with the later passages of the other trajectories of the same cell line gives a higher log score than 5e-5 (development +0.18, confirmation +0.44 per cell). But two things change with β at the same time:
- **The β used in fitting (denoted L):** changes the fitness values of the landscape;
- **The missegregation rate used in prediction (denoted F):** changes how dispersed the predicted distribution is.

This study separates the two, to answer: does the gain come from the fitting β itself, or only from a more dispersed prediction? It also extends the grid to ×256, to see whether the choice still hits the upper bound.

## Rules

- Code in `dev/realdata/beta/decomp/`, not committed. The package is not changed; cd342db is used.
- HPC: xxlarge, one fit or one prediction per task; the add-hours monitor uses the small QOS, 3 days.
- After the runs finish and `report.md` is written, stop and wait for the user's decision.

## Method

**Grid:** multiples k ∈ {1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128, 256}, 11 levels in total, ×1 = 5e-5.

**New fits:**
- Subjects: the 155 Salehi trajectories with later passages.
- Content: add the two levels ×128 and ×256, with the same settings as `beta_salehi.R`, writing out alfak-format results at minobs 10.
- The E1–E3 criteria are recorded as usual, report only.

**Predictions:** `validation/sparse_joint/salehi_horizon/forecast_score.R` run unchanged, 31 (L, F) combinations per trajectory:
- **Diagonal:** L = F = k, 11 combinations, i.e. the original E4 design;
- **Fit only:** L = k, F = ×1, another 10;
- **Prediction only:** L = ×1, F = k, another 10.

All combinations are recomputed; the common support is the union of the landscapes of the 11 fit levels. Scoring as in E4: per-cell log score `ll_fc_s_e2` on the common support, later passages ≤ 60 days.

## Three leave-one-out estimates

All are leave-one-out within the same cell line: the choice for each trajectory is determined by the mean score of the other trajectories of the same cell line.

| Estimate | Line on which it chooses | What it chooses |
|---|---|---|
| J | diagonal L = F | k (= E4 on the extended grid) |
| L | F fixed at ×1 | the fitting β |
| F | L fixed at ×1 | the prediction missegregation rate |

- **Gain:** the mean, over all later passages in the group, of the difference in `ll_fc_s_e2` between the chosen combination and (×1, ×1), the same as the original S1.
- **Interval:** 95% interval from a cluster bootstrap by trajectory, 2000 replicates, with the choices held fixed.
- **Groups:** development group SA609, SA039; confirmation group SA906, SA532, SA535, SA1035, as in `../PLAN.md`.

## Verdicts fixed in advance (based on the confirmation group; the development group only checks the direction)

- **A The fitting β itself contributes:** the gain of L > 0 with an interval excluding 0, and the development group's gain of L ≥ 0.
- **B Mainly a more dispersed prediction:** A does not hold, and the gain of F ≥ 75% of the gain of J.
- **C Both:** A holds, and the gain of F > 0 with an interval excluding 0.
- **D Uncertain:** none of the above holds.

## Report only

- Each cell line's curves on the three lines, i.e. the gain of each k relative to (×1, ×1), with mean and median;
- Interaction term: the gain at J's choice − the gain of L at the same k − the gain of F at the same k;
- Whether J's choice on the extended grid is still at the upper bound ×256;
- The E1–E3 criterion values for ×128 and ×256.

## Stopping rules

- A build or script fails for an unknown reason;
- > 5% of the tasks of an array fail, and it is neither a resource problem nor a model-domain error of the package itself;
- A situation the rules do not cover.

## Corrections
- **2026-11-20 07:53 EST (before any prediction was submitted): many ×256 alfak results are missing; the main analysis switches to ×1/4–×128.**
  - **Situation:** of the ×256 fits themselves, 143 succeeded and 12 are model-domain errors, the same as at the other levels. But when the fits are converted to alfak format (`as_alfak`, `correct_efflux = TRUE`), 64 trajectories report "The level of the efflux correction (alfak's constant k) could not be solved for", and 6 report "no frequent karyotypes detected for minobs = 10".
  - **Impact:** of the 137 trajectories with a ×1 result, only 73 have a ×256 landscape. By cell line: SA532 1/21, SA609 28/49, SA906 13/36, the rest complete. ×128 is 137/137.
  - **Nature:** this is a limitation of the package at extremely high β (the missegregation efflux is too large, and alfak's constant cannot be solved for), not a script error.
  - **Handling:**
    - The main analysis (the three estimates, gains, intervals, the A–D verdicts) uses the 10 levels ×1/4–×128, 28 combinations per trajectory, with all 137 trajectories included.
    - ×256 is report only: prediction missegregation rate ×256 (fit fixed at ×1) for all trajectories; fit ×256 only for the 73 with results.
    - The prediction task table is unchanged, still 31 combinations; combinations whose ×256 fit is missing output "no fit".
