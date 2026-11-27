# Splitting E4: results

**Basis:**
- The plan and the decisions fixed in advance are in `PLAN.md`, finalised 2026-11-20 07:26 EST.
- At 07:53 there is one correction: the ×256 alfak results are largely missing, so the main analysis uses ×1/4–×128 instead. This correction was written before any prediction was submitted.
- The numbers come from `decomp_eval.R`; the results are in `results/realdata/beta/decomp/evaluation.rds`.

**Notation:**
- L is the β used for the fit and F is the missegregation rate used for the prediction; both are written as multiples of 5e-5.
- Gain = per-cell log score on the later passages (`ll_fc_s_e2`, common support, ≤ 60 days) minus the value at (L, F) = (×1, ×1); a positive number means better.

## Conclusion: the verdict is "B: mainly the prediction is more dispersed"; the fit β itself contributes nothing

| Estimate | How it is chosen | Development-group gain [95% interval] | Confirmation-group gain [95% interval] |
|---|---|---|---|
| J | L = F chosen together (E4 on the extended grid) | +0.53 [0.38, 0.66] | +0.31 [0.18, 0.46] |
| L | only L chosen, F fixed at ×1 | +0.11 [0.05, 0.19] | +0.03 [−0.04, 0.09] |
| F | only F chosen, L fixed at ×1 | +1.03 [0.87, 1.18] | +0.45 [0.35, 0.55] |

- Development group SA609, SA039: 51 trajectories, 136 later passages; confirmation group SA906, SA532, SA535, SA1035: 85 trajectories, 165 later passages.
- Each estimate is chosen by leave-one-out within the same cell line. Intervals by cluster bootstrap over trajectories, 2000 times.

**Verdict by the rules fixed in advance:**
- **A does not hold:** on the confirmation group the interval of L contains 0.
- **B holds:** on the confirmation group the gain of F, 0.45, is ≥ J's 0.31 × 75%. It also holds when compared with the 0.44 of the previous E4 round.
- The conclusion is "B: mainly the prediction is more dispersed".

**What this shows:**
1. **Changing only the missegregation rate used for the prediction gives all of the gain, and more.** On the landscapes fitted with β = 5e-5, predicting with a missegregation rate of ×16–×128 makes every cell line better. SA906 is no exception: J originally chose ×1/2–×1 on SA906, but with only F changed it gains +0.18 at ×16.
2. **Changing only the fit β brings no benefit.** With F fixed at ×1, L almost always picks the lowest grid level, ×1/4, and the gain is close to 0. The higher L, the worse: at ×64, every cell line except SA039 is worse by 0.5–1.2.
3. **The E4 gain is the result of two effects that partly cancel.** At the level chosen by J:

   | Group | J (L = F) | Only L changed | Only F changed | Interaction term |
   |---|---|---|---|---|
   | Development | +0.53 | −1.26 | +1.03 | +0.76 |
   | Confirmation | +0.31 | −0.85 | +0.39 | +0.77 |

   The landscapes fitted with a high β are themselves worse; only when paired with an equally high prediction missegregation rate do the two partly cancel. A possible reason: the higher β is in the fit, the larger the missegregation outflow that has to be compensated, so the fitness is estimated higher, and this set of fitness values is self-consistent in prediction only at an equally high missegregation rate.
4. **So the "high β" chosen by E4 is not evidence that these cell lines have a higher β.** What it shows is that at 5e-5 the predicted distribution is too concentrated for the later passages, and newly arising karyotypes get too low a probability. The prediction missegregation rate acts here as a "dispersion".
5. **Consistent with O in phase C–E.** O found that 86–95% of the Salehi Pearson excess comes from rare karyotypes that had an expectation < 1 but appeared anyway. Both point to the same phenomenon: the model gives newly arising karyotypes too low a probability.

**The recommendation on β is unchanged, and now better founded:** the β used for fitting in the package stays at 5e-5.

## 1. Curves per cell line

Each cell is the gain relative to (×1, ×1): the mean, with the median in parentheses. The J curve is the same as in section 6 of the previous round's report, with ×128 added.

**Only L changed (F = ×1):**

| Cell line | Later passages | ×1/4 | ×1/2 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 | ×128 |
|---|---|---|---|---|---|---|---|---|---|---|
| SA039 | 2 | −0.01 | −0.02 | +0.02 | +0.02 | +0.02 | +0.01 | −0.01 | −0.04 | −0.04 |
| SA609 | 134 | +0.11 (+0.02) | +0.07 (+0.02) | −0.05 (−0.01) | −0.13 (−0.04) | −0.28 (−0.12) | −0.52 (−0.24) | −0.85 (−0.50) | −1.05 (−0.86) | −1.27 (−1.04) |
| SA1035 | 20 | +0.07 (+0.03) | +0.11 (+0.02) | −0.03 (−0.02) | −0.18 (−0.07) | −0.35 (−0.14) | −0.58 (−0.30) | −0.88 (−0.39) | −0.54 (−0.46) | −1.08 (−0.62) |
| SA532 | 46 | +0.06 (−0.01) | −0.11 (−0.02) | +0.04 (−0.02) | −0.05 (−0.07) | −0.25 (−0.21) | −0.40 (−0.35) | −0.53 (−0.47) | −0.71 (−0.62) | −1.28 (−0.83) |
| SA535 | 49 | +0.13 (+0.09) | +0.01 (+0.01) | −0.00 (−0.05) | −0.18 (−0.10) | −0.34 (−0.25) | −0.61 (−0.52) | −0.83 (−1.16) | −1.18 (−1.44) | −2.47 (−2.13) |
| SA906 | 50 | +0.08 (+0.06) | +0.06 (+0.04) | −0.08 (−0.05) | −0.26 (−0.15) | −0.40 (−0.29) | −0.56 (−0.43) | −0.70 (−0.55) | −1.22 (−0.98) | −2.86 (−2.20) |

**Only F changed (L = ×1):**

| Cell line | Later passages | ×1/4 | ×1/2 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 | ×128 |
|---|---|---|---|---|---|---|---|---|---|---|
| SA039 | 2 | −0.01 | −0.01 | +0.01 | +0.01 | +0.00 | −0.01 | −0.06 | −0.20 | −0.47 |
| SA609 | 134 | −0.34 (−0.29) | −0.18 (−0.16) | +0.21 (+0.19) | +0.40 (+0.38) | +0.61 (+0.56) | +0.79 (+0.73) | +0.94 (+0.84) | +1.03 (+0.92) | +1.04 (+0.98) |
| SA1035 | 20 | −0.14 (−0.14) | −0.07 (−0.06) | +0.08 (+0.08) | +0.15 (+0.14) | +0.23 (+0.19) | +0.29 (+0.24) | +0.30 (+0.26) | +0.23 (+0.27) | +0.05 (+0.19) |
| SA532 | 46 | −0.25 (−0.22) | −0.14 (−0.13) | +0.15 (+0.14) | +0.30 (+0.27) | +0.43 (+0.37) | +0.54 (+0.38) | +0.59 (+0.40) | +0.58 (+0.33) | +0.52 (+0.20) |
| SA535 | 49 | −0.36 (−0.35) | −0.18 (−0.18) | +0.17 (+0.17) | +0.35 (+0.37) | +0.50 (+0.52) | +0.64 (+0.59) | +0.70 (+0.60) | +0.67 (+0.49) | +0.52 (+0.28) |
| SA906 | 50 | −0.10 (−0.07) | −0.05 (−0.04) | +0.06 (+0.04) | +0.11 (+0.09) | +0.16 (+0.13) | +0.18 (+0.17) | +0.16 (+0.18) | +0.04 (+0.10) | −0.24 (−0.07) |

- **The F curve:**
  - on the 5 cell lines other than SA039, it rises monotonically from ×1/4 to ×16;
  - SA906, SA1035, SA532 and SA535 peak at ×16–×32 and decline afterwards;
  - SA609 is still rising at ×64–×128.
- **If all cell lines used a fixed F = ×16:** the gains are SA609 +0.79, SA535 +0.64, SA532 +0.54, SA1035 +0.29, SA906 +0.18, SA039 −0.01. Only SA039 does not improve, and it has only 2 later passages. This is a post hoc observation, with no criterion fixed in advance.
- **The L curve:**
  - except for SA039, L from ×4 upward is worse than ×1 on every cell line, and the higher, the worse;
  - the levels between ×1/4 and ×2 differ little, with means all within ±0.13.

**What each of the three estimates chose (leave-one-out; number of trajectories in parentheses):**

| Cell line | J | L | F |
|---|---|---|---|
| SA039 | ×2 (1), ×8 (1) | ×2 (1), ×8 (1) | ×2 (1), ×8 (1) |
| SA609 | ×128 (49) | ×1/4 (49) | ×128 (49) |
| SA1035 | ×64 (11), ×128 (1) | ×1/4 (1), ×1/2 (11) | ×16 (1), ×32 (11) |
| SA532 | ×64 (21) | ×1/4 (18), ×2 (3) | ×32 (19), ×64 (2) |
| SA535 | ×64 (15), ×128 (1) | ×1/4 (15), ×1/2 (1) | ×32 (16) |
| SA906 | ×1/2 (24), ×1 (12) | ×1/4 (36) | ×16 (36) |

- **Chosen at the new upper bound ×128:**
  - J has 51 trajectories: all 49 of SA609, and 1 each of SA1035 and SA535;
  - F has 49, all of them SA609;
  - L has none; instead, most of its choices sit at the lower bound ×1/4.
- **J's gain on the confirmation group is lower than in the previous E4 round** (0.31 vs 0.44). The reason is that after ×128 was added to the grid, the leave-one-out choice of two trajectories changed from ×64 to ×128 (SA1035T_X5_l_2 and SA535_CISPLATIN_CombinedT_X6_l_2), and they are much worse at ×128 than at ×64. This is noise of the leave-one-out choice and does not affect the verdict.

## 2. ×256 (report only)

Only trajectories that have ×256 results are listed.

| Combination | SA609 | SA535 | SA532 | SA1035 | SA906 | SA039 |
|---|---|---|---|---|---|---|
| L = F = ×256 | +0.02 (28 trajectories) | −3.47 (16) | −2.40 (1) | −0.41 (12) | −1.75 (13) | −1.35 (2) |
| L = ×256, F = ×1 | −3.98 (28) | −5.80 (16) | −5.35 (1) | −2.65 (12) | −5.06 (13) | −0.43 (2) |
| L = ×1, F = ×256 | +0.99 (49) | +0.24 (16) | +0.41 (21) | −0.22 (12) | −0.70 (36) | −1.10 (2) |

- By ×256 the diagonal is already clearly worse. With only F changed, SA609 is flat over ×128–×256, and the other cell lines decline at ×256.
- So the optimum of F lies between ×16 and ×128 and varies by cell line; the upper bound of the grid is already sufficient.

## 3. Criteria E1–E3 on the 10-level grid (report only)

- **Pooled choice over these 137 trajectories:**
  - E1: SA609 ×16, SA532 ×8, SA535 ×4, SA1035 ×4, SA906 ×1/2;
  - E2: still leans towards ×1/4–×1.
- In the individual versions, only 1 trajectory of SA609 picks ×128 under E1 and E3.
- None of these criteria is reliable on simulations (see the previous round's report); they are only recorded here.

## 4. Problems during the run and deviations from the plan

- **The first supplementary fit did nothing:**
  - 22154615 passed the extra multiples through an environment variable, but the container runs with `--cleanenv`, so R could not read it and all tasks skipped straight away. No files were written; it was cancelled.
  - After switching to a script argument (22154777), `submit.sh` passed the arguments through `sbatch --export`, where the comma was treated as a separator, so R received only "128".
  - Switched to colon separation and resubmitted as 22154979, which filled in ×256. The ×128 results come from 22154777, with the same settings.
- **Conversion of ×256 to alfak format largely failed (correction 07:53):**
  - The ×256 fits themselves: 143 succeeded and 12 had model-domain errors, the same as at the other levels.
  - In the conversion to alfak format, 64 trajectories reported that the efflux correction constant could not be solved, and 6 reported no frequent karyotypes. Only 73/137 trajectories have a ×256 landscape, and SA532 only 1/21.
  - This is a limitation of the package at very high β. The main analysis uses ×1/4–×128 instead, with all 137 trajectories and 301 later passages included.
- **Prediction:** all 4805 tasks completed: 4119 succeeded and 686 are "no fit" (all combinations of the 18 trajectories with model-domain errors, plus the combinations with a missing ×256 fit). The 28 combinations are complete on the 301 later passages, and the common support is the same for every combination.
- **The test output** `decomp/test/` is a single-task test of the driver script (22154773) and is excluded from the evaluation.

## 5. What could be done next (for you to decide)

1. **Keep the fit β at 5e-5 (recommended).** The package is not changed. There is now evidence on two sides: β cannot be determined on simulations, and on Salehi raising the fit β by itself brings no benefit either.
2. **If KARST/alfak landscapes are used to predict later passages, the prediction missegregation rate can be raised separately.**
   - This concerns only the prediction step (`p` of `alfakR::predict_evo`); it is unrelated to the fit, and also unrelated to the comparison between KARST and alfak.
   - To become a recommendation, it needs a test fixed in advance, for example a test of a fixed ×16 on the ABM or on another batch of data.
3. **(Done, see `../disp/report.md`: mainly, new karyotypes arise faster and often across several steps, and most of the far ones lie outside the landscape.) Find out why the prediction is too concentrated.** Candidate reasons:
   - bottleneck sampling at passage;
   - karyotypes outside the landscape;
   - multi-chromosome changes beyond the one-step neighbours.
   - This is the same phenomenon as O's excess. The drift term tried in the alfak2 era (state_space) could not be separated from selection; this time it can be tested directly against the dispersion of the later passages.
