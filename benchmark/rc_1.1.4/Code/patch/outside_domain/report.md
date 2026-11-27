# Underestimate of out-of-domain karyotypes: results of steps 0, 1 (stopped at step 1 as per the plan)

2026-11-25 03:00 EST. The plan and the pre-set criteria are in `PLAN.md` (D1a–D4b). All code is in `patch/outside_domain/`; KARST was not modified.

## Conclusions

1. **Step 0 passed:** the deterministic engine can replace the ABM. On codex's 35 Salehi windows:
   - median difference between the daughter angle and the mean of 5 ABM runs: ALFA-K rate map 0.08°, KARST rate map 0.41°, maximum 2.9°;
   - the median difference of the three group fractions is below 0.1 percentage points.
2. **Step 1 judged "not reproduced"; stopped as per the plan, no verdict for step 2.** In the development half, there are 640 inputs in which the near out-of-domain share (≤4 steps from the major karyotypes) actually rises:
   - KARST (native) predicts a fall in only 34.5% of them, less than the pre-set half; after aligning the level, 48.4%.
   - The true landscape gets the direction right in 83.6%; in the missegregation-rate simulations with the true β, 92.5%.
   - The corresponding fraction on Salehi is 28/30 (93%).
3. **In the simulations KARST does not underestimate the fitness of observed out-of-domain karyotypes; rather, it slightly overestimates it.**
   - Median bias of KARST relative fitness minus truth: 2 steps from the major karyotypes +0.022, 3–4 steps +0.042, >4 steps +0.082 (per day).
   - The medians of the true values are −0.017 to −0.060: these karyotypes are on average worse than the major karyotypes, and the common level KARST gives is not lower than theirs.
   - Exceptions: main simulations λ 0.2 bias −0.13, λ 1.6 bias −0.014.
   - **The real weakness is ranking:** the Spearman between observed out-of-domain karyotypes and the truth is about 0 (−0.06 to −0.03). KARST gives them one common level and cannot tell good from bad.
4. **A Salehi-like failure appears in only one kind of simulation: the population is very dispersed and β is too low.** In this kind of simulation the cause is that the flux of new karyotypes produced by missegregation is too small, not fitness (see the table below).
   - At missegregation rates 2e-4 and 8e-4 with β fixed at 5e-5, among the inputs whose near out-of-domain share rises, KARST predicts a fall in 60% and 71%.
   - In the same groups, the true landscape at the fitted β (5e-5) is also often wrong: the direction is right in only 63% and 54%.
   - With the true β instead, the fraction in which the true landscape gets the direction right rises to 98% and 91%.
5. **Implications for the candidate changes:**
   - B1 (raising observed out-of-domain karyotypes to a neutral level) goes in the wrong direction in the simulations: KARST already slightly overestimates them.
   - B2 changes the inheritance of neighbours; likewise there is no "underestimate" to correct.
   - As per the plan, the step 2 results were not analysed. The B arms were computed in the same batch of jobs and are stored in the result files; I have only seen their printed values on the 5 trial-run inputs.
6. **The Salehi pattern (KARST predicts a falling out-of-domain share in 40/44) is stronger than in any simulation group.** It is most like the "highly dispersed + β too low" group, consistent with the conclusions of fig2d_misseg and the β grid: raising β on Salehi lifts KARST's Fig. 2d from 0.22 to 0.74. Whether there are other mismatches as well cannot be told apart here.

## 1. Step 0: engine check

| Rate map | Daughters | Median angle difference | Max angle difference | Difference in change direction between the two forecasts (median) | Median group-fraction difference |
|---|---|---|---|---|---|
| ALFA-K | 44 | 0.08° | 1.39° | 0.50° | 0.06 percentage points |
| KARST native | 44 | 0.41° | 2.94° | 1.55° | 0.09 percentage points |

The angles on the ABM side differ from codex's numbers by at most 2.7e-10 degrees. Results: `results/step0/engine_check.tsv`.

## 2. Step 1: data

- **Main simulations:** 400 × 3 designs = 1,200 inputs.
- **Missegregation-rate simulations:** 375 inputs with β fixed, 300 with β at the true value; 36 of them already had no KARST fit in the original study (domain errors at high missegregation rates).
- 1,839 usable inputs in total; the counts used by the fits and the reconstructed inputs all agree.
- **Forecast setting:** starting from the parent used for the fit (Salehi cell counts), forecast on the E table to the second sampling (about 13 days). The actual values use the 1000-cell samples stored by the simulation.

## 3. Direction of the near out-of-domain share (both halves pooled, by group)

The denominator of both "KARST predicts a fall" and "true landscape direction right" is the inputs in that group in which the share actually rises.

| Group | Actually rising | KARST native predicts fall | KARST aligned predicts fall | True landscape (fitted β) direction right | True landscape (true β) direction right |
|---|---|---|---|---|---|
| Main simulations λ 0.2 | 70/300 | 21% | 21% | 94% | — |
| Main simulations λ 0.4 | 162/300 | 36% | 47% | 74% | — |
| Main simulations λ 0.8 | 240/300 | 31% | 53% | 90% | — |
| Main simulations λ 1.6 | 261/300 | 19% | 40% | 99% | — |
| Missegregation rate 5e-5 | 63/75 | 25% | 48% | 98% | — |
| Missegregation rate 1e-4, β fixed | 62/75 | 35% | 42% | 63% | 81% |
| Missegregation rate 2e-4, β fixed | 57/74 | **60%** | 63% | 63% | 98% |
| Missegregation rate 4e-4, β fixed | 66/75 | 44% | 47% | 30% | 94% |
| Missegregation rate 8e-4, β fixed | 35/58 | **71%** | 69% | 54% | 91% |
| Missegregation rate 1e-4, β at true value | 62/75 | 35% | 52% | 81% | — |
| Missegregation rate 2e-4, β at true value | 57/74 | 54% | 75% | 98% | — |
| Missegregation rate 4e-4, β at true value | 66/75 | 47% | 52% | 94% | — |
| Missegregation rate 8e-4, β at true value | 35/58 | 51% | 57% | 91% | — |

- The results are the same under both near out-of-domain definitions, 3 steps and 6 steps.
- **Magnitude:** the forecast of the true landscape also gives only about half of the actual rise. For example, main simulations λ 0.8: actual median +0.80 percentage points, true landscape +0.39, KARST +0.26. Possible causes:
  - offspring outside the E table are discarded;
  - the parent is a small sample of Salehi size.

  This part has nothing to do with fitness; both KARST and the true landscape are affected.
- **Fig. 2c-style angle** (median):
  - main simulations: KARST 28°–33°, true landscape 14°–18°;
  - missegregation rate 4e-4 / 8e-4, β fixed: KARST 68° / 81°, true landscape 27° / 25°;
  - the same two levels, β at true value: KARST 63° / 67°.

## 4. Fitness bias (A1, both halves pooled)

| Object | Distance from major karyotypes | Count | Median bias (KARST − truth) | Median truth | Spearman |
|---|---|---|---|---|---|
| Observed out-of-domain karyotypes | 2 steps | 55,421 | +0.022 | −0.017 | −0.03 |
| Observed out-of-domain karyotypes | 3–4 steps | 48,783 | +0.042 | −0.031 | −0.06 |
| Observed out-of-domain karyotypes | >4 steps | 29,546 | +0.082 | −0.060 | −0.04 |
| Unseen neighbours | 2 steps | 206,822 | +0.023 | −0.020 | 0.35 |
| Unseen neighbours | 3–4 steps | 3,174,354 | +0.055 | −0.050 | 0.21 |
| Unseen neighbours | >4 steps | 1,747,229 | +0.126 | −0.081 | 0.02 |

The unit is the net growth rate per day; the truth takes the same reference as KARST as its zero point (the count-weighted mean of the major karyotypes).

**Median bias of observed out-of-domain karyotypes by group:**
- Main simulations λ 0.2: −0.13; λ 0.4: +0.03; λ 0.8: 0.00; λ 1.6: −0.01.
- Missegregation-rate simulations, β fixed: +0.02 to +0.09 (larger at higher missegregation rates).
- Missegregation-rate simulations, β at true value: +0.00 to +0.03.

## 5. Options for the next step (waiting on your decision)

1. **Close the "raise out-of-domain fitness" direction.** The simulations do not support an underestimate, and B1 goes the opposite way in the simulations. I recommend this option.
2. **Flux direction:** raise the missegregation rate used in forecasting on highly dispersed data. The simulations show that this is exactly the cause of the failure, but it belongs to the β line (only externally given values are accepted), so you need to decide. The β grid has already given results on Salehi with the default and the best value side by side.
3. **Ranking direction:** there is almost no ranking among observed out-of-domain karyotypes (Spearman ≈ 0); this is a real weakness. The validation README also records "adding trend variation on top of the common level, not tried". A separate plan can be made.
4. **Analyse the B arms as exploration.** They are already in the result files, but this is not a pre-set verdict, only description.

## 6. Cost and files

- **CPU:** HPC 128 hours (1,875 tasks plus trial runs); locally, step 0 took about 40 minutes on 3 cores.
- **Code:** `R/` (engine.R deterministic engine, variants.R candidate rules, common.R shared functions), `scripts/` (step0_engine_check, sim_task, step1_analyse, step2_analyse, make_sim_tasks), `hpc/` (job.sh, submit.sh).
- **Results:**
  - local `results/`: step0, step1 (status, forecast_rows, a1_rows, verdict), sim_out (one rds per input);
  - on HPC in `/share/lab_crd/taoli/Project/karst_patch/outside_domain/`.
