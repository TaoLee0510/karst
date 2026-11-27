# Underestimated out-of-domain karyotypes: results of steps 0 and 1 (stopped at step 1, as planned)

English version of `report.md` (Ver2, 2026-11-25). The content is unchanged;
the Chinese original is kept.

2026-11-25 03:00 EST. The plan and its pre-registered criteria are in
`PLAN.md` (English: `PLAN_Ver2.md`), options D1a–D4b. All code is in
`patch/outside_domain/`; KARST was not modified.

## Conclusions

1. **Step 0 passed:** the deterministic engine can stand in for the ABM. On
   codex's 35 Salehi windows:
   - the median difference between the daughter angle and the mean of 5 ABM
     runs is 0.08° for the ALFA-K rate map and 0.41° for the KARST rate map,
     with a maximum of 2.9°;
   - the median difference in the three group shares is below 0.1
     percentage points.
2. **Step 1's verdict is "not reproduced", so the study stopped as planned,
   and step 2 was not judged.** In the development half, 640 inputs saw the
   near out-of-domain share (≤ 4 steps from a major) actually rise:
   - KARST (native) predicted a fall in only 34.5% of them, below the
     pre-registered half. The aligned version predicted a fall in 48.4%.
   - The true landscape got the direction right in 83.6%, and in 92.5% with
     the true β of the missegregation simulations.
   - The corresponding share on Salehi is 28/30 (93%).
3. **In the simulations KARST does not underestimate the fitness of observed
   out-of-domain karyotypes; it slightly overestimates it.**
   - Median bias of KARST's relative fitness minus the truth: +0.022 at 2
     steps from a major, +0.042 at 3–4 steps, +0.082 at > 4 steps (per day).
   - The true medians are −0.017 to −0.060: these karyotypes are on average
     worse than the majors, and KARST's common level is not below them.
   - Exceptions: main simulations at λ 0.2 have bias −0.13, and λ 1.6 has
     −0.014.
   - **The real weakness is ranking:** the Spearman correlation between the
     observed out-of-domain karyotypes and the truth is about 0 (−0.06 to
     −0.03). KARST gives them all the same level and cannot tell good from
     bad.
4. **A Salehi-like failure appears in only one kind of simulation: highly
   dispersed populations with β too low.** There, the cause is too little
   missegregation flux into new karyotypes, not fitness (see the table
   below).
   - At missegregation rates 2e-4 and 8e-4, with β fixed at 5e-5, KARST
     predicted a fall in 60% and 71% of the inputs where the near
     out-of-domain share rose.
   - In the same groups, the true landscape at the fitted β (5e-5) was often
     wrong too: right in only 63% and 54%.
   - With the true β, the true landscape was right in 98% and 91%.
5. **What this means for the candidate rules:**
   - B1 (raising observed out-of-domain karyotypes to the neutral level)
     goes the wrong way in the simulations: KARST already slightly
     overestimates them.
   - B2 changes what neighbours inherit; there is likewise no
     underestimation for it to correct.
   - As planned, the step-2 results were not analysed. The B arms were
     computed in the same jobs and stored in the result files; I have only
     seen their printed values on the 5 pilot inputs.
6. **The Salehi pattern (KARST predicts a falling out-of-domain share in
   40/44) is stronger than in any simulated group.** It is closest to the
   "highly dispersed + β too low" group, consistent with fig2d_misseg and the
   β grid: on Salehi, raising β lifts KARST's Fig. 2d from 0.22 to 0.74.
   Whether there are other mismatches cannot be told apart here.

## 1. Step 0: engine check

| Rate map | Daughters | Median angle difference | Maximum angle difference | Median difference in change direction between the two forecasts | Median group-share difference |
|---|---|---|---|---|---|
| ALFA-K | 44 | 0.08° | 1.39° | 0.50° | 0.06 percentage points |
| KARST native | 44 | 0.41° | 2.94° | 1.55° | 0.09 percentage points |

The ABM-side angles differ from codex's numbers by at most 2.7e-10 degrees.
Results: `results/step0/engine_check.tsv`.

## 2. Step 1: data

- **Main simulations:** 400 × 3 designs = 1,200 inputs.
- **Missegregation simulations:** 375 β-fixed inputs and 300 β-true inputs.
  36 of them had no KARST fit in the original study either (domain errors at
  high missegregation rates).
- 1,839 inputs are usable. For all of them, the counts used by the fit match
  the rebuilt input.
- **Forecast setting:** start from the parent used for fitting (Salehi cell
  counts) and forecast on list E to the second sampling time (about 13
  days). The observed value is the 1000-cell sample stored by the
  simulation.

## 3. Direction of the near out-of-domain share (both halves, by group)

The denominator of "KARST predicts a fall" and of "true landscape direction
right" is the inputs in that group where the share actually rose.

| Group | Actually rose | KARST native predicts a fall | KARST aligned predicts a fall | True landscape (fitted β) direction right | True landscape (true β) direction right |
|---|---|---|---|---|---|
| Main λ 0.2 | 70/300 | 21% | 21% | 94% | — |
| Main λ 0.4 | 162/300 | 36% | 47% | 74% | — |
| Main λ 0.8 | 240/300 | 31% | 53% | 90% | — |
| Main λ 1.6 | 261/300 | 19% | 40% | 99% | — |
| Missegregation 5e-5 | 63/75 | 25% | 48% | 98% | — |
| Missegregation 1e-4, β fixed | 62/75 | 35% | 42% | 63% | 81% |
| Missegregation 2e-4, β fixed | 57/74 | **60%** | 63% | 63% | 98% |
| Missegregation 4e-4, β fixed | 66/75 | 44% | 47% | 30% | 94% |
| Missegregation 8e-4, β fixed | 35/58 | **71%** | 69% | 54% | 91% |
| Missegregation 1e-4, β true | 62/75 | 35% | 52% | 81% | — |
| Missegregation 2e-4, β true | 57/74 | 54% | 75% | 98% | — |
| Missegregation 4e-4, β true | 66/75 | 47% | 52% | 94% | — |
| Missegregation 8e-4, β true | 35/58 | 51% | 57% | 91% | — |

- The results are the same with the 3-step and 6-step definitions of the
  near out-of-domain group.
- **Magnitude:** forecasts from the true landscape also give only about half
  of the observed rise. For example, at main λ 0.8 the observed median is
  +0.80 percentage points, the true landscape +0.39, KARST +0.26. Possible
  reasons:
  - offspring outside list E are dropped;
  - the parent is a Salehi-sized small sample.

  This part has nothing to do with fitness; it affects KARST and the true
  landscape alike.
- **Fig. 2c-style angle** (medians):
  - main simulations: KARST 28°–33°, true landscape 14°–18°;
  - missegregation 4e-4 / 8e-4, β fixed: KARST 68° / 81°, true landscape
    27° / 25°;
  - the same two levels, β true: KARST 63° / 67°.

## 4. Fitness bias (A1, both halves)

| Object | Distance from majors | Count | Median bias (KARST − truth) | Median truth | Spearman |
|---|---|---|---|---|---|
| Observed out-of-domain karyotypes | 2 steps | 55,421 | +0.022 | −0.017 | −0.03 |
| Observed out-of-domain karyotypes | 3–4 steps | 48,783 | +0.042 | −0.031 | −0.06 |
| Observed out-of-domain karyotypes | > 4 steps | 29,546 | +0.082 | −0.060 | −0.04 |
| Unseen neighbours | 2 steps | 206,822 | +0.023 | −0.020 | 0.35 |
| Unseen neighbours | 3–4 steps | 3,174,354 | +0.055 | −0.050 | 0.21 |
| Unseen neighbours | > 4 steps | 1,747,229 | +0.126 | −0.081 | 0.02 |

Units are net growth rate per day. The truth is zeroed on KARST's own
reference (the count-weighted mean of the majors).

**Median bias of observed out-of-domain karyotypes, by group:**
- main simulations: λ 0.2 −0.13; λ 0.4 +0.03; λ 0.8 0.00; λ 1.6 −0.01;
- missegregation simulations, β fixed: +0.02 to +0.09 (larger at higher
  missegregation rates);
- missegregation simulations, β true: +0.00 to +0.03.

## 5. Options for the next step (the user's decision)

1. **Close the "raise out-of-domain fitness" direction.** The simulations do
   not support underestimation, and B1 points the other way in them. I
   recommend this option.
2. **Flux direction:** raise the missegregation rate used for forecasts on
   highly dispersed data. The simulations show this is the cause of the
   failure. But it belongs to the β line (only externally given values are
   accepted), so it needs the user's decision. The β grid has already shown
   the default and the best value side by side on Salehi.
3. **Ranking direction:** observed out-of-domain karyotypes are hardly ranked
   at all (Spearman ≈ 0), a real weakness. The validation README also lists
   "adding the trend's variation to the common level: untried". This can get
   its own plan.
4. **Analyse the B arms as exploration.** They are already in the result
   files, but this was not a pre-registered judgement and would be
   descriptive only.

## 6. Cost and files

- **CPU:** 128 hours on the HPC (1,875 tasks plus the pilot); step 0 locally
  took about 40 minutes on 3 cores.
- **Code:**
  - `R/`: engine.R (the deterministic engine), variants.R (the candidate
    rules), common.R (shared functions);
  - `scripts/`: step0_engine_check, sim_task, step1_analyse, step2_analyse,
    make_sim_tasks;
  - `hpc/`: job.sh, submit.sh.
- **Results:**
  - locally in `results/`: step0; step1 (status, forecast_rows, a1_rows,
    verdict); sim_out (one rds per input);
  - on the HPC in `/share/lab_crd/taoli/Project/karst_patch/outside_domain/`.
