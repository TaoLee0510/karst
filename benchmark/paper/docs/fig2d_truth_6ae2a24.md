# Fig. 2d on simulations with known truth: round 6ae2a24

2026-11-23. Numbers come from the paper's 400 ABM simulations, the fits of round 6ae2a24 (step 1), and the results refitted at Salehi cell numbers (step 2). KARST uses this round's library (1.1.3), alfak uses dev2@f93b2d2.

## Conclusion

1. **On simulations with known truth, the ranking in Fig. 2d follows landscape accuracy: true landscape > KARST > alfak.** This is the order at every cell number and every wavelength, both for raw values and with noise held fixed.
   - Step 2 (both fitting and observation use Salehi-like cell numbers): true landscape 0.89, KARST 0.74, alfak 0.58. KARST minus alfak is +0.16 (0.12 to 0.21), +0.20 (0.15 to 0.25) with noise held fixed.
   - Landscape accuracy of the same batch of fits: median r KARST 0.59, alfak 0.37; KARST is higher in 91% of the fits.
2. **The fewer the cells, the lower the ceiling for all methods, but the order does not change.** The true landscape drops from 0.96 at 1000 cells to 0.82 at 100 cells; ranking by noise alone rises from 0.67 to 0.77.
3. **The alfak lead seen on Salehi is not reproduced in the simulations.** At Salehi cell numbers, ranking by noise alone reaches 0.72, comparable to Salehi's own 0.73–0.77; but the true landscape still has 0.89, KARST 0.74. So "few cells per passage, large sampling noise" alone does not make alfak overtake KARST on Fig. 2d.
4. **alfak's forecast step size does not follow noise more than KARST's does.** Spearman of forecast change with noise: true landscape 0.52, KARST 0.49, alfak 0.42. The mechanism conjectured on Salehi does not appear in the simulations.
5. **Therefore the claim "alfak's advantage on Fig. 2d comes from sampling noise" is not supported by this test.** What still holds are two facts on Salehi:
   - with two time points and no selection, the two methods rank equally badly (0.29 vs 0.32);
   - with noise held fixed, there is no reliable difference.
   But on Salehi both methods reach only about 0.3, far below the 0.58–0.74 in the simulations at the same noise. The gap more likely comes from how far the fitted landscape is from the real situation in Salehi, rather than from sampling noise. This point is an inference and has not been checked.
6. **Impact on report Ver10:** it describes Fig. 2d as "follows sampling noise" and "cannot tell a better landscape from a worse one on these data". This wording needs to change; see section 6. The report has not been changed before your confirmation.

## 1. Question

- The paper's Fig. 2d makes a rank comparison (Spearman) on the Salehi passages: on one side "how far the forecast moved from the parent passage", on the other "how far the daughter passage actually moved from the parent", both measured with the Wasserstein distance.
- On this measure alfak is higher than KARST: under the paper's rule, two time points 0.45 vs 0.22, all time points 0.69 vs 0.28.
- The earlier analysis (`scripts/fig2d_noise.R`) found that ranking by noise alone (splitting the cells of each passage randomly into two halves and taking the distance between the halves) already reaches 0.73–0.77, higher than any forecast; with noise held fixed, there is no reliable difference between the two.
- This test redoes Fig. 2d on simulations with known truth, at Salehi cell numbers, to see how the true landscape, KARST and alfak each rank.

Three questions to answer:
- **Q1:** how low a ceiling does sampling noise set for Fig. 2d?
- **Q2:** can the phenomenon on Salehi be reproduced?
- **Q3:** does each method's forecast step size follow noise?

## 2. Method

**Simulations:** the paper's ABM, 400 in total, 100 each at λ = 0.2, 0.4, 0.8, 1.6, missegregation rate 5e-5. The true landscape (GRF) is known. What the paper stored are already samples of 1000 cells per time point, not the complete population, so 1000 cells is the low-noise end, not noise-free.

**Parent and daughter:**
- Parent = the last time point used in the fit, i.e. the last sampling before day 120.
- Daughters = the three samplings after it, with median distances from the parent of 7.1, 14.2 and 21.2 days.
- The main results use the sampling in each simulation closest to 15 days (the paper's evaluation day on Salehi); the next sampling serves as a second set.

**Three forecasts, present for every row:**
- KARST, alfak: their own landscapes;
- true landscape: the same karyotype table (alfak's table at that minobs; in version A KARST's table is identical to it, checked in all 400 simulations), with the fitness replaced by the truth.
- All three use the same forecast program (`R/forecast.R`, alfakR's ABM forecast, 5 replicates, same seeds).
- Forecast change = distance from the parent sample to the forecast; observed change = distance from the parent sample to the daughter sample. Both sides count only the karyotypes in the karyotype table and are renormalised, as in the paper.

**Step 1 (no refitting):**
- Uses the existing fits of this round (input 1000 cells).
- Noise is added only to the observation axis: n = 100, 200, 400 cells are drawn from both the parent and the daughter, 20 draws for each n, with Spearman averaged over the 20 draws; the stored 1000 cells are added as well.

**Step 2 (refitting):**
- The cell number at each time point is drawn from the cell numbers of the two-time-point Salehi inputs: median 389, quartiles 243–529. Seeds are fixed per simulation (`abm_subsample_x` in `R/abm.R`).
- Both methods are refitted on the subsampled inputs: three two-time-point designs w2, w4, w8, fitted at every minobs. Forecasts start from the subsampled parent, as on Salehi.
- Two observation axes: the subsampled samples (i.e. the Salehi situation), and the original 1000 cells.

**Statistics:**
- Raw Spearman, and partial Spearman with noise held fixed. Noise is measured the same way as in `fig2d_noise.R`: the parent and the daughter are each split randomly into two halves, and the mean distance over 20 splits is taken.
- "Ranking by noise alone": Spearman of observed change with noise.
- "Forecast step size follows noise": Spearman of forecast change with noise.
- Intervals: resampling of simulations within each wavelength, 1000 bootstrap replicates.
- Each row = simulation × design × minobs. Only rows where all three forecasts exist are used.
- Two kinds of selection: no selection (every minobs counted); the paper's rule (for each simulation and each design, take the minobs with the highest cross-validation score; the true landscape uses the same rows).

## 3. Checks

- **Step 0:** 4 simulations × two methods × 5 replicates, forecasting only 3 sampling times, are identical to the first rows of the stored forecasts (10 sampling times), maximum difference 0 (`engine_check.tsv`).
- **Subsampling is reproducible:** the fit tasks and the per-row computation each subsample anew, with the same function and seeds. Spot check of 7 simulations × 2 designs × two methods, 28 forecasts in total: the starting point of the method forecasts (the subsampled parent) is identical to the starting point of the true-landscape forecasts, maximum difference 0.
- **Karyotype table:** in every row the karyotype tables of KARST and alfak are identical (`same_list` all true).
- **Failures:**
  - In step 2, 144 / 3,600 alfak fits failed, with the error "Regression matrix for fixed part of model is colinear"; KARST had no failures.
  - Step 1 inherits the failures of this round: 79 / 3,600 for alfak.
  - These rows do not enter the comparison: 3,521 rows remain for step 1, 3,456 rows for step 2.
- **Tasks:** in step 1, 2 tasks on node hpctpa3pc0061 could not start the container ("unknown userid") and wrote nothing; they completed after resubmission.

## 4. Step 1: noise added only to the observation axis

All wavelengths together, no selection, daughter at the sampling closest to 15 days.

| Cells used for observation | True landscape | KARST | alfak | Noise only |
|---|---|---|---|---|
| 1000 | 0.96 | 0.88 | 0.63 | 0.67 |
| 400 | 0.92 | 0.86 | 0.62 | 0.69 |
| 200 | 0.89 | 0.82 | 0.61 | 0.73 |
| 100 | 0.82 | 0.77 | 0.58 | 0.77 |

**With noise held fixed:**

| Cells | True landscape | KARST | alfak |
|---|---|---|---|
| 1000 | 0.95 | 0.83 | 0.47 |
| 400 | 0.90 | 0.79 | 0.45 |
| 200 | 0.85 | 0.74 | 0.43 |
| 100 | 0.73 | 0.64 | 0.37 |

- **KARST minus alfak:** 1000 cells +0.25 (0.22 to 0.29), 100 cells +0.18 (0.15 to 0.22). Above 0 at every cell number.
- **True landscape minus KARST:** only +0.05 to +0.08.
- **Split by wavelength:** same order. The largest difference is at λ = 0.8: at 400 cells, true landscape 0.87, KARST 0.73, alfak 0.37.
- **Answer to Q1:** the fewer the cells, the lower the ceiling. At 100 cells, ranking by noise alone (0.77) is already close to the true landscape (0.82).
- Figure: `analysis/step1/figures/fig2d_truth_cells.pdf`.

## 5. Step 2: fitting also uses only Salehi-like cell numbers

**All wavelengths together, no selection, daughter at the sampling closest to 15 days:**

| | Observation: Salehi-sized samples | With noise held fixed | Observation: 1000 cells | With noise held fixed |
|---|---|---|---|---|
| True landscape | 0.89 | 0.87 | 0.95 | 0.93 |
| KARST | 0.74 | 0.64 | 0.80 | 0.70 |
| alfak | 0.58 | 0.44 | 0.61 | 0.46 |
| Noise only | 0.72 | – | 0.67 | – |

- **KARST minus alfak** (Salehi-sized samples): raw +0.16 (0.12 to 0.21), +0.20 (0.15 to 0.25) with noise held fixed.
- **True landscape minus KARST:** +0.15 (0.12 to 0.18). The KARST fit used only a few hundred cells and is further from the true landscape than in step 1.

**Split by wavelength** (Salehi-sized samples, raw values; in brackets, with noise held fixed):

| Wavelength | True landscape | KARST | alfak | Noise only |
|---|---|---|---|---|
| 0.2 | 0.96 (0.89) | 0.91 (0.78) | 0.57 (0.44) | 0.79 |
| 0.4 | 0.91 (0.86) | 0.72 (0.57) | 0.59 (0.39) | 0.70 |
| 0.8 | 0.79 (0.76) | 0.47 (0.39) | 0.30 (0.21) | 0.60 |
| 1.6 | 0.85 (0.84) | 0.66 (0.68) | 0.51 (0.49) | 0.58 |

**Paper's rule** (the rows each method selects, Salehi-sized samples):
- KARST 0.70 (0.65 to 0.75), alfak 0.64 (0.59 to 0.68); the true landscape on each method's rows is 0.87 in both cases.
- Selection narrows the gap from 0.16 to about 0.07. On Salehi too it is selection that determines the two-time-point gap, but in the simulations the direction does not reverse.

**Daughter at the next sampling (about 7 days):** KARST minus alfak is +0.07 (0.03 to 0.11), +0.14 (0.09 to 0.18) with noise held fixed.

**Spearman of forecast step size with noise** (answer to Q3): true landscape 0.52, KARST 0.49, alfak 0.42. alfak's step size does not follow noise more.

**Landscape accuracy of the same batch of fits** (median r, KARST / alfak, and the fraction where KARST is higher):

| Wavelength | r | Frequent karyotypes | One-step neighbours | Two-step neighbours | KARST higher (r) |
|---|---|---|---|---|---|
| All | 0.59 / 0.37 | 0.95 / 0.88 | 0.67 / 0.46 | 0.58 / 0.37 | 91% |
| 0.2 | 0.21 / 0.08 | 0.97 / 0.97 | 0.25 / 0.17 | 0.11 / 0.08 | 98% |
| 0.4 | 0.48 / 0.36 | 0.95 / 0.91 | 0.59 / 0.50 | 0.47 / 0.37 | 83% |
| 0.8 | 0.79 / 0.62 | 0.95 / 0.86 | 0.84 / 0.68 | 0.79 / 0.62 | 91% |
| 1.6 | 0.80 / 0.50 | 0.91 / 0.69 | 0.82 / 0.53 | 0.80 / 0.50 | 91% |

- **Answer to Q2:** the order in Fig. 2d agrees with the order of landscape accuracy; the alfak lead seen on Salehi is not reproduced.
- Figure: `analysis/step2/figures/fig2d_truth_salehi_cells.pdf`.

## 6. Comparison with Salehi

| | Salehi (two time points) | Simulations, Salehi cell numbers (step 2) |
|---|---|---|
| Ranking by noise alone | 0.73–0.77 | 0.72 |
| KARST, no selection | 0.29 | 0.74 |
| alfak, no selection | 0.32 | 0.58 |
| KARST, paper's rule | 0.22 | 0.70 |
| alfak, paper's rule | 0.45 | 0.64 |

The Salehi numbers come from `analysis/fig2d_noise/` and `analysis/salehi/`; "no selection" on Salehi is the 471 rows common to both.

- The size of the noise (measured by "ranking by noise alone") is about the same on both sides. But in the simulations both methods rank much better than on Salehi, and KARST leads clearly.
- So the main reason that both methods rank poorly on Salehi is not sampling noise. More likely reasons (conjecture, not checked):
  - the real data do not follow the GRF/ABM model;
  - the missegregation rate is higher than 5e-5, and the population is more dispersed;
  - the rate of change differs a lot between cell lines.
- One of these can be tested separately: the paper itself has another set of simulations that vary the missegregation rate (`ABM_varyMisseg`). But this needs a separate plan.

**Wording in report Ver10 that needs to change (waiting for your confirmation):**
- the summary table "Fig. 2d ... follows sampling noise";
- "Fig. 2d mostly measures sampling noise";
- the title of the `#fig2dnoise` page and "cannot tell a better landscape from a worse one on these data".

Suggested replacement:
- On Salehi both methods rank poorly, and noise alone ranks better; with noise held fixed there is no reliable difference.
- On simulations at Salehi cell numbers, this test can tell landscapes apart, with the order true landscape > KARST > alfak.
- So alfak being higher on Salehi cannot be attributed to noise; the two-time-point gap comes from selection, and the cause of the all-time-point gap is unknown.

## 7. Scale and cost

- **CPU:** 540 CPU hours in total.
  - Step 2 fits 331 hours (trial run 21 + the rest 310);
  - step 1 per-row computation 118 hours;
  - step 2 per-row computation 91 hours;
  - the per-row computation is spent almost entirely on the true-landscape forecasts.
- **Wall clock:** 10:20 to 11:52.
- **Memory:** fits at most 5.3 GB (trial run at most 8.0 GB); per-row computation at most about 0.5 GB.
- **Jobs:** setup 22344954; trial runs 22345164 / 22345166 / 22345168; full 22346354 / 22346356 / 22346852; resubmission 22355028. Every job had an add-hours monitor.

## 8. Files

**Code (not committed):**
- `scripts/fig2d_truth.R`: cell-number table, task table, step 0 check, per-row computation per simulation.
- `scripts/fig2d_truth_analyse.R`: analysis and figures.
- `R/abm.R`: `ABM_CELLS` subsampling, `abm_subsample_x`, `abm_sim_seed`.
- `R/config.R`: `ABM_FORECAST_HORIZON`.
- `R/metrics_abm.R`: forecast scoring is computed by the number of forecast rows; results at the default 10 sampling times are unchanged.
- `hpc/env.sh`: passes the new variables into the image.
- `hpc/f2t_job.sh`: job wrapper.

**Results:** `results_6ae2a24/fig2d_truth/`, both on HPC and locally:
- `salehi_cells.tsv`, `tasks/`, `engine_check.tsv`, `jobs.tsv`;
- `step1/rows/`, `step2/rows/`: one file per simulation;
- `analysis/step{1,2}/`: `fig2d.tsv`, `differences.tsv`, `accuracy.tsv`, `figures/`.

On HPC there are also `step2/fits/` (the fits and forecasts of step 2), `step{1,2}/truth/` (the true-landscape forecasts), `code/` (staged code) and `logs/`.
