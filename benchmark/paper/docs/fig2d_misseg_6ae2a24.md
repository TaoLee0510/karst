# Fig. 2d on simulations with different missegregation rates (ABM_varyMisseg): round 6ae2a24

2026-11-23. The data are the paper's 125 `ABM_varyMisseg` simulations (λ = 0.8, missegregation rate 5e-5 to 8e-4). The input at each time point was subsampled to Salehi-like cell counts and then refit: KARST with this round's library (1.1.3), alfak with dev2@f93b2d2. This is a follow-up check to `docs/fig2d_truth_6ae2a24.md`, and the method is the same as its step 2.

## Conclusions

1. **With β fixed at 5e-5 (as done on Salehi), the more dispersed the population, the faster KARST's Fig. 2d falls compared with alfak. At 2e-4 and above, KARST is no longer above alfak.**
   All of the following are without selection, with the daughter at about 15 days, and with Salehi-sized samples as observations:

   | Missegregation rate | KARST | alfak | KARST minus alfak |
   |---|---|---|---|
   | 5e-5 | 0.72 | 0.31 | +0.41 (0.18 to 0.56) |
   | 1e-4 | 0.58 | 0.51 | +0.08 (−0.17 to 0.26) |
   | 2e-4 | 0.24 | 0.27 | −0.03 (−0.26 to 0.14) |
   | 4e-4 | 0.05 | 0.18 | −0.13 (−0.29 to 0.04) |
   | 8e-4 | −0.10 | 0.17 | −0.27 (−0.50 to −0.05) |

   - The median dispersion of the Salehi passages falls between 1e-4 and 2e-4 (Section 3).
   - At 2e-4, the values of the two methods (0.24 / 0.27) are close to those on Salehi (0.30 / 0.32).
   - So the pattern on Salehi, "both methods reach only about 0.3, with alfak slightly higher", is reproduced in simulations with similar dispersion and β fixed at 5e-5.

2. **The reversal comes from KARST's fit: when β is lower than the true missegregation rate, the fitted landscape gets worse. It does not come from using the wrong rate in the forecast, nor from sampling noise.**
   - **When fitting and forecasting with the true missegregation rate**, KARST is not below alfak at any level: 0.27 vs 0.23 at 4e-4, 0.30 vs 0.14 at 8e-4. However, at 1e-4 and above, the intervals of the difference all include 0.
   - **Landscape accuracy** (median r, 1e-4 / 4e-4 / 8e-4):

     | | KARST | alfak |
     |---|---|---|
     | β fixed 5e-5 | 0.81 / 0.59 / 0.20 | 0.64 / 0.50 / 0.32 |
     | β = true missegregation rate | 0.81 / 0.69 / 0.51 | 0.64 / 0.52 / 0.37 |

     With β fixed, KARST is below alfak at 8e-4. With β at the true value, KARST is higher at every level. alfak is almost unaffected by β.
   - **Changing only the forecast rate costs little.** The true landscape forecast at 5e-5 is at most 0.08 lower than when forecast at the true missegregation rate (4e-4).
   - **The 4e-4 level has to be looked at separately.** With β fixed, KARST's r is still higher than alfak's (0.59 vs 0.50), yet its raw Fig. 2d value is lower (0.05 vs 0.18). With noise held fixed, KARST is again above alfak (+0.11, interval includes 0).
     - The reason is that at this level alfak's predicted step size follows the noise, while KARST's does not: their Spearman with noise is 0.18 and −0.04 respectively.
     - And here, noise alone already ranks the observed changes at 0.84.
     - So part of alfak's raw advantage at this level comes from noise.

3. **The more dispersed the population, the closer Fig. 2d comes to measuring only sampling noise; this is also true for the true landscape.**
   - Ranking by noise only: rises from 0.66 at 5e-5 to 0.84 at 4e-4.
   - True landscape: falls from 0.79 to 0.59. With noise held fixed, only 0.27 and 0.33 remain at 4e-4 and 8e-4.
   - When observations use 1000 cells instead, the true landscape still reaches 0.87 at 4e-4. So the low values come from Salehi-like small samples, not from the landscape itself.
   - Therefore, in the dispersed cell lines on Salehi, Fig. 2d can hardly distinguish even the true landscape; it is a weak test.

4. **When the population is dispersed, the karyotype list covers fewer cells, and fits also fail more often.**
   - Fraction of parent cells on the karyotype list: 100% at 5e-5, 93% at 2e-4, 69% at 4e-4, 33% at 8e-4.
   - At 8e-4, about 30% of fits report "no frequent karyotypes" (minobs 10 and 20), the same for both methods.
   - In the dispersed cell lines on Salehi, too, only 26%–73% of parent cells are on alfak's list.

5. **Relation to the previous check (`fig2d_truth`):**
   - "On simulations, Fig. 2d follows landscape accuracy" still broadly holds: the fall in KARST's Fig. 2d goes together with the fall in its r.
   - The previous check said that alfak's predicted step size does not follow noise more. This holds in simulations that are not dispersed; but at 4e-4 and 8e-4 with β fixed, alfak's step size does follow noise more (0.18 vs −0.04, 0.18 vs −0.08). So "noise helped alfak" partly holds in dispersed populations.
   - The previous check guessed that both methods reach only about 0.3 on Salehi because "the fitted landscapes are far from the true situation on Salehi". Here a case that reproduces it was found in simulations: the population is more dispersed than at 5e-5, while β is fixed at 5e-5.
   - This is still simulation. The dispersion of the Salehi cell lines differs widely, and no check was done cell line by cell line.

6. **The β question:** the "true missegregation rate" here is only a control within the simulations. It shows that KARST is affected more than alfak when β is too low, but it does not provide a way to determine β on real data. Per the earlier decision, β stays fixed at 5e-5, and only externally given values are accepted.

7. **Report:** neither Ver10 nor Ver11 has been changed; the decision waits until you have read this document (Section 8).

## 1. Question

The previous check (`docs/fig2d_truth_6ae2a24.md`) found on the paper's main simulation set (missegregation rate 5e-5) that the Fig. 2d order is true landscape > KARST > alfak; the pattern on Salehi of alfak being slightly higher was not reproduced.

But the Salehi passages are much more dispersed than the main simulation set: the most common karyotype takes a lower share, and there are more karyotypes seen only once (Section 3). The paper also has a set of simulations that change only the missegregation rate, `ABM_varyMisseg`, with which the effect of dispersion can be tested.

Questions to answer:
- **Q1:** When the population is more dispersed, does Fig. 2d fall? Does the order of KARST and alfak reverse?
- **Q2:** If it reverses, what is the cause: sampling noise, landscape accuracy, a mismatch between β and the true missegregation rate, or fewer cells covered by the karyotype list?

## 2. Method

**Simulations:**
- The paper's `ABM_varyMisseg`, 125 in total: 25 each at missegregation rate 5e-5, 1e-4, 2e-4, 4e-4, 8e-4, all with λ 0.8.
- Each time point stores a sample of 1000 cells; the sampling times are the same as in the main simulation set.
- The 5e-5 files here are different simulations from the files of the same name in the main simulation set.

**Fits (same as step 2 of `fig2d_truth`):**
- The cell count at each time point is drawn from the cell counts of the Salehi two-time-point inputs, median 389. The seed is fixed by the simulation and now also includes the missegregation rate (`abm_sim_seed` in `R/abm.R`; the seeds of the main simulation set are unchanged).
- Designs: the three two-time-point designs w2, w4, w8, each fit at minobs 5, 10 and 20.
- Forecasts: starting from the subsampled parent, to the next 3 sampling times, 5 replicates.

**Two β arms:**
- **β fixed 5e-5:** both fit and forecast use 5e-5, as on Salehi. 125 simulations.
- **β = true missegregation rate:** both fit and forecast use the simulation's true missegregation rate. 100 simulations; the 5e-5 level is the same as in the previous arm and is reused directly.

**Forecasts present in every row:**
- KARST, alfak;
- the true landscape, forecast with the true missegregation rate;
- the true landscape, forecast with 5e-5;
- noise only.

The true-landscape forecasts use the same karyotype list (alfak's list at that minobs) and are shared by both β arms.

**Statistics:**
- Same as the previous check: raw Spearman; partial Spearman with noise held fixed; ranking by noise only; Spearman of predicted step size with noise.
- Intervals: resampling by simulation within each missegregation-rate level, 1000 bootstrap replicates.
- Each row = simulation × design × minobs. Only rows that have all forecasts are used.
- Two kinds of selection: no selection, and the paper rule (for each simulation and each design, take the minobs with the highest cross-validation score).

**Dispersion:**
- Draw 300 cells each from parent and daughter, 10 draws, averaged. Parent = the last sample before day 120; daughter = the second sample after it.
- Four quantities:
  - Simpson diversity;
  - the share of the most common karyotype;
  - the fraction of cells in karyotypes seen only once;
  - the fraction of daughter cells in karyotypes the parent does not have.
- For Salehi, non-overlapping passage pairs in which both passages have at least 300 cells are used, 109 pairs in total.

## 3. Dispersion

| Data | Simpson | Most common karyotype | Seen only once | New in daughter |
|---|---|---|---|---|
| Main simulation set λ 0.8 (5e-5) | 0.80 | 0.38 | 0.08 | 0.11 |
| varyMisseg 5e-5 | 0.78 | 0.42 | 0.08 | 0.11 |
| varyMisseg 1e-4 | 0.92 | 0.21 | 0.17 | 0.19 |
| varyMisseg 2e-4 | 0.96 | 0.14 | 0.31 | 0.36 |
| varyMisseg 4e-4 | 0.99 | 0.056 | 0.57 | 0.61 |
| varyMisseg 8e-4 | 0.995 | 0.023 | 0.80 | 0.83 |
| Salehi (109 pairs, median) | 0.90 | 0.21 | 0.22 | 0.29 |

- varyMisseg 5e-5 has the same dispersion as the main simulation set at λ 0.8.
- The Salehi median falls between 1e-4 and 2e-4, but the cell lines differ widely:

| Cell line | Simpson | Most common karyotype | Seen only once | Parent cells on alfak list | Closest level |
|---|---|---|---|---|---|
| SA039U | 0.42 | 0.65 | 0.03 | 0.99 | below 5e-5 |
| SA1035U | 0.43 | 0.75 | 0.08 | 0.95 | below 5e-5 |
| SA1035T | 0.77 | 0.47 | 0.19 | 0.95 | 5e-5 to 1e-4 |
| SA532 | 0.84 | 0.35 | 0.26 | 0.88 | 5e-5 to 1e-4 |
| SA906a / SA906b | 0.89 / 0.83 | 0.21 / 0.30 | 0.19 / 0.17 | 0.92 / 0.94 | about 1e-4 |
| SA535 | 0.99 | 0.053 | 0.63 | 0.65 | about 4e-4 |
| SA609 | 0.985 | 0.074 | 0.58 | 0.73 | about 4e-4 |
| SA609H2, R2, R2b | 0.99 | 0.06–0.07 | 0.77–0.78 | 0.48–0.83 | 4e-4 to 8e-4 |
| SA000, SA001, SA609H1 | 0.995 | 0.016 | 0.81–0.82 | 0.34–0.60 | about 8e-4 |
| SA609TIIDT | 0.997 | 0.003 | 1.00 | 0.26 | above 8e-4 |

- What is used here is the shape of the dispersion, not an estimate of the Salehi missegregation rate. The Salehi dispersion may also have other causes, such as weaker selection or a mixture of different subpopulations.
- Figure: `analysis/figures/dispersion_misseg.pdf`.

## 4. Checks

**Fit failures** (within each β arm, out of the number of fits of KARST and of alfak):

| | β fixed 5e-5 | β true value |
|---|---|---|
| KARST failed | 72 / 1,125 | 72 / 900 |
| alfak failed | 85 / 1,125 | 82 / 900 |

- The vast majority are "no frequent karyotypes": 69 per method in each arm, almost all at 8e-4 with minobs 10 and 20. It is the same for both methods, because here almost no karyotype appears minobs times.
- Other failures:
  - alfak reports "Regression matrix ... colinear": 16 with β fixed, 13 with β at the true value;
  - KARST reports "A timepoint has no cell inside the model domain": 3 in each arm.
- One alfak fit (1e-4 rep 19, w8, minobs 5, β fixed) exceeded 20 GB of memory and was killed. Rerun with 48 GB, it succeeded, with a peak of 27.5 GB and a run time of 55 minutes; the per-row results of this simulation were then recomputed.
- Only rows with all forecasts enter the comparison; the fractions are in Section 7.

**A bug in the distance computation (fixed):**
- The default networkflow algorithm of `transport::wasserstein` never finishes on a few inputs. Example: 111 parent karyotypes against a KARST forecast of 13,692 karyotypes did not finish in 1.5 hours and gave no error either.
- Fix (`wdist` in `scripts/fig2d_truth.R`): larger problems are first computed with networkflow in a child process; if this takes more than 60 seconds, revsimplex is used instead.
- revsimplex is another exact solver for the same optimal transport problem: the stuck example finished in 3 seconds and gave the same value (0.30297). For calls that did not get stuck, the results are bit-for-bit identical to before.
- Across all results, revsimplex was used 7 times, in 7 simulations. Results completed earlier are not affected.
- The same default algorithm is also used in this round's scoring scripts (`R/metrics_abm.R`, `R/metrics_salehi.R`) and in `scripts/fig2d_noise.R`. These all finished earlier, so they did not run into this problem, and they were not changed this time.

**Tasks:**
- On node hpctpa3pc0003, the containers of 5 tasks failed to start ("unknown userid"); there was also 1 on hpctpa3pc0002. They wrote nothing and completed after resubmission.
- hpctpa3pc0003 has been added to the exclusion list.

## 5. Results: β fixed 5e-5 (as on Salehi)

**No selection, daughter at about 15 days, Salehi-sized samples as observations** (values in parentheses are with noise held fixed):

| Missegregation rate | True landscape | True landscape, forecast at 5e-5 | KARST | alfak | Noise only |
|---|---|---|---|---|---|
| 5e-5 | 0.79 (0.81) | 0.79 (0.81) | 0.72 (0.67) | 0.31 (0.23) | 0.66 |
| 1e-4 | 0.84 (0.88) | 0.83 (0.88) | 0.58 (0.58) | 0.51 (0.45) | 0.56 |
| 2e-4 | 0.65 (0.73) | 0.63 (0.72) | 0.24 (0.29) | 0.27 (0.25) | 0.67 |
| 4e-4 | 0.59 (0.27) | 0.52 (0.23) | 0.05 (0.16) | 0.18 (0.05) | 0.84 |
| 8e-4 | 0.61 (0.33) | 0.63 (0.42) | −0.10 (−0.06) | 0.17 (0.06) | 0.79 |

**KARST minus alfak:**

| Missegregation rate | Raw | Noise held fixed |
|---|---|---|
| 5e-5 | +0.41 (0.18 to 0.56) | +0.44 (0.18 to 0.62) |
| 1e-4 | +0.08 (−0.17 to 0.26) | +0.14 (−0.08 to 0.32) |
| 2e-4 | −0.03 (−0.26 to 0.14) | +0.04 (−0.22 to 0.26) |
| 4e-4 | −0.13 (−0.29 to 0.04) | +0.11 (−0.11 to 0.33) |
| 8e-4 | −0.27 (−0.50 to −0.05) | −0.12 (−0.38 to 0.09) |

- With noise held fixed, the difference at 4e-4 turns positive again. The reason is that alfak's predicted step size follows noise more than KARST's does: at 4e-4, the Spearman of predicted step size with noise is 0.18 for alfak and −0.04 for KARST.
- **True landscape minus KARST** grows with the missegregation rate: +0.07 at 5e-5, +0.25 at 1e-4, +0.41 at 2e-4, +0.55 at 4e-4, +0.70 at 8e-4.

**Paper rule** (rows selected by each method, raw values, KARST / alfak):

| Missegregation rate | 5e-5 | 1e-4 | 2e-4 | 4e-4 | 8e-4 |
|---|---|---|---|---|---|
| KARST / alfak | 0.71 / 0.31 | 0.58 / 0.56 | 0.22 / 0.24 | 0.08 / 0.23 | −0.09 / 0.33 |

**Daughter taken at the next sample (about 7 days), KARST minus alfak:**
- +0.27 (0.01 to 0.44) at 5e-5;
- −0.14 (−0.36 to 0.05) at 2e-4;
- −0.38 (−0.59 to −0.17) at 4e-4;
- −0.17 (−0.38 to 0.02) at 8e-4.

**Observations with 1000 cells instead** (raw values, true landscape / KARST / alfak):

| Missegregation rate | 5e-5 | 1e-4 | 2e-4 | 4e-4 | 8e-4 |
|---|---|---|---|---|---|
| True landscape / KARST / alfak | 0.96 / 0.78 / 0.45 | 0.96 / 0.68 / 0.55 | 0.86 / 0.40 / 0.32 | 0.87 / 0.22 / 0.36 | 0.63 / −0.02 / 0.21 |

- Figures: `analysis/figures/fig2d_misseg.pdf` (left column), `fig2d_misseg_1000_cells.pdf`.

## 6. Results: β = true missegregation rate (control within the simulations)

**No selection, daughter at about 15 days, Salehi-sized samples as observations** (values in parentheses are with noise held fixed; the two true-landscape columns are the same as in Section 5):

| Missegregation rate | KARST | alfak | KARST minus alfak (raw) |
|---|---|---|---|
| 5e-5 | 0.72 (0.67) | 0.31 (0.23) | +0.41 (0.18 to 0.56) |
| 1e-4 | 0.62 (0.62) | 0.50 (0.44) | +0.12 (−0.11 to 0.28) |
| 2e-4 | 0.38 (0.42) | 0.28 (0.25) | +0.10 (−0.15 to 0.32) |
| 4e-4 | 0.27 (0.16) | 0.23 (0.08) | +0.05 (−0.17 to 0.24) |
| 8e-4 | 0.30 (0.15) | 0.14 (0.01) | +0.16 (−0.11 to 0.44) |

- Compared with β fixed, KARST is higher by 0.14, 0.22 and 0.40 at 2e-4 to 8e-4 respectively; alfak barely changes.
- With observations of 1000 cells, KARST is 0.78 / 0.73 / 0.55 / 0.54 / 0.33 and alfak is 0.45 / 0.54 / 0.32 / 0.35 / 0.17. KARST is higher at every level.
- Paper rule (KARST / alfak): 0.71 / 0.31, 0.64 / 0.56, 0.38 / 0.28, 0.36 / 0.24, 0.31 / 0.42. At 8e-4, alfak selected only 66 rows.
- With the daughter at the next sample (about 7 days), KARST minus alfak at 2e-4 to 8e-4 is −0.12, −0.09, −0.04, and all intervals include 0.
- Figure: `analysis/figures/fig2d_misseg.pdf` (right column).

## 7. Landscape accuracy and the karyotype list

**Landscape accuracy of the same fits** (median r, KARST / alfak; in parentheses, the fraction in which KARST is higher):

| Missegregation rate | β fixed 5e-5 | β true value |
|---|---|---|
| 5e-5 | 0.76 / 0.58 (91%) | same as left |
| 1e-4 | 0.81 / 0.64 (91%) | 0.81 / 0.64 (96%) |
| 2e-4 | 0.77 / 0.63 (91%) | 0.79 / 0.64 (97%) |
| 4e-4 | 0.59 / 0.50 (77%) | 0.69 / 0.52 (90%) |
| 8e-4 | 0.20 / 0.32 (38%) | 0.51 / 0.37 (63%) |

- On frequent karyotypes, KARST is higher in both arms: at 8e-4, 0.79 vs 0.66 (β fixed) and 0.83 vs 0.62 (β true value).
- A too-low β mainly harms the part outside the frequent karyotypes.
- Figure: `analysis/figures/accuracy_misseg.pdf`.

**Fits, complete rows, cells covered by the karyotype list** (β fixed; almost the same with β at the true value):

| Missegregation rate | KARST succeeded | alfak succeeded | Complete rows | Parent cells on list | Daughter cells on list | List size (median) |
|---|---|---|---|---|---|---|
| 5e-5 | 100% | 98.7% | 98.7% | 100% | 99.7% | 7,491 |
| 1e-4 | 100% | 98.7% | 98.7% | 98.6% | 98.3% | 9,207 |
| 2e-4 | 99.1% | 97.8% | 98.7% | 93.0% | 90.5% | 9,587 |
| 4e-4 | 98.7% | 98.2% | 99.5% | 68.5% | 64.6% | 6,864 |
| 8e-4 | 70.2% | 68.9% | 96.2% | 32.7% | 31.0% | 2,741 |

"Complete rows" means the fraction, among rows in which both methods fit successfully, that have all forecasts.

## 8. Comparison with Salehi

| | Salehi (two time points) | Simulation 2e-4, β fixed | Simulation 4e-4, β fixed | Simulation 2e-4, β true value |
|---|---|---|---|---|
| Ranking by noise only | 0.75 | 0.67 | 0.84 | 0.67 |
| KARST, no selection | 0.30 | 0.24 | 0.05 | 0.38 |
| alfak, no selection | 0.32 | 0.27 | 0.18 | 0.28 |
| KARST, noise held fixed | 0.47 | 0.29 | 0.16 | 0.42 |
| alfak, noise held fixed | 0.44 | 0.25 | 0.05 | 0.25 |

The Salehi numbers come from `analysis/fig2d_noise/fig2d_noise.tsv` (two time points, per trajectory).

- At 2e-4, where the dispersion is close to the Salehi median, with β fixed the raw values of both methods are close to Salehi's, and alfak is slightly higher.
- But with noise held fixed, both are 0.44–0.47 on Salehi, higher than the 0.25–0.29 in the simulations. So the simulations do not fully reproduce Salehi.
- Salehi is a mixture of different cell lines (Section 3): some look like 5e-5, some like 8e-4. Comparing it with a single level shows only the rough picture.

**Places in report Ver10 that may need additions or changes** (waiting for your decision; nothing has been changed yet):
- The change suggested by the previous document (`fig2d_truth`) was: "On simulations at Salehi cell counts, Fig. 2d can tell good landscapes from bad ones, with the order true landscape > KARST > alfak". A condition now has to be added: this holds only when the population is not very dispersed, or when β matches the true missegregation rate.
- A sentence could be added: on simulations with dispersion similar to Salehi and β fixed at 5e-5, the Fig. 2d of both methods falls to about 0.25 and KARST no longer leads; with the true missegregation rate KARST still leads. So alfak being slightly higher on Salehi is consistent with "β is too low for the dispersed cell lines".
- This points in the same direction as the earlier Fig. 4d finding: using the same β for all Salehi cell lines biases KARST's results to one side.

## 9. Scale and cost

- **CPU:** 447 CPU hours in total.
  - Fits 368 hours: trial run 37, the rest 330, reruns 1;
  - Per-row computation 79 hours. The β-fixed arm had to compute the true-landscape forecasts, which took most of this; the β-true arm reused these forecasts directly and took only 4 hours.
- **Wall clock:** 13:20 to 18:11.
- **Memory:**
  - Fits at most 20 GB (trial run 11.9 GB); one task exceeded 20 GB and used 27.5 GB on rerun (Section 4);
  - Per-row computation at most 1.1 GB.
- **Jobs:**
  - setup: 22355424;
  - trial fits: 22355426, 22355428;
  - trial per-row: 22355430, 22355432;
  - remaining fits: 22355823, 22355825;
  - remaining per-row: 22355827, 22355829;
  - reruns: 22371022, 22371024, 22371026, 22371028, 22371286, 22371373.
  - Every job had an add-hours monitor.

## 10. Files

**Code (not committed):**
- `scripts/fig2d_truth.R`:
  - `--source=misseg`;
  - a new `dispersion` stage;
  - two β arms (`--arm`);
  - true-landscape forecasts shared by both arms;
  - the timeout and revsimplex fallback in `wdist`.
- `scripts/fig2d_truth_analyse.R`: with `F2T_SOURCE=misseg`, groups by missegregation rate, adds a "true landscape, forecast at 5e-5", and draws four figures.
- `scripts/run_task.R`: with `MISSEG_FORECAST=TRUE`, the missegregation-rate simulations are also forecast.
- `R/abm.R`: the seed includes the missegregation rate.
- `hpc/env.sh`, `hpc/f2t_job.sh`: pass the new variables and provide the misseg steps.

**Results:** `results_6ae2a24/fig2d_misseg/`, both on HPC and locally:
- `dispersion.tsv`, `jobs.tsv`;
- `beta_fixed/rows/`, `beta_true/rows/`: one file per simulation;
- `analysis/beta_fixed/`, `analysis/beta_true/`: `fig2d.tsv`, `differences.tsv`, `accuracy.tsv`, `fits_and_list.tsv`;
- `analysis/figures/`: `fig2d_misseg.pdf`, `fig2d_misseg_1000_cells.pdf`, `accuracy_misseg.pdf`, `dispersion_misseg.pdf`.

Also on HPC:
- `salehi_cells.tsv`, `tasks/`;
- `beta_{fixed,true}/fits/` (fits and forecasts);
- `truth/` (true-landscape forecasts, shared by both arms);
- `code/` (staged code, including the submit script `submit_f2m.sh`) and `logs/`.
