# Round report for KARST 6ae2a24 on the own benchmark (100 held-out landscapes per wavelength)

**Basic information:**
- Date: 2026-11-21.
- Method versions: KARST 6ae2a24 (1.1.3); alfak dev2@f93b2d2 (minobs 5, 10, 20).
- Code for this round: 942182c.
- Results folder: `benchmark/results/sparse_two_timepoint_6ae2a24/` inside the karst clone on HPC.
- Previous round: e7e32f8 (KARST 1.0.0, 3 held-out landscapes per wavelength), folder `benchmark/results/sparse_two_timepoint/`.

> **Terms:**
> - **"Dataset"**: one sampling design of the scored patient (patient 1) on one held-out landscape.
> - **"Design"**: the cell counts at the two time points and the interval between them. There are 8 core designs and 6 stress designs.
> - **"Distance"**: the number of steps to the nearest observed karyotype. Distances 0–2 form the main table, 3–4 the far-distance table.
> - **"Error ratio"**: RMSE ÷ the RMSE of the null prediction. The null prediction is "all karyotypes have the same fitness"; a value below 1 means better than the null prediction.
> - **"Direction accuracy"**: among pairs of karyotypes that differ by one copy, the fraction in which the predicted direction of the fitness change is the same as the truth. Ties count as 1/2; random guessing gives 1/2.
> - **"Coverage"**: the fraction of true values that fall inside the reported 95% interval.

## 0. What this round did

**Larger held-out set:** from 3 held-out landscapes per wavelength to 100, 500 in total over 5 wavelengths.
- L16–L500 are newly generated (seeds 82016–82500);
- the datasets and nodes of L01–L15 are copied from e7e32f8.

**Fits:**
- **KARST:** both initial states, default (balance prior) and independent, were all refit, on both core and stress designs: 14000 for distances 0–2 and 8000 for distances 3–4.
- **alfak:** fit with minobs 5, 10, 20 on the core designs of L16–L500; the alfak fits of L01–L15 are copied from e7e32f8.
- **Null prediction:** one per dataset.

**Scoring and summaries:**
- Scoring and the two summaries (`own_summary.json`, `own_coverage.json`) were already completed in this round's folder; report pages Ver8 and Ver9 use exactly these.
- This report adds: the completeness check, the one-by-one comparison with e7e32f8, the independent method, direction accuracy and coverage split by design, and the stress designs.
- The additional computations only read this round's folder; their results and scripts are not in the repository, see Section 11.
- The recomputed default and alfak error ratios, Spearman and direction accuracy are identical to `own_summary.json` (maximum difference 0); coverage differs from `own_coverage.json` by at most 4e-16.

**Scope:**
- All truth simulations and fits use missegregation rate 5e-5, so this round measures accuracy when the missegregation rate is known.
- The data contain no dynamics of the clone-replacement kind.

## 1. Conclusions

1. **This round is complete (Section 2).**
   - Every expected dataset × method has a fit file and complete score rows: 33000 in the main table, 24000 in the far-distance table.
   - KARST has not a single failed fit (default and independent 7000 + 4000 each).
   - The alfak failures are all its own:
     - minobs 5 / 10 / 20 fail 110 / 145 / 253 each on every table;
     - the causes are "no high-frequency karyotype" (1 / 31 / 135) and "collinear regression matrix" (109 / 114 / 118);
     - by convention, failed fits are scored as the null prediction.
2. **From 1.0.0 to 1.1.3 only 7 fits changed (Section 3).**
   - The 210 datasets of L01–L15 are identical to e7e32f8, and so are the alfak fits.
   - The independent fits are all identical.
   - 13 default fits differ from before:
     - 6 differ only below 1e-7;
     - 7 have real differences; the largest is L06 turned_3000_100, with a prediction difference of 0.074.
3. **Numeric accuracy: KARST is better than the null prediction and alfak in every cell (Section 4).**
   - The default error ratio is 0.32–0.75 at distances 0–2 and 0.28–0.59 at distances 3–4; the upper limits of all 95% intervals are below 0.76.
   - alfak is worse than the null prediction in every wavelength × distance cell, error ratio 1.37–81.
   - By dataset, default is better than the null prediction in 95–100% of datasets and better than alfak's best setting in 92–100% of datasets.
4. **Direction accuracy: KARST 0.67–0.99 (distances 0–2) and 0.57–0.81 (distances 3–4); alfak 0.29–0.55 (Section 5).**
   - alfak is below 1/2 in 12 of the 15 distance 0–2 cells, that is, it calls more directions wrong than right.
5. **Numeric accuracy and direction accuracy part ways in two kinds of situation (Section 10):**
   - **Rugged landscapes (λ 0.2), far distances:** the error ratio is the best (0.28), but direction is close to random (0.57 at distance 4), and Spearman is only 0.02.
   - **Smooth landscapes (λ ≥ 1.2), observed karyotypes:** both are on the weak side, error ratio 0.67–0.75, direction 0.67–0.74.
6. **default is the better of the two KARST methods (Sections 4, 5).**
   - Direction accuracy is the same as independent (differences at most 0.02).
   - Numeric accuracy is clearly better at λ 0.4–1.2, for example 0.56 vs 0.77 at λ 0.8 distance 1; by design cell, better in 110 / 120 cells.
   - This round shows no trade-off between the two accuracies.
7. **Coverage (Section 7):**
   - core designs: 0.954–0.998 in every wavelength × distance cell; 0.42–0.77 with the model sd alone.
   - The misses outside the interval are almost all on the outer side (the true value is farther from the reference than the estimate); inner misses are at most 0.6% in every cell.
   - Over core and stress designs together there are 24.9 million intervals (distances 0–2), and **not a single one excludes 0**.
   - In the stress designs, when the landscape changes between the two time points (shift, turned), coverage at distance 0 drops to 0.88–0.91, and to only 0.75 for the worst wavelength.
8. **Cost (Section 9):** the median default fit takes 32 seconds (main table) and 24 seconds (far distance), 4–5 times independent. The KARST fits of the whole round took 250 CPU hours in total; alfak minobs 5 took 590 CPU hours in total.

## 2. Completeness

| Table | Method | Expected | With file | Scores complete | Succeeded | Failed | Median per fit (s) | CPU hours |
|---|---|---|---|---|---|---|---|---|
| Main | default | 7000 | 7000 | 7000 | 7000 | 0 | 31.7 | 136.8 |
| Main | independent | 7000 | 7000 | 7000 | 7000 | 0 | 6.1 | 25.4 |
| Main | alfak 5 / 10 / 20 | 4000 each | 4000 | 4000 | 3890 / 3855 / 3747 | 110 / 145 / 253 | 17.8 / 6.9 / 4.5 | 309.8 / 44.5 / 15.3 |
| Main | Null prediction | 7000 | 7000 | 7000 | 7000 | 0 | 0 | 0 |
| Far distance | default | 4000 | 4000 | 4000 | 4000 | 0 | 24.1 | 73.6 |
| Far distance | independent | 4000 | 4000 | 4000 | 4000 | 0 | 5.5 | 14.8 |
| Far distance | alfak 5 / 10 / 20 | 4000 each | 4000 | 4000 | 3890 / 3855 / 3747 | 110 / 145 / 253 | 20.8 / 8.2 / 5.3 | 280.1 / 53.8 / 17.9 |
| Far distance | Null prediction | 4000 | 4000 | 4000 | 4000 | 0 | 0 | 0 |

- alfak fits only the core designs; the far-distance table also has only core designs.
- For the main table, "scores complete" means all three layers, distances 0, 1, 2, are present; for the far-distance table it means both layers, distances 3, 4, are present.

**Jobs:** the final state of every job in the job register was checked.
- **Out of memory (OUT_OF_MEMORY):** on the large datasets of alfak minobs 5 (100/3000, 3000/100), 22 tasks in total ran out of memory at 32G:
  - 13 in the main table, 22173673; 9 in the far-distance table, 22173679.
  - All were rerun to completion at 64G: 22226038, 22236094 and 22226040, 22236096.
  - Of these, task 602 of 22236094 failed at container start ("unknown userid", node hpctpa3pc0022); it was rerun with that node excluded (22236306) and completed.
- **Cancelled jobs:** only the scoring jobs 22173681, 22173683, 22226042, 22226044; their dependencies could never be satisfied, and they were replaced one after another by the split scoring.
- **KARST library:** the logs of all 14000 main-table tasks state `karst 1.1.3 …/lib/karst_6ae2a24` and code 942182c. The far-distance logs do not print the library, but that table uses the same `hpc/env.sh` as the main table, which sets the library path from `karst.lock` (6ae2a24); all 8000 logs state code 942182c.

## 3. One-by-one comparison with e7e32f8 (L01–L15)

- **Datasets:** the count matrices of the 210 scored datasets are identical.
- **alfak:** the fits are identical one by one (copied), and so are the failures.
- **independent:** the 168 (main table) and 96 (far distance) fits that have a counterpart are identical in status, karyotypes and predictions.
  - e7e32f8 has no independent fits for L13–L15 (λ 0.2), and no stress designs for these three landscapes.
- **default:**
  - of the 192 main-table fits with a counterpart, 184 are identical;
  - of the 120 in the far-distance table, 115 are identical;
  - the statuses are all identical.
  - 6 differ by less than 1e-7, which is numerical rounding.
  - The 7 with real differences:

| Table | Dataset | Max prediction difference | Max sd difference |
|---|---|---|---|
| Main | L06 turned_3000_100_360d | 0.074 | 0.084 |
| Main | L01 shift_3000_100_360d | 0.0058 | 0.0048 |
| Main | L11 turned_500_700_360d | 0.0043 | 0.0036 |
| Main | L11 long_500_700_1080d | 0.0023 | 0.0020 |
| Main | L09 sparse_200_200_180d | 0.0018 | 0.0017 |
| Far distance | L09 sparse_200_200_180d | 0.0021 | 0.0019 |
| Far distance | L11 long_500_700_1080d | 0.0020 | 0.0016 |

These differences come from the changes between 1.0.0 and 1.1.3. This report did not trace, one by one, which change caused each.

## 4. Numeric accuracy: error ratio (core designs)

Averaged first within the 8 designs, then across the 100 landscapes. The brackets after default are 95% intervals from a landscape cluster bootstrap (2000 resamples, seed 1).

**default:**

| λ | Distance 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0.2 | 0.537 [0.521, 0.554] | 0.389 [0.382, 0.396] | 0.317 [0.311, 0.323] | 0.283 [0.276, 0.291] | 0.276 [0.267, 0.285] |
| 0.4 | 0.528 [0.517, 0.540] | 0.486 [0.478, 0.495] | 0.481 [0.474, 0.487] | 0.463 [0.456, 0.470] | 0.452 [0.445, 0.459] |
| 0.8 | 0.611 [0.603, 0.620] | 0.558 [0.550, 0.568] | 0.570 [0.563, 0.578] | 0.575 [0.567, 0.583] | 0.589 [0.581, 0.596] |
| 1.2 | 0.670 [0.663, 0.677] | 0.559 [0.552, 0.566] | 0.558 [0.553, 0.564] | 0.568 [0.562, 0.574] | 0.583 [0.577, 0.588] |
| 1.6 | 0.746 [0.735, 0.757] | 0.594 [0.581, 0.608] | 0.546 [0.538, 0.555] | 0.530 [0.523, 0.536] | 0.529 [0.522, 0.535] |

**independent:**

| λ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0.2 | 0.536 | 0.408 | 0.335 | 0.293 | 0.277 |
| 0.4 | 0.687 | 0.598 | 0.562 | 0.533 | 0.511 |
| 0.8 | 0.808 | 0.766 | 0.754 | 0.742 | 0.740 |
| 1.2 | 0.724 | 0.651 | 0.670 | 0.686 | 0.701 |
| 1.6 | 0.761 | 0.614 | 0.599 | 0.605 | 0.614 |

**alfak (in each cell the best of minobs 5, 10, 20, chosen after the fact, which favours alfak):**

| λ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0.2 | 4.80 | 3.20 | 2.89 | 2.84 | 2.87 |
| 0.4 | 3.63 | 2.46 | 2.07 | 1.92 | 1.82 |
| 0.8 | 2.45 | 1.82 | 1.55 | 1.44 | 1.37 |
| 1.2 | 2.69 | 1.96 | 1.66 | 1.52 | 1.44 |
| 1.6 | 21.8 | 14.7 | 10.8 | 8.84 | 7.41 |

- **Each alfak setting on its own:** error ratio 1.37–21.8 for minobs 5, 1.46–28.5 for minobs 10, 1.94–81.0 for minobs 20. The lower limits of all bootstrap intervals are above 1.2.
- **By design (core, 120 design × wavelength × distance 0–2 cells):**
  - default is better than the null prediction in all 120 cells (0.29–0.86), and better than all three alfak settings in all 120 cells;
  - independent is better than the null prediction in 116 cells;
  - each alfak setting is better than the null prediction in only 1–2 cells.
- **By dataset (core, 800 per wavelength × distance):**
  - the fraction in which default beats the null prediction is 95.0% (λ 0.2, distance 0) to 100%;
  - the fraction in which it beats alfak's best setting is 92.1% (λ 1.6, distance 0) to 100%;
  - alfak's best setting beats the null prediction in only 18–39%.
- **stress designs (KARST only):** default is better than the null prediction in all 90 design × wavelength × distance cells, error ratio 0.29–0.90. The weakest is λ 1.6, shift_3000_100, distance 0 (0.90).

## 5. Direction accuracy

Brackets are bootstrap 95% intervals; for default the interval width is at most ±0.008 in every cell.

**default:**

| λ | Distance 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0.2 | 0.990 [0.986, 0.995] | 0.850 | 0.739 | 0.640 | 0.574 |
| 0.4 | 0.958 [0.952, 0.963] | 0.877 | 0.836 | 0.789 | 0.755 |
| 0.8 | 0.812 [0.806, 0.819] | 0.825 | 0.832 | 0.804 | 0.785 |
| 1.2 | 0.743 [0.738, 0.747] | 0.820 | 0.821 | 0.810 | 0.796 |
| 1.6 | 0.674 [0.666, 0.681] | 0.789 | 0.800 | 0.804 | 0.795 |

**Difference between independent and default:**
- At most 0.017 at distances 0–2, at most 0.012 at distances 3–4.
- Almost identical at λ 0.2, 0.4; at λ ≥ 0.8 independent is slightly lower in most cells. The only exception is λ 1.6, distance 0, where independent is 0.003 higher.

**alfak (in each cell the best of the three settings):**

| λ | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 0.2 | 0.502 | 0.541 | 0.516 | 0.547 | 0.535 |
| 0.4 | 0.460 | 0.496 | 0.486 | 0.536 | 0.542 |
| 0.8 | 0.394 | 0.433 | 0.416 | 0.456 | 0.467 |
| 1.2 | 0.449 | 0.443 | 0.429 | 0.464 | 0.476 |
| 1.6 | 0.522 | 0.498 | 0.491 | 0.510 | 0.517 |

- **alfak below 1/2:** each of the three settings is below 1/2 in 12 of the 15 distance 0–2 cells. Of these, the cells where the upper limit of the bootstrap interval is also below 1/2 number 12 / 11 / 8 for minobs 5 / 10 / 20.
- **By design (core, 120 cells):**
  - default and independent are both above 1/2 in all 120 cells;
  - alfak minobs 5 / 10 / 20 only in 20 / 33 / 47 cells.
  - default at distance 1 is 0.745–0.895 by design.
- **stress designs:** default is 0.58–0.99. The weakest is λ 1.6, shift_3000_100, distance 0 (0.58).

## 6. Spearman (overall similarity, not a measure of direction)

| λ | default: distance 0 | 1 | 2 | 3 | 4 | alfak minobs 10: distance 0 | 2 |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.461 | 0.185 | 0.090 | 0.045 | 0.021 | 0.052 | 0.005 |
| 0.4 | 0.466 | 0.266 | 0.192 | 0.165 | 0.141 | 0.042 | 0.011 |
| 0.8 | 0.461 | 0.308 | 0.229 | 0.208 | 0.176 | 0.098 | 0.042 |
| 1.2 | 0.499 | 0.441 | 0.361 | 0.342 | 0.305 | 0.248 | 0.166 |
| 1.6 | 0.465 | 0.615 | 0.556 | 0.527 | 0.492 | 0.385 | 0.308 |

independent is lower than default in most cells, by at most 0.07 (λ 0.4, distance 0). The only exception is λ 1.6, distance 0: independent 0.471, default 0.465.

## 7. Coverage (default, 95% intervals, distances 0–2)

**core designs (averaged first within designs, then across landscapes):**

| λ | Distance 0: calibrated sd / model sd | 1 | 2 |
|---|---|---|---|
| 0.2 | 0.960 / 0.681 | 0.994 / 0.719 | 0.998 / 0.767 |
| 0.4 | 0.957 / 0.705 | 0.985 / 0.614 | 0.978 / 0.548 |
| 0.8 | 0.972 / 0.712 | 0.977 / 0.525 | 0.954 / 0.423 |
| 1.2 | 0.996 / 0.733 | 0.997 / 0.527 | 0.992 / 0.436 |
| 1.6 | 0.991 / 0.770 | 0.997 / 0.578 | 0.998 / 0.542 |

- **Inner and outer misses:**
  - inner misses (the true value falls outside the interval on the side towards 0) are at most 0.6% per cell;
  - outer misses are at most 4.5% (λ 0.8, distance 2);
  - that is, the misses are mainly cases where the estimate underestimates the distance from the reference.
- **Whether intervals contain 0:** the calibrated sd of 1.1.3 contains the term 0.6 × |estimate| (observed karyotypes) or 1.1 × |estimate| (unobserved karyotypes), so the interval always contains 0 by construction. In this round all 24.9 million intervals of core and stress together (13.9 million in core) contain 0.
- **All karyotypes pooled (core):**
  - default: observed karyotypes 0.987 (model sd 0.735), unobserved 0.987 (0.543);
  - independent: 0.982 and 0.962.
- **Worst wavelength by design × distance:**
  - The lowest in the core designs is 0.89: unbalanced_100_3000, distance 2.
  - In the stress designs, when the landscape changes between the two time points, coverage at distance 0 drops clearly:

| Design (distance 0) | Coverage (landscape mean) | Worst wavelength | Inner misses | Outer misses |
|---|---|---|---|---|
| turned_3000_100_360d | 0.881 | 0.748 | 2.6% | 9.3% |
| turned_500_700_360d | 0.895 | 0.818 | 0.9% | 9.6% |
| shift_3000_100_360d | 0.911 | 0.850 | 2.2% | 6.7% |
| shift_500_700_360d | 0.891 | 0.858 | 0.8% | 10.1% |

  In addition, the worst wavelength of short_500_700_30d at distance 2 is 0.876.

## 8. Weak points by design (default)

- **Smooth landscapes, observed karyotypes are the weakest.**
  - λ 1.6, distance 0: very_sparse_50_50 error ratio 0.858, direction 0.623; unbalanced_3000_100 direction 0.651.
  - λ 1.2, short_500_700_30d: error ratio 0.888, direction 0.630.
- **Landscapes that change between the two time points (shift, turned) are the weakest.** At λ 1.6, distance 0:
  - error ratio 0.81–0.90;
  - direction 0.58–0.69;
  - for coverage see Section 7.
- **On rugged landscapes, observed karyotypes are the strongest.** At λ 0.2, 0.4, direction at distance 0 is 0.88–1.00 for all designs.

## 9. Cost

| Method | Median per fit (main / far distance) | CPU hours (main + far distance) |
|---|---|---|
| default | 32 s / 24 s | 137 + 74 |
| independent | 6 s / 6 s | 25 + 15 |
| alfak minobs 5 | 18 s / 21 s | 310 + 280 |
| alfak minobs 10 | 7 s / 8 s | 45 + 54 |
| alfak minobs 20 | 5 s / 5 s | 15 + 18 |

## 10. Where 1.1.3 is good and where it is weak

**Good:**
- Better than the null prediction and alfak in every wavelength × distance cell and every design cell, in both numeric value and direction.
- On rugged landscapes, the direction for observed karyotypes is almost always right (0.96–0.99 at λ 0.2, 0.4).
- The calibrated intervals have coverage of at least 0.95 in the core designs.

**Where the two accuracies part ways:**
- **λ 0.2, distances 3–4:**
  - The error ratio is the lowest in the whole table (0.28), but direction is only 0.57–0.64 and Spearman only 0.02–0.05.
  - There, the karyotypes have almost no ordering to speak of, so the numeric advantage does not come from ordering but from predicting the average level correctly.
  - "Low error ratio" does not mean "accurate direction".

**Weak:**
1. **Observed karyotypes on smooth landscapes (λ ≥ 1.2):** error ratio 0.67–0.75, direction 0.67–0.74, worse even than the unobserved karyotypes at distances 1–2 on the same landscape.
2. **Landscape changes between the two time points:** error ratio up to 0.90, direction down to 0.58, coverage at distance 0 only 0.75 for the worst wavelength.
3. **Intervals always contain 0:** users cannot use the interval to say "this karyotype differs from the reference".
4. **Coverage with the model sd alone is only 0.42–0.77:** coverage relies mainly on the calibration term.

**default vs independent:**
- The same in direction; numerically default is clearly better at λ 0.4–1.2, at the cost of fits 4–5 times slower.
- This round shows no "numeric value for direction" trade-off: default is no worse than independent on either accuracy.

**Not covered by this round:**
- a wrongly set missegregation rate (the truth is always 5e-5);
- dynamics of the clone-replacement kind. In another development study that is not in the repository, coverage on such data fell as low as 0.73;
- real data.

## 11. Runs and files

**Timeline:**
- 2026-11-20 15:42 generation and fits submitted;
- around 22:00 all fits and scoring finished;
- 11-21 01:02 the two summaries finished;
- the f ≤ 0 check was rerun on this round and committed (6b2b5a5, `doc/notes/nonpositive_fitness/round_6ae2a24/`).

**Jobs:** recorded in `round_jobs.tsv` in this round's folder.
- generation 22173660, index 22173662, null prediction 22173663;
- KARST main table 22173665, far distance 22173667;
- alfak 22173669 / 71 / 73 (main table) and 22173675 / 77 / 79 (far distance);
- 64G reruns 22226038 / 40, 22236094 / 96, make-up run for task 602 22236306;
- scoring 22236451 / 53 / 55 / 57;
- summaries 22236111 (`own_summary.json`), 22236113 (`own_coverage.json`).

**Additional computations for this report:**
- 2026-11-21, run in karst_dev on HPC, jobs 22299063 (trial run), 22299075 (20 chunks).
- They only read this round's folder and the e7e32f8 folder; every array was paired with an add-hours monitor.
- Scripts and results are not in the repository:
  - HPC `/share/lab_crd/taoli/Project/karst_dev/results/ownround/`;
  - including `check.rds`, `numeric.rds`, `compare.rds`, and `direction_*.csv`, `coverage_*.csv` (per dataset).

**Definitions:** the same as the report page's `doc/report/build/own_benchmark_summary.R` and `own_coverage.R`. Direction accuracy uses exactly the functions in them; scoring uses `score.R`.
