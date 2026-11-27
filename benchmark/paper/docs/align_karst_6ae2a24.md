# Paper tests with alfak aligned to KARST: version B, round 6ae2a24

2026-11-23. All numbers come from the fits of round 6ae2a24; B only changes the karyotype table and the cross-validation target, and nothing was refitted.

## Conclusions

1. **Prediction: KARST is clearly better than alfak in both versions.**
   - ABM: the prediction results of A and B are almost the same: the mean difference over 10 passages is < 0.01, and at a single passage at most 0.017.
   - Salehi paired, wasserstein: A 0.51 vs 0.20, B 0.54 vs 0.17.
2. **Landscape accuracy (simulations): on the whole of E, both methods are lower than on alfak's table.** KARST still leads; alfak drops more.
3. **In B, alfak's landscape is mostly its Kriging extrapolation:** ABM median 80%, Salehi 98%.
4. **In B, the paper's selection rule keeps a different set of fits.**
   - With the like-for-like target, KARST has fewer fits with score > 0: on Salehi 42–45 trajectories per method, against 48–50 in A.
   - In the results that depend on this rule (Fig. 1g, 2c, 2d, 2f, 3, 5), part of the difference between A and B comes from the different fits that were selected.
5. **Number of statements in Fig. 4–6 that agree with the paper's conclusions:** alfak rises from 12–13 to 13–15, KARST falls from 7–10 to 6–7. The reference is the paper's own alfak fit, so this only shows closeness to the paper, not greater accuracy.
6. **In both versions, the steady states of Fig. 6 are determined by extrapolated fitness at the edge of the table:** the dominant karyotype at the lowest missegregation rate is almost always the edge karyotype with the highest fitness on the table. B's table is larger, which amplifies this.

## 1. Question

The paper tests of this round (version A) are **KARST aligned to alfak**:
- For each minobs, both methods report landscapes on alfak's karyotypes. alfak's karyotypes = the frequent karyotypes with at least minobs cells + their one-step neighbours + the neighbours' neighbours.
- Prediction, scoring and Fig. 3–6 all use these karyotypes.
- So KARST's results also carry minobs: one landscape per minobs, and the paper's rule then picks one among the minobs by cross-validation score.

Version B is the reverse, **alfak aligned to KARST**:
- The test methods are unchanged; only the karyotypes and the minobs-dependent sets are replaced by KARST's.
- All of alfak's minobs fits are kept.
- Karyotypes beyond alfak's own landscape are filled in with alfak's own Kriging.
- Run with `ALIGN=karst`; the code is in `R/align.R`; the results are in `results_6ae2a24/align_karst/`, with the same layout as the round.
- The fits themselves are unchanged: both methods use this round's fits, and nothing was refitted.

## 2. Method

**Extended karyotype table E.**
- The same approach as in `dev/realdata/extend`.
- E = all karyotypes observed at the time points used in the fit + their neighbours within two steps. One step = one more or one fewer copy of a single chromosome, with the copy number not below 1.
- E does not depend on minobs, and contains alfak's karyotype table for every minobs (checked in every task).
- Each karyotype records its distance class to the observed karyotypes: 0, 1, 2.

**alfak on E (one per minobs).**
- On alfak's own karyotype table, the mean of its landscape is used.
- The other karyotypes are predicted with its stable Kriging model: `krig_stable_mean` of alfakR `fitKrig()`, i.e. `fields::Krig` fitted to the bootstrap means of the frequent karyotypes and their neighbours.
- The model is rebuilt from each fit's `bootstrap_res.Rds`.

**KARST on E (one per method, shared by all minobs folders).**
- Fitness on E is computed with `predict()`.
- The scale is taken from its minobs 10 landscape: the level of `as_alfak()` (read back from the stored landscape, one value per fit) and the efflux correction of each karyotype.

**The minobs-dependent sets are replaced by KARST's, the same for both methods:**
- frequent karyotypes `fq` = the major karyotypes of KARST's default fit;
- `nn` = their one-step neighbours;
- landscape accuracy on simulated data is reported by distance class 0, 1, 2. The three columns of the same names in version A are alfak's frequent karyotypes, one-step neighbours and two-step neighbours, so the columns `rfq`, `rnn`, `rd2` cannot be compared directly between the two versions; only the whole-table `r` can.

**Cross-validation score (the paper's rule picks fits by it).**
- alfak keeps its own.
- KARST uses a like-for-like target that does not depend on alfak: the folds are centred on KARST's major karyotypes, and the targets come from each karyotype's own counts, resampled once by cell, seed 1. That is `as_alfak(xval_target = "own")`; see `xv_target_6ae2a24_Ver2.md`.

**Prediction and scoring.**
- The prediction runs alfak's ABM on E in the same way as `predict_evo()`: the initial allocation, time step, recording interval and seed are all the same.
- The results store only the karyotypes with a frequency greater than 0 at each recorded time.
- Before scoring, they are expanded to all karyotypes that appeared in the prediction or the observation, which is equivalent to scoring on the whole of E.

**Fig. 6.**
- When building its matrix, alfakR's `find_steady_state(Nmax = 1)` loops over all pairs of karyotypes and reserves n² elements; on E (up to 1.6 million karyotypes) there is not enough memory.
- `align_find_steady_state()` (`R/align.R`) is used instead. With Nmax = 1 the matrix only connects each karyotype to itself and its one-step neighbours, so the neighbours are found by table lookup; each element is the same product of alfakR's per-chromosome probabilities, in the same order; the eigenvector step is copied from alfakR unchanged.
- Only karyotypes with a share greater than 0 at some missegregation rate are kept. The Fig. 4–6 summaries use only the dominant karyotype at each missegregation rate, so no information is lost.

## 3. Checks

- **Sparse prediction = `predict_evo()`** (`scripts/align_check.R`, log `xv_target_code/logs/align_check_22300070.out`): 5 ABM and 3 Salehi inputs, 1962–17314 karyotypes; the maximum difference of the frequencies at every recorded time is 0.
- **Kriging rebuild = the model stored by alfak** (same log): 5 Salehi trajectories, 3909–8171 points, maximum difference 0. The Pearson r between Kriging and the mean on alfak's own karyotype table is 0.97–1.00.
- **Fig. 6 steady states = alfakR's** (`scripts/align_check_steady.R`, `align_karst/align_checks_steady.Rds`):
  - this round's 6 Salehi landscapes (2757–16065 karyotypes) × 8 missegregation rates: in 48/48 the matrices and steady states are `identical()`, and the dominant karyotypes are the same.
  - on B's longest table (1.64 million karyotypes): the matrix has 10.4 million non-zero elements, about 30 seconds per missegregation rate, peak 2.6 GB.
- **The karyotype sets of the two methods are the same** (`scripts/check_ksets.R`): for every fit folder in which both methods have a landscape, the karyotype table E, `fq` and `nn` were compared. ABM and missegregation rate 11,747, Salehi 2,366, all identical.
- **The tasks are complete**: ABM and missegregation rate 11,925 and Salehi 2,676 fit folders, none with missing output.

## 4. Size and time

`scripts/align_summary.R`, recorded in `align_karst/align_tasks.Rds`. One task = one landscape of one input (one minobs of alfak, or one method of KARST); 13,688 tasks in total, all completed.

**Length of E (number of karyotypes, minimum / median / maximum).**

| Data | Design | Observed karyotypes (median) | E |
|---|---|---|---|
| ABM | Two time points | 105 (all ABM) | 9805 / 49,405 / 113,043 (window 2) |
| ABM | All time points | | 12,567 / 62,924 / 148,341 (window 4); 13,182 / 77,869 / 189,256 (window 8) |
| Missegregation rate | Two time points | 629 (all missegregation rate) | 23,196 – 28,754 / about 250,000 / about 1.1 million |
| Missegregation rate | All time points | | 36,461 / 362,968 / 1,867,749 (window 4); 54,993 / 487,100 / 3,122,551 (window 8) |
| Salehi | Two time points | 463 (all Salehi) | 18,904 / 256,121 / 861,895 |
| Salehi | All time points | | 37,502 / 552,752 / 1,643,257 |

**Composition of E:** the vast majority are two-step neighbours. 8 landscapes each, drawn at random from Salehi and from ABM (24,000–790,000 karyotypes): observed karyotypes make up 0.1–0.2%, one-step neighbours 5–6%, two-step neighbours 94–95%.

**How much of alfak's landscape relies on Kriging (fraction of E filled in by Kriging, 10% / median / 90%).**

| Data | Fraction | Pearson r between Kriging and the mean on alfak's own table |
|---|---|---|
| ABM | 0.67 / 0.80 / 0.88 | 0.88 / 0.98 / 1.00 (minimum 0.29) |
| Missegregation rate | 0.80 / 0.91 / 0.99 | 0.91 / 0.97 / 0.99 (minimum 0.66) |
| Salehi | 0.91 / 0.98 / 1.00 | 0.92 / 0.99 / 1.00 (minimum 0.30) |

That is, in B alfak's landscape is mostly its Kriging extrapolation; on its own table, Kriging is mostly very close to its mean.

**KARST's landscape on its own table is the same as in this round**: maximum difference ≤ 7e-16 for every fit.

**Time.**
- All tasks together took 3434 CPU hours: ABM 1095, missegregation rate 98, Salehi 2241.
- Most of it was spent on prediction. The slowest were the two trajectories with the longest E in Salehi all time points (1.27 million and 1.36 million karyotypes): 14–18 hours per prediction.
- Peak memory at most 9.5 GB (Salehi).

## 5. Results: simulated data (ABM and missegregation rate)

The summary tables are in `analysis_align_compare/` (`scripts/align_report.R`); each table has a column `version` (A, B).

**Prediction (paired, KARST two time points vs alfak two time points, mean over 10 passages): the two versions are almost the same.**

| Window | Metric | Fraction beating the baseline, KARST A / B | alfak A / B |
|---|---|---|---|
| 2 | euclidean | 0.781 / 0.772 | 0.433 / 0.433 |
| 4 | euclidean | 0.901 / 0.897 | 0.530 / 0.530 |
| 8 | euclidean | 0.900 / 0.899 | 0.485 / 0.479 |

The other three metrics and the all-time-points results are the same: averaged over 10 passages, the differences are all within 0.01; at a single passage at most 0.017 (KARST) and 0.011 (alfak). E is much longer than alfak's table, yet the prediction results hardly change.

**Fig. 1g (picked by cross-validation, score > 0).**
- alfak: averaged over 10 passages, the differences between A and B are all within 0.01; at a single passage at most 0.04.
- KARST: averaged over 10 passages, a few cells differ by up to ±0.06; at a single passage at most 0.12. The reason is that in B KARST's score uses the like-for-like target: fewer inputs have a score > 0, 34–92 per cell, against 52–100 in A. So in B the paper's rule keeps a different set of inputs. This is not a change in accuracy, and should not be read as KARST getting better or worse in B.

**Landscape accuracy (Pearson r, paired, KARST two time points vs alfak two time points).**

| λ | Window | KARST A / B | alfak A / B | Fraction where KARST is better, A / B |
|---|---|---|---|---|
| 0.2 | 2 | 0.201 / 0.134 | 0.084 / 0.055 | 0.887 / 0.891 |
| 0.4 | 4 | 0.508 / 0.437 | 0.384 / 0.299 | 0.946 / 0.977 |
| 0.8 | 8 | 0.855 / 0.827 | 0.724 / 0.670 | 0.997 / 1.000 |
| 1.6 | 8 | 0.869 / 0.857 | 0.608 / 0.553 | 0.993 / 1.000 |

- On E, both methods have a lower r.
- alfak drops more, especially at large λ; outside its own table it relies entirely on Kriging.
- The fraction where KARST is better is similar in the two versions: A 0.84–1.00, B 0.79–1.00, with some rising and some falling (λ 0.2 window 4: 0.964 → 0.924; λ 0.4 window 2: 0.835 → 0.791).
- By distance class, which only B has (columns `rfq`, `rnn`, `rd2`, paired, KARST two time points vs alfak two time points, 12 cells = 4 λ × 3 windows):

  | Distance class | KARST median r | alfak median r | Fraction where KARST is better |
  |---|---|---|---|
  | 0 observed karyotypes | 0.54–0.88 | 0.14–0.82 | 0.79–0.98 |
  | 1 one-step neighbours | 0.14–0.87 | 0.07–0.76 | 0.75–0.99 |
  | 2 two-step neighbours | 0.06–0.86 | 0.05–0.67 | 0.61–1.00 |

  - At λ = 0.2, the farther from the observed karyotypes, the closer both methods get to 0 (two-step neighbours: KARST 0.06–0.08, alfak 0.05).
  - At λ ≥ 0.8, KARST still has 0.72–0.86 on the two-step neighbours.

**Supp. Fig. 3 (missegregation-rate data, median Pearson r): both methods are lower in B, and KARST is still above alfak.**
- KARST two time points minus alfak two time points, 15 cells (3 windows × 5 missegregation rates): in A all positive, 0.01–0.35; in B 14 positive, −0.05–0.35, with the only negative value at window 8, missegregation rate 8e-4.
- All time points: A has 3 negative cells, B has 2.

## 6. Results: Salehi

The paper's own fits ("paper (alfak main@4ad82b0)") use the results stored by the paper; they are not aligned in B, are the same in both versions, and serve only as a reference. Neither version refits, so the number of failed fits is exactly the same.

**The fits selected by the paper's rule changed.**
- **KARST:** after the cross-validation switched to the like-for-like target, fewer trajectories have a score > 0: two time points 50 → 44, all time points 49 → 42, pairwise 50 → 43, state space 48 → 45.
- **alfak:** the scores did not change, but all time points gained 2 (46 → 48).
  - These two fits are SA000_X4_l_4_d1_3_d2_5 and SA609TIIDT_X4_l_4_d1_1_d2_0, both with minobs 5.
  - In A their predictions failed, because "none of the cells of the last passage are in the landscape": alfak's table does not contain these karyotypes.
  - E contains all observed karyotypes, so in B they can be predicted.
- In the paired comparison without picking, one more fit is added for the same reason: SA535_CISPLATIN_CombinedT_X7_l_3_d1_2_d2_2, minobs 20, two time points.

**Paired prediction (no picking, all fits): KARST is clearly better than alfak in both versions.**

Fraction beating the baseline, KARST / alfak:

| Comparison | Metric | A | B |
|---|---|---|---|
| KARST two time points vs alfak two time points | wasserstein | 0.512 / 0.204 | 0.539 / 0.166 |
| | euclidean | 0.440 / 0.241 | 0.468 / 0.203 |
| KARST all time points vs alfak all time points | wasserstein | 0.498 / 0.238 | 0.532 / 0.194 |
| | euclidean | 0.457 / 0.238 | 0.387 / 0.195 |

- McNemar p ≤ 1e-9 for every row.
- alfak's fraction beating the baseline is lower in B throughout.
- KARST has some rising and some falling. wasserstein rises throughout. For all time points, euclidean falls for all three methods: anchored 0.457 → 0.387, pairwise 0.471 → 0.439, state space 0.458 → 0.392.

**Fig. 2f (paper rule, summarised over the evaluation passages, fraction beating the baseline).**

| Fit | wasserstein A / B | euclidean A / B |
|---|---|---|
| KARST two time points | 0.48 / 0.59 | 0.34 / 0.48 |
| alfak two time points | 0.25 / 0.20 | 0.25 / 0.30 |
| KARST all time points | 0.53 / 0.55 | 0.47 / 0.43 |
| alfak all time points | 0.30 / 0.25 | 0.24 / 0.27 |
| KARST pairwise | 0.62 / 0.60 | 0.48 / 0.51 |
| KARST state space | 0.63 / 0.53 | 0.50 / 0.40 |
| The paper's own fits (reference) | 0.38 / 0.38 | 0.36 / 0.36 |

On these two metrics, KARST is above alfak in both versions. For cosine, in A alfak two time points (0.32) is above KARST (0.28), and in B it is the reverse (0.18 vs 0.43). KARST has fewer trajectories in B, see above.

**Fig. 2c (angle between the predicted and the actual change; a smaller median is better): almost all smaller in B.**

| Fit | A | B |
|---|---|---|
| KARST two time points | 72.1° | 49.3° |
| alfak two time points | 83.2° | 64.1° |
| KARST all time points | 66.0° | 56.6° |
| alfak all time points | 79.9° | 65.6° |
| KARST pairwise | 64.7° | 66.1° |
| KARST state space | 73.7° | 63.0° |

- All became smaller except pairwise.
- p < 0.05 for all fits (in A, alfak two time points 0.010, the rest ≤ 0.001).
- In both versions, the angle of KARST anchored is smaller than that of alfak.

**Fig. 2d (Spearman correlation between the predicted and the actual Wasserstein distance from the parental passage).**

| Fit | A | B |
|---|---|---|
| alfak two time points | 0.45 | 0.54 |
| KARST two time points | 0.22 | 0.35 |
| alfak all time points | 0.70 | 0.65 |
| KARST all time points | 0.28 | 0.30 |
| KARST pairwise | 0.57 | 0.47 |
| KARST state space | 0.56 | 0.51 |

- In both versions alfak is above KARST anchored; the gap narrows slightly in B: two time points 0.23 → 0.19, all time points 0.41 → 0.36.
- Pearson changes in mixed directions: alfak two time points 0.76 → 0.37, KARST two time points 0.47 → 0.66.

**Fig. 3e–f (newly appearing karyotypes).**
- Candidate karyotypes for the regression = karyotypes on the landscape that are not in the parental passage.
- In B all methods use the same table E, so the regressions of the two methods can be compared directly.

| Fit | Fraction with significant fitness, paper rule A / B | All fits A / B |
|---|---|---|
| KARST two time points | 0.51 / 0.50 | 0.57 / 0.56 |
| alfak two time points | 0.23 / 0.34 | 0.21 / 0.23 |
| KARST all time points | 0.46 / 0.48 | 0.53 / 0.57 |
| alfak all time points | 0.15 / 0.23 | 0.23 / 0.27 |

- Paired (all fits, KARST two time points vs alfak two time points): B 0.562 vs 0.233, McNemar p = 5.8e-26; A 0.567 vs 0.208.
- **The distance term becomes stronger:** in B, the fraction in which "the proportion of cells in the parental population one step away from it" (d1) is significant rises from about 0.6 to about 0.9, the same for all methods.
  - The reason is consistent with the definition of the regression: most of the extra karyotypes in E are two steps away from the observed karyotypes, have d1 = 0, and never appeared.
  - So the fraction in which fitness is the most significant variable falls: KARST two time points 0.40 → 0.19, alfak 0.17 → 0.06.

## 7. Results: Fig. 4–6

- The reference is the same as in A: the conclusions obtained from the paper's own fits with the same code (Fig. 4 uses the bootstrap stored by the paper).
- The paper's own fits are not recomputed in B; they are copied over unchanged.
- "Agree" only means that the conclusion is the same as the paper's, not which method is more accurate: the paper's conclusions themselves come from alfak's landscapes.

**Number of statements that agree with the paper's conclusions (statements tested, 16 in total).**

| Fit | A | B |
|---|---|---|
| alfak two time points | 12 | 13 |
| alfak all time points | 13 | 15 |
| KARST two time points | 10 | 7 |
| KARST all time points | 7 | 6 |
| KARST pairwise | 8 | 7 |
| KARST state space | 8 | 7 |
| The paper's own fits | 16 | 16 |

In B alfak is closer to the paper and KARST farther from it.

**Statements for KARST two time points that changed from A to B.**
- **Became disagreeing:**
  - Fig. 4d cisplatin group variance: A significant, 0.36; B not significant, 0.25.
  - Fig. 4f main effect of distance: A not significant, B significantly negative, −0.056.
  - Fig. 5b the three WGD terms: A agrees; in B all three are significant but in the direction opposite to the paper, for example dk × WGD+ is −0.042, against +0.016 in the paper.
- **Became agreeing:** the two terms of the Fig. 4f parallel pairs.
- **Disagreeing in both versions:** the two PDX terms of Fig. 4d. In both A and B, KARST gives a significant PDX_v opposite to the paper: A −1.15, B −1.18. This is the same as the alfak2-era finding, and is not given a biological interpretation.
- **Fewer fits in Fig. 5:** KARST has only 10–14 in B, against 16–28 in A. Fig. 5 uses only fits with a score > 0, and B's like-for-like target leaves fewer fits with a score > 0.
- **Fig. 5c (KS permutation test):** in both versions, all methods disagree with the paper. The paper has p = 0.03; none of the others is significant.

**Fig. 6 (descriptive, not counted in the numbers above).**

| Fit | Fraction where the dominant karyotype changes with the missegregation rate, A / B | Distance of the dominant karyotype at the highest missegregation rate from that at the lowest (median), A / B |
|---|---|---|
| alfak two time points | 0.77 / 0.97 | 2 / 9 |
| alfak all time points | 0.75 / 0.92 | 1 / 12 |
| KARST two time points | 0.68 / 0.47 | 1 / 0 |
| KARST all time points | 0.82 / 0.60 | 2 / 3 |
| KARST pairwise | 0.80 / 0.45 | 2 / 0 |
| KARST state space | 0.84 / 0.64 | 2 / 3 |
| The paper's own fits | 0.82 | 2 |

- Fraction where the dominant karyotype is the same as for the paper's own fits: A 0.03–0.38, B ≤ 0.06.
- **Where the dominant karyotypes come from** (`scripts/align_check_fig6.R`):
  - **B, lowest missegregation rate:** in almost every landscape, the dominant karyotype is the karyotype with the highest fitness on the whole table, located at the edge of the table. For alfak, 31 of 35 are Kriging values on two-step neighbours; for KARST, 31 of 32 are two-step neighbours.
  - **A is the same:** the dominant karyotype at the lowest missegregation rate is also the one with the highest fitness on the table, and almost never among the frequent karyotypes and their one-step neighbours (alfak 34 of 35, KARST 36 of 40).
  - So the steady states of Fig. 6 are determined mainly by extrapolated fitness at the edge of the table. B's table is larger and its edge farther out, so alfak's dominant karyotype at high missegregation rates is farther away (median 9–12 steps).
- Fig. 6 should therefore not be read, in either version, as a prediction of biological steady states. B amplifies this.

## 8. Files

- Code: `R/align.R`, `scripts/align_task.R`, `align_check.R`, `align_check_steady.R`, `align_check_fig6.R`, `align_summary.R`, `align_report.R`, `check_ksets.R` (parallel, merging rows from other sources), `fig6_steady_states.R` (uses `align_find_steady_state()` when `ALIGN=karst`).
- Results: `results_6ae2a24/align_karst/` (fit folders, scores, `analysis/`, `tables/ksets.Rds`); comparison tables in `results_6ae2a24/analysis_align_compare/`.
