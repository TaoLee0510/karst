# Cross-validation "compared with what": three targets, round 6ae2a24 (Ver2)

2026-11-21. Ver1 (`xv_target_6ae2a24.md`) had only "the target aligned to alfak". This version adds a like-for-like target that does not depend on alfak, corresponding to KARST 1.1.4's `as_alfak(xval_target = "own")`. All numbers in both versions come from the same round of fits; nothing was refit.

## 1. The question

The cross-validation of both methods uses the same folds (each frequent karyotype together with its one-step neighbours) and reports the same statistic, R2R, but what the held-out prediction is "compared with" differs:

| Target | How it is obtained | Used in |
|---|---|---|
| alfak's target | one bootstrap draw of the karyotype's first-step fitness: mainly from its own counts, with sampling noise | alfak's own cross-validation |
| KARST default target (`"model"`) | KARST's estimate from all the data, from the same smoothing model as the held-out prediction | KARST 1.1.3 and earlier |
| alfak's target lent to KARST (`alfak_target`) | the target values stored in the alfak fit | diagnostic of step 1; depends on the alfak fit |
| KARST's own like-for-like target (`own_target`) | the other latent variables fixed at the fit's mode, the karyotype's fitness prior removed, estimated from its own counts only; the counts are first resampled once by cell, corresponding to alfak's single bootstrap draw | 1.1.4's `xval_target = "own"`; does not depend on alfak |
| Same, without resampling (`own_point`) | same as above, using the observed counts directly | sensitivity analysis (option b of D2) |

The paper's tests first keep only fits with score > 0, and then pick minobs and window length by score, so "compared with what" changes which fits of each method enter the evaluation. This report compares the scores under these targets, and the results of the paper tests.

## 2. Method

**Step 0: karyotype sets.** As in Ver1: for all fits of the three KARST methods, the landscape's karyotypes, the frequent karyotypes and the one-step neighbour sets are all the same as alfak's.

**Step 1: alfak's target** (`scripts/xv_alfak_target.R`). As in Ver1.

**Step 3: KARST's own like-for-like target.**
- **In the package:**
  - New file `R/sparse-joint-own.R`: for each of the four models (two time points, anchored, pairwise, state space) an R version of the "single karyotype's own terms", including its counts and the initial abundance term (balance prior or free prior), and for state space also its own transition term; neighbours appear only through the inflow at the fit's mode.
  - Check: at the fit's mode, the intensity of each karyotype compared with the value reported by TMB, maximum difference 3e-13 (test `tests/testthat/test-sparse-joint-own.R`).
  - `as_alfak()` gains the arguments `xval_target = "model" | "own"`, `xval_resample`, `xval_seed`. Default `"model"`, with output unchanged byte for byte. Only karyotypes with at least 1 cell are tested (option a of D1).
- **Benchmark:** `scripts/xv_own_target.R`.
  - The held-out predictions do not depend on "compared with what", so the stored ones are used directly; only the scores are recomputed against the new target.
  - The target is converted to the landscape's scale: `(own net growth rate − reference value + level) / efflux correction`. The level is read back from the stored `landscape.Rds`; each fit has only one value, deviation ≤ 2e-16.
- **Check:** 213 fits were recomputed with the full `as_alfak(xval_target = "own")`, including refitting each fold; the scores are exactly the same as the post-hoc computation (difference 0), landscape differences ≤ 9e-16.
- **Coverage:** 5546 of the 5617 fit folders have scores; the rest are fits that failed.
- **Rerun of the paper tests:** `XV_SOURCE=own_target`, done as in step 2 of Ver1, results in `analysis_xv_own_target/`. `own_point` computes only the scores and Fig. 2d, without rerunning all tests.
- **Comparison tables:** `scripts/xv_target_report.R`, in `analysis_xv_targets/`. The script checks that the rows for alfak and for the paper's fits are exactly the same across all runs. This time they are all the same.

## 3. The scores themselves

KARST anchored, fits with status ok; inputs excluded by the analysis are not counted. In parentheses, the fraction with score > 0.

| Data | Design | Default | alfak target | Own target (resampled) | Own target (not resampled) | alfak itself |
|---|---|---|---|---|---|---|
| ABM | Two time points | 0.51 (70%) | 0.001 (51%) | 0.004 (57%) | 0.013 (51%) | −0.03 (41%) |
| ABM | All time points | 0.57 (73%) | −0.003 (41%) | 0 (50%) | 0.13 (58%) | −0.01 (44%) |
| Missegregation rate | Two time points | 0.86 (93%) | 0.055 (80%) | 0.003 (58%) | 0.23 (83%) | 0 (50%) |
| Salehi | Two time points | 0.75 (72%) | 0.016 (61%) | 0.003 (54%) | 0.098 (72%) | −0.15 (30%) |
| Salehi | All time points | 0.80 (77%) | 0.074 (70%) | −0.008 (36%) | 0.25 (74%) | 0.04 (56%) |

- **Both like-for-like targets bring KARST's scores down to around 0, on the same scale as alfak.**
- **The target that does not depend on alfak (resampled) is the most similar to alfak's scores:** the median is about 0.
- **Without resampling the scores are higher:** the sampling noise of the target is removed.

## 4. How close the targets are to the true fitness (simulations)

For each fit, on the karyotypes that KARST tested and that have cells, compute the Spearman between each target and the true fitness, then take the median across fits (`truth_targets.tsv`, ABM two time points):

| λ | Own target (resampled) | Own target (not resampled) | KARST default target | alfak's target |
|---|---|---|---|---|
| 0.2 | 0.41 | 0.40 | 0.59 | 0.17 |
| 0.4 | 0.56 | 0.71 | 0.76 | 0.38 |
| 0.8 | 0.60 | 0.85 | 0.89 | 0.48 |
| 1.6 | 0.59 | 0.83 | 0.88 | 0.31 |

- **KARST's estimate from its own counts only is closer to the true fitness than alfak's target**, also after resampling.
- **KARST's default target is closest to the truth,** because it uses all the data and the landscape prior. So the default scores are high not because the target is poor, but because the target and the prediction come from the same model.

## 5. Can the scores reflect the true accuracy of the landscape (simulations)

Within each window × wavelength cell, the Spearman between the score and the landscape's true Pearson r, with the median taken across cells (`truth_scores.tsv`, KARST anchored):

| Data | Design | Default | alfak target | Own target (resampled) | Own target (not resampled) | alfak itself |
|---|---|---|---|---|---|---|
| ABM | Two time points | 0.30 | 0.36 | 0.08 | 0.18 | 0.31 |
| ABM | All time points | 0.36 | 0.35 | 0.03 | 0.31 | 0.39 |
| Missegregation rate | Two time points | 0.43 | 0.44 | 0.07 | 0.26 | 0.43 |

- **Although the once-resampled target is closer to the truth, the scores computed with it hardly reflect the accuracy of the landscape** (0.03–0.08). The noise of a single resample determines R2R. This is the same flaw as in alfak's own cross-validation: it too compares with one noisy draw.
- **Without resampling it is somewhat better** (0.18–0.31), but no better than the default score.
- **On Salehi, none of the scores predicts how good the forecast is:** the Spearman between the score and "the fraction of forecasts beating the parent on daughters" lies between −0.06 and 0.24 (`forecast_salehi.tsv`).

## 6. Fig. 2d: ranking of the magnitude of change (Salehi)

Spearman (actual change, predicted change), evaluation day, Wasserstein (`fig2d_swap.tsv`):

| KARST score | Two time points: KARST / alfak under the paper's rule | On the rows KARST picked: KARST / alfak | All time points: KARST / alfak under the paper's rule |
|---|---|---|---|
| Default | 0.218 / 0.447 | 0.218 / 0.344 | 0.282 / 0.695 |
| alfak target | 0.372 / 0.447 | 0.372 / 0.376 | 0.464 / 0.695 |
| Own target (resampled) | 0.298 / 0.447 | 0.298 / 0.275 | 0.289 / 0.695 |
| Own target (not resampled) | 0.423 / 0.447 | 0.423 / 0.487 | 0.271 / 0.695 |

Without selection, using all rows shared by both: two time points KARST 0.295, alfak 0.322 (difference −0.03, interval −0.19 to 0.15); all time points 0.300 vs 0.556. These two rows do not depend on the score.

- **Two time points:** switching to either like-for-like target narrows the gap under the paper's rule (from 0.229 down to 0.02–0.15). On the rows KARST itself picked, KARST is even slightly higher with the resampled target (+0.02).
- **All time points:** most of the gap does not change with the target, because anchored already ranks lower on all rows, independent of selection.

## 7. Other tests that select by score

KARST two time points (all time points); the numbers for alfak and the paper are exactly the same across runs.

| Test | Default | alfak target | Own target | alfak |
|---|---|---|---|---|
| Fig. 2f beats parent (Wasserstein) | 0.48 (0.53) | 0.56 (0.54) | 0.54 (0.54) | 0.25 (0.30) |
| Fig. 2f beats parent (Euclidean) | 0.34 (0.47) | 0.48 (0.52) | 0.44 (0.50) | 0.25 (0.24) |
| Fig. 2c angle (smaller is better) | 72° (66°) | 81° (78°) | 56° (62°) | 83° (80°) |
| Fig. 3 new karyotypes: fraction with significant fitness | 0.51 (0.46) | 0.50 (0.44) | 0.56 (0.50) | 0.23 (0.15) |
| Fig. 2b median cross-validation score (> 3 passages) | 0.89 (0.92) | 0.14 (0.16) | 0.01 (−0.01) | 0.05 (0.47) |
| Fig. 4–6 agree with the paper (16 in total) | 10 (7) | 9 (10) | 11 (10) | 12 (13) |

- **The other two groups in Fig. 4–6:** pairwise and state space are 8/9/8 and 8/9/10 respectively under the three scores.
- **ABM:** after selection, the landscape accuracy and Fig. 1g hardly change under the three scores, and KARST still clearly leads. For example, for two time points, window 4, the fraction beating "no change" at passage 10: default 0.83, alfak target 0.87, own target 0.82, alfak 0.51.

## 8. Conclusions

1. **"Compared with what" does affect KARST's ranking in the paper tests.** On two time points, after switching to a like-for-like target, the Fig. 2d gap narrows clearly, and most results of Fig. 2c, 2f, 3 and Fig. 4–6 improve. On all time points the Fig. 2d gap is mainly unrelated to selection.
2. **KARST's estimate from its own counts only is closer to the true fitness than alfak's target**, but after resampling once, the scores computed with it hardly reflect landscape accuracy. This shows that cross-validation with "one noisy draw" as the target cannot by itself select good fits, and the same holds for alfak.
3. **Recommendation: keep the 1.1.4 default at `"model"` (option a of D3).** How well it reflects accuracy is no worse than any of the others; `"own"` serves as the option equivalent to alfak, for use when a comparison on the same scale as alfak is needed.

## 9. Limitations

- **Different karyotype ranges:** the target from own counts only tests only karyotypes with cells; alfak also tests neighbours without cells, whose values come mainly from the prior. The sets of karyotypes tested by the two are not exactly the same.
- **A single resample is random:** the seed is fixed at 1; changing the seed changes the score of an individual fit a lot.
- **Outflow and the neighbours' balance prior are not counted:** the single-karyotype estimate does not count the part that flows to its neighbours, nor its role in the neighbours' balance prior. This is the definition of "own counts only", not an omission.

## 10. Files and reproduction

- **Package:**
  - `R/sparse-joint-own.R`;
  - `R/alfak-compat.R` (new arguments of `as_alfak()`, `alfak_compat()`);
  - `man/as_alfak.Rd`, `man/alfak_compat.Rd`;
  - `tests/testthat/test-sparse-joint-own.R`.

  All tests pass: 18 files, 113 tests, 0 failures.
- **Benchmark:**
  - `scripts/xv_own_target.R`, `scripts/xv_target_report.R`;
  - `R/xv_source.R` (`own_target`, `own_point`);
  - `hpc/env.sh` (`KARST_LIB` can be set from the environment, and variables such as `OWN_CHUNK_DIR`, `CHECK_ALL`, `ALIGN` are passed into the image).
- **Results (`results_6ae2a24/` on HPC, with a local copy):**
  - `tables/xv_own_target.Rds`;
  - `analysis_xv_own_target/`, `figures_xv_own_target/`;
  - `analysis_xv_targets/` (comparison tables of the three targets).
- **Jobs:** the code is not committed yet; it was run from the staged copy `results_6ae2a24/xv_target_code/`, with the library in `lib/karst_wip/`:
  - full run 22299789 (200 chunks), merge 22299791;
  - paper-test reruns 22300293, 22300295, 22300297, 22300299.

  All COMPLETED.

  The first full submission (22299566) wrote to a wrong path because an environment variable was passed as an empty string; it was cancelled and resubmitted.
