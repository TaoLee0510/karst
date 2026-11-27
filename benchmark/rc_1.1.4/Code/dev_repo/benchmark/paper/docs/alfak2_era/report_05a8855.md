# Comparison report of alfak2 05a8855 (three multi-time-point methods) and alfakR on the ALFA-K paper tests

Date: 2026-11-15. alfak2 version 05a8855 (`hpc/alfak2.lock`), alfakR version alfakR_dev2 f93b2d2.

> Terminology: the input columns are cells taken at different times; this report calls them "time points" ("time points" in the figures), not "samples".
> Samples in the statistical sense (bootstrap etc.) and the Salehi source lines (SA609 etc.) are not covered by this.
> The figures and tables have been redrawn accordingly (without refitting); the versions before the change are in `results_05a8855/pre_time_points/`, with identical numbers except for the labels and the auto Fig. 4–6.

- Results are in `results_05a8855/`; the results and reports of earlier rounds are where they were:
  - 8327282: `results/`, `docs/report.md`
  - ff13cf4: `results_ff13cf4/`, `docs/report_ff13cf4.md`
  - 8512694: `results_8512694/`, `docs/report_8512694.md`
  - c7f809c: `results_c7f809c/`, only a cost test (54 tasks), superseded by 05a8855
- Tables are in `results_05a8855/analysis/{abm,salehi,summary,fig456}/`, figures in `results_05a8855/figures/`, written as `figures/…` in the text.
  The summary figure is `figures/summary_comparison.png`; the Fig. 4–6 figures are `figures/fig4_effects.png`, `fig5_wgd.png`, `fig6_dominance.png`.

## 0. What this round tests

| commit | Change | Effect on this benchmark |
|---|---|---|
| c7f809c | `fit_sparse_joint(series = …)`: besides the default `"anchored"`, adds `"pairwise"` (pair-by-pair fits of adjacent time points) and `"state_space"` (fitness can drift between time points, `diagnostics$drift`) | Two more fits for each "all time points" input; pairwise requires every pair of adjacent time points to have coverage, so 42 Salehi sweep inputs cannot be fitted |
| 05a8855 | When pairwise / state_space meet a pair without coverage, time points are removed until every pair has coverage (`sparse_joint_pair_left_out()`); when only the first and last remain, it falls back to the two-time-point fit | Each of the 42 inputs loses 1 time point, 4 of them fall back to the two-time-point fit; anchored, two-time-point and fully covered fits are bitwise identical to c7f809c |

- The alfakR fits are the same as in round 8327282, reused through symbolic links.
- All 5617 alfak2 tasks were rerun:
  - two time points (ends) and anchored all time points: 3107, identical to 8512694 (checked in section 6);
  - pairwise and state_space all time points: 1255 each.

## 1. Conclusions

1. **All three alfak2 all-time-point fits beat the alfakR all-time-point fit on almost every accuracy test** (section 3),
   consistent with the conclusion of 8512694. The only exception is pairwise on the 8-passage window of the missegregation-rate experiment (see below).
2. **Within alfak2, the three methods each have their strengths; none wins across the board:**
   - **state_space has the most accurate landscape**: ABM paired median r 0.66 / 0.72 (4- / 8-passage windows), anchored 0.66 / 0.71;
     state_space is higher on 65% of the inputs. The same on the missegregation-rate experiment (0.77 / 0.83 vs 0.76 / 0.82, 67–71%).
     But the ABM forecast is slightly worse: beating the baseline at passage 10, 0.84 / 0.89 vs 0.86 / 0.91 (Euclidean, p ≤ 5e-4).
     On Salehi the Wasserstein forecast is also worse: 0.49 vs 0.54 (p = 0.007).
   - **pairwise is clearly the worst on ABM**: r 0.57 / 0.63, lower than anchored on 68–70% of the inputs;
     forecast 0.76 / 0.78. On the 8-passage window of the missegregation-rate experiment it is even below alfakR all time points (0.74 vs 0.80).
   - **pairwise is the best on Salehi new karyotypes**: fraction where fitness predicts the appearance of new karyotypes 0.66 vs anchored 0.52 (p = 4e-9, 351 pairs);
     the forecast is not significantly different from anchored.
3. **Sweep subset** (Salehi trajectories where anchored removed time points, 48 forecast pairs, 92 new-karyotype regressions):
   - Forecast beats parent: Wasserstein pairwise 0.66, anchored 0.64, state_space 0.63;
     Euclidean 0.48, 0.44, 0.39. None of the differences is significant (p ≥ 0.29).
   - New karyotypes: pairwise 0.58 vs anchored 0.47 (p = 0.03), state_space 0.45.
   - alfakR on this subset still has only 0.11.
4. **Biological conclusions of Fig. 4–6** (section 5):
   - The paper's fits, recomputed with this repository's code, reproduce all 19 conclusions; alfakR dev2 all time points reproduces 15, two time points 14.
   - The alfak2 arms reproduce only 9–11 (the auto arm, not adopted, 13). The largest disagreement is **Fig. 4d: alfak2 gives a conclusion in the opposite direction to the paper**:
     the paper says that the fitness effect of copy-number changes is larger in PDX than in vitro (PDX_v +1.23), while all alfak2 arms (including auto) are negative (−0.38 to −1.40),
     because alfak2's |Δf| in vitro is 4–8 times the paper's, while in PDX they are similar.
   - Fig. 5c (WGD+ and WGD− have different Δf distributions): only the paper's fits are significant (p = 0.03); alfak2 arms p 0.10–0.20, alfakR 0.38–0.60.
   - Fig. 6 (dominant karyotype changes with the missegregation rate): anchored 0.80, state_space 0.82, the same as the reference 0.80; auto 0.77, pairwise 0.69, two time points 0.60.
5. **Cost**: state_space is the most expensive. Median time of all-time-point fits on ABM 804 seconds (anchored 313, pairwise 419, alfakR 3525);
   Salehi longest 8.5 hours, 13.6 GB (SA609_X10, 10 passages). All tasks finished, with no timeout or out-of-memory.
6. **Fraction not accepted by the convergence check**: state_space is highest on long windows (ABM w8_all 49%, Salehi all 38%), pairwise lowest (24%, 14%).
7. **auto (rule-based choice of one of the three methods, thresholds 0.3 / 0.5) did not pass the criteria registered in advance** (section 9); 3 of 13 failed:
   - missegregation-rate experiment 8-passage window r 0.804 vs 0.815 (allowed difference 0.01);
   - Salehi Wasserstein forecast 0.506 vs 0.541 (p = 0.004);
   - sweep Euclidean forecast 0.417 vs 0.438 (a difference of 1 of 48 daughters, allowed difference 0).

   auto was therefore not adopted; the alfak2 package removed this option, and the default is still anchored.
   On Salehi new karyotypes auto is better than anchored (0.607 vs 0.521, p = 6e-5).

## 2. Setup

| Fit | Input |
|---|---|
| alfak2 2 time points | The first and last time points of the window or trajectory |
| alfak2 all time points | All time points, default `series = "anchored"` |
| alfak2 pairwise all time points | All time points, `series = "pairwise"` |
| alfak2 state space all time points | All time points, `series = "state_space"` |
| alfak2 auto all time points | For each input, one of the three all-time-point fits above chosen by rule, no separate fit; evaluated, not adopted (section 9) |
| alfakR 2 time points / all time points | First and last time points / all time points (the paper's approach) |
| paper alfakR | The paper's own fits, reference only |

- Data, windows, minobs and parameters are the same as in 8512694 (`docs/report_8512694.md` section 2).
- On inputs with only two time points (ABM 2-passage windows, Salehi 2-passage trajectories), the three methods are the same two-time-point fit.
- Accuracy is compared only on the same inputs where all arms have results (paired); failures, time and memory are counted over all inputs.
- Sweep subset: Salehi trajectories where the anchored fit removed intermediate time points ("all time points (sweeps)" in `summary_comparison.R`).

## 3. Accuracy (`figures/summary_comparison.png`, `analysis/summary/summary_accuracy.tsv`)

All time points, paired medians or win fractions; in brackets, the fraction of pairs better than anchored (ABM landscape) or the p of the paired test.

| Test | anchored | pairwise | state_space | alfakR all |
|---|---|---|---|---|
| ABM landscape r, 4-passage window (1177 pairs) | 0.66 | 0.57 (32%) | 0.66 (65%) | 0.40 |
| ABM landscape r, 8-passage window (1190 pairs) | 0.71 | 0.63 (30%) | 0.72 (65%) | 0.45 |
| ABM rescaled R², 4 / 8 | 0.41 / 0.49 | 0.30 / 0.38 | 0.42 / 0.51 | 0.12 / 0.18 |
| Missegregation-rate experiment r, 4 / 8 (125 pairs) | 0.76 / 0.82 | 0.69 / 0.74 (19% / 10%) | 0.77 / 0.83 (67% / 71%) | 0.65 / 0.80 |
| ABM beats baseline at passage 10, Euclidean, 4 / 8 | 0.86 / 0.91 | 0.76 / 0.78 | 0.84 / 0.89 | 0.50 / 0.53 |
| ABM beats baseline at passage 10, Wasserstein, 4 / 8 | 0.87 / 0.92 | 0.72 / 0.79 | 0.85 / 0.91 | 0.45 / 0.49 |
| Salehi forecast beats parent, Wasserstein (283 pairs) | 0.54 | 0.52 (p = 0.24) | 0.49 (p = 0.007) | 0.27 |
| Salehi forecast beats parent, Euclidean | 0.42 | 0.45 (p = 0.11) | 0.42 (p = 0.9) | 0.26 |
| Salehi new karyotypes: fitness predicts appearance (351 pairs) | 0.52 | 0.66 (p = 4e-9) | 0.55 (p = 0.24) | 0.24 |
| Sweeps: forecast beats parent, Wasserstein (48 pairs) | 0.64 | 0.66 (p = 0.67) | 0.63 (p = 0.86) | 0.11 |
| Sweeps: forecast beats parent, Euclidean | 0.44 | 0.48 (p = 0.29) | 0.39 (p = 0.43) | 0.11 |
| Sweeps: new karyotypes (92 pairs) | 0.47 | 0.58 (p = 0.03) | 0.45 (p = 0.62) | 0.11 |

- Paired tests of the three alfak2 methods against alfakR all time points give p ≤ 0.006, always in the direction of alfak2 being better, with only one exception:
  on the 8-passage window of the missegregation-rate experiment, pairwise is better than alfakR in only 35% (p = 6e-5, alfakR better).
- The two-time-point results are identical to 8512694 (same fits, section 6).

## 4. The paper's tests, by arm (two designs)

### 4.1 ABM and missegregation-rate experiment (`analysis/abm/`)

**Supp. Fig. 2a**: median landscape Pearson r, minobs 10, by λ = 0.2 / 0.4 / 0.8 / 1.6.

| Fit | 4-passage window | 8-passage window |
|---|---|---|
| alfak2 2 time points | 0.20 / 0.50 / 0.82 / 0.85 | 0.21 / 0.61 / 0.85 / 0.86 |
| alfak2 all (anchored) | 0.18 / 0.47 / 0.80 / 0.86 | 0.15 / 0.54 / 0.84 / 0.90 |
| alfak2 pairwise | 0.20 / 0.45 / 0.77 / 0.80 | 0.18 / 0.52 / 0.80 / 0.79 |
| alfak2 state space | 0.20 / 0.50 / 0.80 / 0.86 | 0.18 / 0.57 / 0.85 / 0.88 |
| alfakR 2 time points | 0.09 / 0.39 / 0.66 / 0.46 | 0.10 / 0.48 / 0.72 / 0.56 |
| alfakR all time points | 0.12 / 0.38 / 0.65 / 0.63 | 0.15 / 0.48 / 0.75 / 0.74 |
| paper alfakR | 0.09 / 0.38 / 0.66 / 0.59 | 0.11 / 0.46 / 0.76 / 0.73 |

- Frequent-karyotype r (`rfq`) on the 8-passage window: state_space 0.98 / 0.96 / 0.95 / 0.92, close to two time points;
  anchored only 0.85 at λ = 0.2, pairwise 0.84 at λ = 0.4.

**Supp. Fig. 2b, 2f** (`supp2f_flows.tsv`, 8-passage window): fraction of frequent karyotypes with R² > 0 at λ = 0.2:
anchored 0.58, pairwise 0.81, state_space 0.86, alfakR all time points 0.77, paper 0.83; at λ ≥ 0.4 all alfak2 arms are ≥ 0.98.

**Supp. Fig. 2c–e** (`supp2c_example.tsv`, `supp2de_example.tsv`, w_1p6_m_5e-05_rep_74):
- 2c, 8-passage window, minobs 10, r between true and estimated fitness of frequent karyotypes: anchored 0.955, state_space 0.957, pairwise 0.925,
  alfak2 two time points 0.891, alfakR all time points 0.831, paper 0.849. 2-passage window: alfak2 0.46, alfakR −0.24, paper −0.46.
- 2d–e, frequency trajectories of the 8 most common karyotypes, mean absolute error between the median projection and the observations: on the 8-passage window alfakR 0.004,
  anchored 0.011, state_space 0.012, pairwise 0.026. alfak2's projection error is larger than alfakR's.

**Supp. Note 3.2** (`note32_initial_frequencies.tsv`): the estimate of the initial frequencies has a median correlation ≥ 0.9998 on all arms,
median RMSE 0.002–0.003.

**Supp. Fig. 3** (`supp3_misseg_pearson.tsv`, 8-passage window, missegregation rate 5e-5 / 1e-4 / 2e-4 / 4e-4 / 8e-4):

| Fit | Median r |
|---|---|
| alfak2 2 time points | 0.84 / 0.88 / 0.88 / 0.85 / 0.83 |
| alfak2 all (anchored) | 0.83 / 0.85 / 0.84 / 0.77 / 0.70 |
| alfak2 pairwise | 0.80 / 0.77 / 0.77 / 0.68 / 0.62 |
| alfak2 state space | 0.85 / 0.85 / 0.85 / 0.79 / 0.77 |
| alfakR 2 time points | 0.68 / 0.76 / 0.83 / 0.80 / 0.81 |
| alfakR all time points | 0.72 / 0.80 / 0.81 / 0.79 / 0.81 |

  At high missegregation rates, alfak2's all-time-point fits are all worse than two time points; at missegregation rate 8e-4 all three all-time-point fits are below alfakR (0.70 / 0.62 / 0.77 vs 0.81).

**Supp. Fig. 4** (`supp4b_cv_flags_accuracy.tsv`, 8-passage window): fraction with landscape R > 0 when CV > 0; at λ ≥ 0.4
all alfak2 arms are 1.00; at λ = 0.2 anchored 0.36, pairwise 0.30, state_space 0.41 (higher, in fact, when CV < 0: 0.47 / 0.41 / 0.54),
CV has no diagnostic value for alfak2 at λ = 0.2; alfakR all time points 0.52 vs 0.46, paper 0.10 vs 0.21.

**Fig. 1g / Supp. Fig. 5b** (`fig1g_beat_baseline.tsv`, selection by CV > 0, passage 10, 8-passage window, λ = 0.4 / 0.8):
anchored 0.83 / 0.98, state_space 0.85 / 0.96, pairwise 0.71 / 0.81, alfak2 two time points 0.78 / 0.91,
alfakR all time points 0.70 / 0.79, paper 0.63 / 0.73.

### 4.2 Salehi (`analysis/salehi/`)

| Test (paper's selection rule) | alfak2 2 | anchored | pairwise | state_space | alfakR 2 | alfakR all | Paper |
|---|---|---|---|---|---|---|---|
| Fig. 2b median CV score (each on its own only) | 0.86 | 0.85 | 0.88 | 0.92 | 0.05 | 0.47 | 0.48 |
| Fig. 2c prediction angle, median of per-cluster medians (°) | 72 | 55 | 52 | 74 | 83 | 80 | 71 |
| Fig. 2d Wasserstein distance of prediction and observation, Spearman | 0.25 | 0.33 | 0.57 | 0.56 | 0.45 | 0.70 | 0.69 |
| Fig. 2f beats baseline, Wasserstein | 0.49 | 0.56 | 0.62 | 0.61 | 0.25 | 0.30 | 0.38 |
| Fig. 2f beats baseline, Euclidean | 0.36 | 0.50 | 0.53 | 0.51 | 0.25 | 0.24 | 0.36 |
| Number of fits (Fig. 2c–f) | 53 | 52 | 53 | 49 | 44 | 46 | 45 |
| Fig. 3e–f fitness significantly predicts new karyotypes, paper's rule | 0.50 | 0.49 | 0.55 | 0.50 | 0.23 | 0.15 | 0.32 |
| Fig. 3e–f same, all fits | 0.57 | 0.53 | 0.63 | 0.55 | 0.21 | 0.23 | 0.38 |

- Fig. 2c: the angles of all arms are significantly below 90° (cluster bootstrap p ≤ 0.011); unrelated lineages 90.4° (p = 0.80), diverging sister lineages 79°.
- Fig. 3 paired (`fig3_paired.tsv`, 503 regressions): pairwise vs anchored, significant only for pairwise 59, significant only for anchored 10 (McNemar p = 7.5e-9);
  state_space vs anchored 33 vs 24 (p = 0.29).
- Paired forecasts over all minobs (`paired_forecast.tsv`, 401–408 pairs): the three alfak2 methods vs alfakR all time points p ≤ 2e-7;
  pairwise vs anchored is significant only on the overlap metric (worse, p = 0.002); state_space vs anchored is significant only on Wasserstein
  (worse, 0.47 vs 0.51, p = 0.007), consistent with section 3.

## 5. Fig. 4–6: can the paper's biological conclusions be obtained (`analysis/fig456/conclusions.tsv`)

- **Method**: each arm uses its own landscapes, selects fits by the paper's rule (xv > 0, the best minobs per trajectory), and redoes the paper's analyses:
  - Fig. 4: 200 bootstrap draws of non-overlapping lineages, Δf of the 44 single-chromosome neighbours of the frequent karyotypes, glmmTMB and lmer models;
  - Fig. 5: piecewise lmer and permutation KS test (1e4 permutations) on SA906;
  - Fig. 6: the best landscape per clean_id, steady states from `find_steady_state()` at 8 missegregation rates.
- **Reference**:
  - Fig. 4 is the paper's published bootstrap result;
  - Fig. 5 and 6 are the paper's fits recomputed with the same code of this repository ("paper alfakR, every input").
    Fig. 6 uses the recomputation instead of the published values because the `find_steady_state()` of alfakR dev2 differs from the version used in the paper:
    with the same fits and the same karyotypes, the steady states differ by up to 0.18, and only 26 of 35 landscapes have the same dominant karyotype at every missegregation rate (`fig6_reproduction.tsv`).
- **Exclusion**: 24 trajectories that cannot be analysed are removed for all arms (`FIG456_EXCLUDE=1`); the reproduction of the paper's fits on all inputs is computed separately (`*_all_inputs`).
- **"Agree"**: same sign as the reference, and the same significance.

| Arm | Conclusions that agree (19 in total) |
|---|---|
| paper alfakR (recomputed with this code) | 19 |
| paper alfakR, all inputs | 19 |
| alfakR all time points | 15 |
| alfakR 2 time points | 14 |
| alfak2 auto (not adopted) | 13 |
| alfak2 2 time points | 11 |
| alfak2 all (anchored) | 10 |
| alfak2 state space | 10 |
| alfak2 pairwise | 9 |

**Fig. 4d** (`figures/fig4_effects.png`; bootstrap medians, * means the 95% interval excludes 0):

| Coefficient | Paper | alfak2 2 | anchored | pairwise | state_space | auto | alfakR 2 | alfakR all |
|---|---|---|---|---|---|---|---|---|
| PDX_v (PDX vs in vitro, \|Δf\|) | +1.23* | −1.2* | −1.3* | −1.4* | −0.38* | −1.2* | +0.69* | +1.5* |
| PDX_kv (same, variance) | +1.72* | −0.56 | −0.78* | −1.1* | +0.58 | −0.71* | +1.2* | +1.7* |
| cisplatin_v | +0.92* | +0.36* | +0.24 | +0.26 | +0.28* | +0.35* | +0.42 | +0.94* |

- The cause lies in vitro (`fig4_deltaf.tsv`, median |Δf|):
  - in vitro: alfak2 0.027–0.051 (auto 0.024), paper 0.006, alfakR 0.004–0.009;
  - PDX control: alfak2 0.013–0.017, paper 0.012, alfakR 0.013–0.021.
  - In vitro, alfak2 gives the single-chromosome neighbours of the frequent karyotypes very large fitness differences (and mostly positive, the neighbours being worse than the frequent karyotypes), while in PDX it is comparable to the paper.
- Data differences between in vitro and PDX: median interval between adjacent passages 37.5 days vs 15 days, median cells per time point 602 vs 316.
  So it is not a "short interval" problem; the cause is to be investigated on the alfak2 side.
- Fig. 4f (similarity model): the alfak2 arms agree with the paper on `pair_typesame_traj` and its interaction terms;
  `pair_typeparallel` and its interaction terms are significant in alfak2 but not in the paper, counted as disagreement.

**Fig. 5** (`fig5_tests.tsv`, `fig5_model.tsv`):

| Arm | Number of fits | KS D | Permutation p |
|---|---|---|---|
| paper alfakR | 48 | 0.39 | **0.030** |
| alfak2 2 time points | 49 | 0.36 | 0.20 |
| alfak2 all (anchored) | 49 | 0.36 | 0.13 |
| alfak2 pairwise | 52 | 0.32 | 0.10 |
| alfak2 state space | 36 | 0.40 | 0.15 |
| alfak2 auto (all of SA906 chose anchored, identical to anchored) | 49 | 0.36 | 0.13 |
| alfakR 2 time points | 34 | 0.23 | 0.60 |
| alfakR all time points | 59 | 0.23 | 0.38 |

- alfak2's D is close to the paper's, but none is significant; alfakR dev2's D is much smaller.
- Fig. 5b piecewise model: each arm's dk has the same sign as the paper and is significant, and so is dk:cut3 (only alfakR two time points is not significant).
  `dk:wgd_strWGD+:cut3TRUE` has the opposite sign to the paper in anchored, pairwise, state_space and alfakR all time points;
  `dk:wgd_strWGD+` has the opposite sign in pairwise and state_space. Intercepts cannot be compared across methods (the zero point of fitness differs).

**Fig. 6** (`fig6_dominance.tsv`): fraction of landscapes whose dominant karyotype changes with the missegregation rate:

| Arm | Landscapes | Fraction changed |
|---|---|---|
| Reference (paper's fits recomputed) | 35 | 0.80 |
| Paper, published | 35 | 0.77 |
| alfak2 all (anchored) | 39 | 0.80 |
| alfak2 state space | 38 | 0.82 |
| alfak2 pairwise | 39 | 0.69 |
| alfak2 auto | 39 | 0.77 |
| alfak2 2 time points | 40 | 0.60 |
| alfakR all / 2 time points | 36 / 35 | 0.75 / 0.77 |

  All arms reach the conclusion "the dominant karyotype of most landscapes changes". But landscape by landscape, the fraction whose dominant karyotype agrees with the reference at all missegregation rates is:
  alfak2 arms 0%, alfakR all time points 22%, two time points 3% (`fig6_agreement.tsv`).

## 6. Failures, time, memory (`analysis/summary/summary_performance.tsv`)

| Data | Group | Method | Failures | Median time, seconds (90% quantile; max) | Median peak memory (max), MB |
|---|---|---|---|---|---|
| ABM | All time points | anchored | 0/2400 | 313 (971; 2441) | 1802 (6582) |
| ABM | All time points | pairwise | 0/2400 | 419 (1146; 3257) | 1932 (6261) |
| ABM | All time points | state_space | 0/2400 | 804 (3816; 18538) | 2371 (8538) |
| ABM | All time points | alfakR | 33/2400 | 3525 (8728; 25142) | 12577 (88804) |
| Missegregation rate | All time points | anchored | 0/250 | 310 (1237; 2589) | 2536 (6714) |
| Missegregation rate | All time points | pairwise | 0/250 | 550 (1355; 2869) | 2694 (7209) |
| Missegregation rate | All time points | state_space | 0/250 | 1045 (7806; 22514) | 3472 (12104) |
| Missegregation rate | All time points | alfakR | 0/250 | 4426 (18839; 30045) | 25734 (147156) |
| Salehi | All time points | anchored | 15/561 | 228 (622; 1508) | 1587 (4624) |
| Salehi | All time points | pairwise | 15/561 | 244 (741; 1838) | 1725 (5926) |
| Salehi | All time points | state_space | 15/561 | 663 (4156; 30720) | 2327 (13622) |
| Salehi | All time points | alfakR | 15/561 | 1586 (5111; 12465) | 8003 (41358) |

- The two-time-point rows are identical to 8512694.
- Failures are counted per input × minobs, with excluded inputs removed; the 15 on Salehi are all "no frequent karyotypes at minobs 20", the same for all four methods.
- All 5617 tasks COMPLETED, with no TIMEOUT or OUT_OF_MEMORY. The longest is state_space on SA609_X10_l_10: 8.5 hours,
  13 GB in the c7f809c cost test.
- Scoring (277 Salehi trajectories) takes ≤ 8 minutes per trajectory.

## 7. Removed time points (`analysis/summary/summary_dropped.tsv`)

| Method | Analysable Salehi multi-passage trajectories | Trajectories with time points removed | Time points removed | Fell back to two time points |
|---|---|---|---|---|
| anchored | 187 | 53 | 103 | 10 |
| pairwise | 187 | 42 | 42 | 4 |
| state_space | 187 | 42 | 42 | 4 |

- On ABM and the missegregation-rate experiment, none of the three methods removed time points.
- The 24 excluded trajectories are the same for all methods (anchored's rule). 4 of them could in fact be fitted under pairwise's rule,
  but they are still excluded so that all arms use the same inputs.
- `drift` of state_space: ABM median 0 (max 0.0025); Salehi median 0, 90% quantile 0.010, max 0.16.
  On ABM state_space hardly uses drift.

## 8. Identical to 8512694 (`analysis/compare_rounds/results_8512694_summary.tsv`)

05a8855 did not change the anchored and two-time-point code paths. `scripts/compare_rounds.R` checks the 8071 alfak2 fits one by one (two time points and anchored all time points):

- Status is identical; landscapes are identical within numerical rounding (maximum difference of the means 2.2e-16, of the standard deviations 1.8e-15).
- Cross-validation output, history and forecast replicates are all identical.
- The largest score difference is 0.039, like the 0.038 of the previous round (8512694 vs ff13cf4); it is the order of tied values in the Spearman correlation, not different fits.

## 9. The auto arm (`analysis/summary/auto_criterion.tsv`, `tables/auto_picks.tsv`)

**Rule**: fixed by the multi-time-point session before seeing any E4 result. The thresholds were chosen on e1–e3, and the selection criterion was also registered before scoring:
among rules that pass the e1 gate and the e3 neutral test, take the one with the lowest e2 ratio, and within 0.01 take the simpler one.
The provisional 0.3 / 0.8 rule has an e2 ratio of 0.822; the adopted 0.3 / 0.5 has 0.810.

- `ends_major`: the fraction of cells of the intermediate time points that fall on the major karyotypes of the first and last time points, taking the minimum;
- `weight`: `intermediate_weight` of the anchored fit;
- `ends_major` ≥ 0.3 → anchored; otherwise if `weight` is NA (no intermediate time points left) or ≥ 0.5 → pairwise; otherwise state_space.

**Method**: auto does not fit separately.

1. `scripts/auto_picks.R 0.3 0.5` chooses a method for each multi-time-point input;
2. `scripts/compose_auto.R` links the chosen fit directories as method `alfak2_auto`, 1255 links in total, none missing;
3. scored and summarised like the other arms, and compared paired with anchored on the same inputs.

**Choices**:

| Data | anchored | pairwise | state_space |
|---|---|---|---|
| ABM (w4_all, w8_all) | 800 | 0 | 0 |
| Missegregation-rate experiment | 157 | 93 | 0 |
| Salehi | 80 | 102 | 23 |

- 28 inputs have `weight` NA.
- The intermediate time points of ABM all fall on the first and last major karyotypes, so all choose anchored, and on ABM auto is bitwise identical to anchored.
- All 36 SA906 trajectories also choose anchored, so auto in Fig. 5 is identical to anchored.

**Criteria** (registered by the multi-time-point session before E4):

| Criterion | Pairs | auto | anchored | Difference | Allowed | Result |
|---|---|---|---|---|---|---|
| ABM landscape r, 4-passage window | 1177 | 0.660 | 0.660 | 0 | −0.01 | Pass |
| ABM landscape r, 8-passage window | 1190 | 0.709 | 0.709 | 0 | −0.01 | Pass |
| Missegregation-rate experiment landscape r, 4-passage window | 125 | 0.754 | 0.764 | −0.010 | −0.01 | Pass (exactly at the boundary) |
| Missegregation-rate experiment landscape r, 8-passage window | 125 | 0.804 | 0.815 | −0.011 | −0.01 | **Fail** |
| ABM beats baseline at passage 10 (Euclidean / Wasserstein, 4 / 8) | 1177 / 1190 | Same as anchored | | 0 | −0.02 | Pass (4 criteria) |
| Salehi forecast beats parent, Euclidean | 283 | 0.410 | 0.422 | −0.012 | −0.02 | Pass |
| Salehi forecast beats parent, Wasserstein | 283 | 0.506 | 0.541 | −0.035 | −0.02 | **Fail** |
| Sweeps: forecast beats parent, Euclidean | 48 | 0.417 | 0.438 | −0.021 | 0 | **Fail** |
| Sweeps: forecast beats parent, Wasserstein | 48 | 0.639 | 0.639 | 0 | 0 | Pass |
| Salehi new karyotypes: fitness predicts appearance | 351 | 0.607 | 0.521 | +0.086 | −0.02 | Pass |

**Conclusion: auto fails** (3 of 13 fail), **evaluated, not adopted**.

- The user's rule is: auto is implemented only if the criteria pass. The multi-time-point session therefore removed `series = "auto"` from the alfak2 package, and the default is still anchored:
  - 4ff8906 added auto;
  - 94cef55 ("Call the columns time points, and drop the automatic choice") removed it.

  From 94cef55 on, alfak2 has no auto option. The anchored, pairwise and state_space code of 94cef55 is identical to 05a8855; only text was changed.
- The rule and the results here are recorded as a negative result in alfak2's `validation/sparse_joint/series_methods/README.md` and in the methods document (tex).
- No new round is needed: removing auto does not change the anchored, pairwise and state_space fits, and the terminology commit only changes text.

- **Missegregation-rate experiment, 8-passage window**: the difference −0.011 is the difference of medians rounded to two decimals, only 0.001 beyond the allowed value.
  Paired, it is a real difference: auto is worse on 32% of the inputs and better on 6% (p = 4e-5). The cause is the 93 inputs that chose pairwise;
  pairwise is the weakest on the missegregation-rate experiment anyway (section 3).
- **Salehi, Wasserstein**: worse pairs 4.6%, better 1.4% (p = 0.004), also a real difference.
  auto's 0.506 is not ahead of pairwise (0.52) or state_space (0.49) used alone.
- **Sweeps, Euclidean**: the difference −0.021 is exactly 1 of 48 daughters (p = 0.85). It fails only because the allowed value of this criterion is 0.
- **The good side**: Salehi new karyotypes auto 0.607, better than anchored 0.521 (p = 6e-5); on the sweep subset 0.554 vs 0.47 (p = 0.08).
  This comes from choosing pairwise, consistent with pairwise being the best on new karyotypes in section 3.

**The auto arm in Fig. 4–6** (tables of section 5; reruns 21953562, 21953563):

- Reproduces 13 of the paper's 19 conclusions, the most of the alfak2 arms (the others 9–11).
  The 3 more than anchored (10) are the two cisplatin coefficients of Fig. 4d (cisplatin_v +0.35*, cisplatin_kv +0.62*, same sign as the paper and significant)
  and the intercept of Fig. 4f (not significant, agreeing with the paper).
- PDX vs in vitro in Fig. 4d is still opposite to the paper (PDX_v −1.22*); median in-vitro |Δf| 0.024, 4 times the paper's.
- Fig. 5: identical to anchored.
- Fig. 6: fraction of landscapes whose dominant karyotype changes 0.77 (30 / 39), reference 0.80; landscape by landscape, agreement with the reference is still 0 / 31.
- These are not among the registered criteria and do not change the "fail" conclusion.

## 10. Convergence check and CV scores

**Fraction not accepted by the convergence check** (minobs 10, `tables/{abm,salehi}_runs.Rds`):

| Data | anchored | pairwise | state_space | Two time points |
|---|---|---|---|---|
| ABM w4_all / w8_all | 21% / 41% | 20% / 24% | 23% / 49% | w4 20%, w8 36% |
| Missegregation rate w4_all / w8_all | 11% / 25% | 16% / 21% | 14% / 31% | w4 10%, w8 19% |
| Salehi all | 20% | 14% | 38% | 17% |

**CV scores** (`analysis/summary/summary_cv.tsv`): R²R is relative to each method's own full-data estimate and cannot be compared across methods; only the median and the fraction above 0 are reported.

| Data | anchored | pairwise | state_space | alfakR all |
|---|---|---|---|---|
| ABM all time points | 0.56, 72% | 0.12, 58% | 0.51, 73% | −0.01, 43% |
| Missegregation rate all time points | 0.89, 96% | 0.52, 93% | 0.84, 98% | 0.07, 78% |
| Salehi all time points | 0.80, 77% | 0.84, 82% | 0.86, 77% | 0.04, 56% |

## 11. Run record and decisions

| Job | Content |
|---|---|
| 21945965, 21945966 | c7f809c build and cost test (54 tasks), superseded by 05a8855 |
| 21946021 | Build alfak2 05a8855 |
| 21946022–21946042 (even) | Memory tiers for ABM, missegregation rate and Salehi, 5617 tasks in total; each with an add-hours monitor (id + 1) |
| 21951672 | Build glmmTMB 1.1.11 into `lib/extra_r` (not in the image) |
| 21951796, 21952078 | Salehi scoring of 277 trajectories; the state_space fit of SA609_X10 finished only after scoring, rescored separately |
| 21952074 | ABM summary, Supp. Fig. 2c–e, k-set check, comparison with 8512694 |
| 21952127 | Salehi summary, Fig. 5 |
| 21952128–21952131 | Fig. 4 (7 arms × 10 batches, plus 10 batches for the paper's fits on all inputs), Fig. 6 (7 arms plus 1) |
| 21952227, 21952262 | Summary figure, resources, CV; Fig. 4–6 summary |
| 21953043, 21953044 | auto: `auto_picks.R 0.3 0.5`, `compose_auto.R` |
| 21953046 | Salehi scoring for auto (277 trajectories), add-hours monitor 21953049 |
| 21953047, 21953048 | auto arm of Fig. 4 and 6, first attempt: read a scoring table that did not yet have the auto rows (collect_salehi ran only after them), discarded |
| 21953051 | Full summary (including auto_criterion, Fig. 5, Fig. 4–6 summary, comparison with 8512694) |
| 21953562, 21953563 | Rerun of the auto arm of Fig. 4 and 6, add-hours monitor 21953564 |
| 21953565 | Figures and analysis tables of each round redrawn with the "time points" terminology, without refitting: 8512694 and 05a8855 redrawn (old ones kept in `pre_time_points/`); ff13cf4 changed in content when redrawn with the current code, so the original figures were restored and the redrawn version put in `redrawn_2026-11-15/`; 8327282 lacks the input coverage table and was not redrawn |
| 21953580 | Fig. 5 tables (`fig5_tests*.tsv`) rewritten with the new labels, numbers unchanged |

Decisions (the multi-time-point session authorised an automatic run, decided according to my suggestions):

- The round of record is 05a8855; the pair check of c7f809c would block 42 Salehi sweep inputs; it was only a cost test.
- No separate test round: the c7f809c cost test already covered the two new methods, and 05a8855 only changes the time-point selection of pairwise / state_space.
- Memory tiers: pairwise 16 GB, state_space 32 GB, 12 hours per task, the add-hours monitor checks every hour.
- The reference for Fig. 6 was changed to the recomputation of the paper's fits (reason in section 5); the published values are kept as a row.
- When choosing the best minobs for Fig. 4, ties in xv take the smallest minobs (the paper's fits have no ties; the results are unchanged).
- The auto thresholds were set by the multi-time-point session (0.3 / 0.5) and sent before auto was scored; before that, no E4 accuracy number was disclosed to it.
