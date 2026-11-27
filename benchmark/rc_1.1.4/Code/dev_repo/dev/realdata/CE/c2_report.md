# Phase C: simulations matched to Salehi (C2)

The plan is in section "C" of "Phases C–E" in `../PLAN.md` (fixed 2026-11-20 00:35 EST). Results are on HPC in `karst_dev/results/realdata/CE/c2/`; the acceptance tables are local in `dev/results/realdata/CE/c2/accept_*.rds`. "a / b / c" is the minimum / median / maximum across datasets.

## Conclusions

| regime | corresponds to | result | settings adopted |
|---|---|---|---|
| H | SA609 | **accepted** (round 1) | ABM missegregation rate 5e-4, bottleneck 3000, pre-evolution 720 days |
| M | confirmation group | **accepted** (round 1) | ABM missegregation rate 1.5e-4, bottleneck 3000, pre-evolution 720 days |
| W | SA906a (WGD) | **failed in both rounds, stopped** | — |

H and M have 48 datasets each (12 landscapes × 4 intervals); the cd342db base fits (ce_base) all succeeded. W stopped and was reported as the plan requires; the criterion that depends on W (W1 of D3) is not automatically judged as failed; the user decides.

## Simulator and design

- **ABM**: `dev/abm/abm_c2.cpp`, a copy of R10's `abm_wgd.cpp` (exported names `c2_abm_cpp`, `c2_landscape_fitness`). Three changes, none of which changes the division and missegregation processes:
  1. Passage by time: every 15 days (150 steps), after that step's division and after recording, the population is binomially diluted to a mean of B cells;
  2. Between two passages the population is capped at 2e5; above that it is diluted to the cap (culling_survival is no longer used);
  3. R10's ploidy-aware GRF is turned on only in regime W; H and M use the GRF of the dev benchmark.
- **landscape**: the GRF of the dev benchmark (10 waves, `stt_generate_centroids`, diploid founder division rate > 0.4), 6 each at wavelengths 0.8 and 1.2 (seeds 64001–64012).
- **Trajectories**: one per landscape. Starting from B diploid cells, pre-evolved with the passages above for 720 days (t₀); the whole population is recorded before each passage from t₀ to t₀ + 120 days (9 in total).
- **Datasets**: 300 cells taken at t₀, 450 cells at t₀ + 15, 30, 45, 60 days (`stt_sample_counts`). The true frequencies of later passages are in the trajectory files, for distant tests.
- **Truth**: the expected daily growth rate of a karyotype in the ABM, minus the part that leaves through missegregation at the simulated missegregation rate: log(1 + f(1 − 2p_any)Δt)/Δt, p_any = 1 − (1 − β_sim)^(copy number) (in W, 2.5β_sim for ploidy 4). All fits use beta 5e-5 (the same as the Salehi benchmark), which differs from the simulated missegregation rate.
- **Fits**: `c2_fit.R`, `fit_sparse_joint(counts, dt, beta = 5e-5, birth_rate = 1, min_cn = 1, max_cn = 8)` (10 for W), the rest default. Scoring nodes the same as in the dev benchmark (observed karyotypes, plus at most 2000 neighbours each at distances 1 and 2).

## Round 1 (one grid shared by H and M)

Grid: ABM missegregation rate {1.5e-4, 5e-4, 1.5e-3} × bottleneck {300, 3000}, pre-evolution 720 days. 48 datasets per setting.

| setting | fraction of single-cell karyotypes at the last passage | model-domain coverage (worst time point) | Pearson dispersion | failures (model domain) | falls in H | falls in M |
|---|---|---|---|---|---|---|
| 1.5e-4, B 300 | 0.12 | 0.95 | 2.40 | 0 | no | no (coverage > 0.92) |
| **1.5e-4, B 3000** | **0.23** | **0.82** | **2.11** | 0 | no | **yes** |
| 5e-4, B 300 | 0.39 | 0.60 | 2.26 | 0 | no | yes |
| **5e-4, B 3000** | **0.64** | **0.31** | **1.50** | 0 | **yes** | no |
| 1.5e-3, B 300 | 0.79 | 0.08 | 1.81 | 5 | no (coverage < 0.10) | no |
| 1.5e-3, B 3000 | 0.95 | 0.02 | 5.40 | 13 | no | no |

- Ranges: H single-cell [0.34, 1.0], coverage [0.10, 0.40], dispersion [0.86, 2.5]; M single-cell [0.10, 0.50], coverage [0.55, 0.92], dispersion [2.0, 6.0].
- Among the settings that fall within the ranges, the one closest to the target (H: 0.64, 0.19, 1.65; M: 0.25, 0.85, 4.18) is taken (sum of |log ratio| of the three medians): H has only one; M has two, and 1.5e-4/B 3000 is taken (distance 0.79, the other 1.40).
- **The dispersion of M, 2.11, is only slightly above the lower limit 2.0**, and below the confirmation group's 4.18. With bottleneck 300 the dispersion is not higher, because once single-cell karyotypes become more numerous, the Pearson dispersion is diluted.

**Details of the two regimes adopted**

| regime | interval (days) | fraction of single-cell karyotypes | model-domain coverage | dispersion | number of karyotypes | karyotypes shared by the two time points |
|---|---|---|---|---|---|---|
| H | 15 / 30 / 45 / 60 | 0.64 / 0.64 / 0.62 / 0.63 | 0.34 / 0.30 / 0.31 / 0.27 | 1.37 / 1.43 / 1.58 / 1.54 | 525–552 | 53 / 44 / 37 / 35 |
| M | 15 / 30 / 45 / 60 | 0.26 / 0.22 / 0.24 / 0.22 | 0.81 / 0.84 / 0.84 / 0.84 | 1.92 / 2.13 / 2.16 / 2.33 | 230–246 | 57 / 54 / 48 / 46 |

- Ranges across datasets: H single-cell 0.37 / 0.64 / 0.85, coverage 0.01 / 0.31 / 0.67, dispersion 0.27 / 1.50 / 6.82; M single-cell 0.11 / 0.23 / 0.40, coverage 0.59 / 0.82 / 0.95, dispersion 1.20 / 2.11 / 3.40.
- Landscapes with wavelength 1.2 are more diverse: in H, coverage is 0.47 at wavelength 0.8 but only 0.10 at 1.2.
- Refitting with `regimes/H` and `regimes/M` (ce_base, 96 tasks) gives medians identical to round 1.

**Performance of base on C2 (report only, not an acceptance criterion)**: per dataset, RMSE ÷ the RMSE of predicting 0 everywhere, median:

| regime | distance 0 | distance 1 | distance 2 |
|---|---|---|---|
| H | 1.22 | 1.78 | 1.64 |
| M | 1.04 | 1.02 | 1.01 |

In a regime as diverse as SA609, the cd342db base fit is worse than "predict 0 everywhere", most clearly for unobserved neighbours; in a regime like the confirmation group, it is level with 0.

## W (WGD): failed in both rounds

- Settings: M's missegregation rate 1.5e-4, bottleneck 3000, pre-evolution 720 days; ω = 1e-4; missegregation rate ×2.5 for ploidy 4; ploidy-aware GRF; copy number cap 10. Gain g in the order of the plan: round 1 0.005, round 2 0.
- Acceptance: median fraction of ploidy-4 cells at the first time point in [0.02, 0.5], and a rise of ≥ 0.1 in ≥ 25% of datasets.

| round | g | ploidy-4 fraction at the first time point (min / 25% / median / 75% / max) | datasets rising ≥ 0.1 | verdict |
|---|---|---|---|---|
| W1 | 0.005 | 0.003 / 0.007 / **0.013** / 0.063 / 0.197 | **6%** (90% quantile of the rise 0.036, maximum 0.25) | failed |
| W2 | 0 | 0 / 0.003 / **0.010** / 0.028 / 0.98 | **0%** | failed |

**Cause (recorded as it is)**: the direction of adjustment in the plan, "first 0.005, then 0", came from the lesson of R10 (in R10's ABM, with β = 5e-5 and g = 0.02, ploidy 4 filled the first time point). C2's missegregation rate is 3 times that of R10, and ploidy 4 is multiplied by a further 2.5, so tetraploids lose more through missegregation, and lowering g to 0 can only make tetraploids fewer. To get rising tetraploids like SA906a, the needed direction is to increase g (or reduce the missegregation multiplier for ploidy 4); this is not among the adjustments the plan allows, so by the rules it stopped, with no further tuning.

## Deviations and notes

1. **One trial run (local) was made before the round 1 grid**: one trajectory (L01, missegregation rate 1.5e-3, bottleneck 300) and fits of two of its datasets, single-cell karyotypes 0.77–0.85, coverage 0.03–0.16, showing that diversity was too high. On this basis the not yet submitted grid was changed from {5e-4, 1.5e-3, 4e-3} to {1.5e-4, 5e-4, 1.5e-3}. This trial run is not counted as a round. Before W1 was submitted, L01 was also trial-run locally (same settings as W1); nothing was changed on that basis.
2. The pre-evolution length (adjustable under the plan) was fixed at 720 days; the adjustment was not used.
3. Fits use `min_cn = 1`: zero copies are not viable in the ABM (the plan says "0 or 1, as appropriate").
4. The W datasets have been generated (`grid/rW1_g0.005`, `grid/rW2_g0`), and the task tables `tasks/c2_fit_W1.tsv` and `c2_fit_W2.tsv` are written, but they were not fit because they failed acceptance.

## Files and jobs

- Code: `dev/abm/abm_c2.cpp`; `dev/realdata/CE/c2_sim.R` (design, truth, statistics), `c2_make.R` (trajectories and datasets), `c2_fit.R` (fits and scoring), `c2_tasks.R`, `c2_accept.R`, `c2_promote.R`, `c2_regime_tasks.R`, `c2_check_regimes.R`; settings table `tasks/c2_settings.tsv`.
- Data: `c2/regimes/{H,M}/{datasets,nodes}/` (48 each, `SETTING` records the setting adopted); trajectories `c2/grid/<setting>/truth/`.
- base fits: `c2/fits/ce_base/{H,M}/` (and for each grid setting).
- Jobs: round 1 simulation 22139890, fits 22140147; W1 simulation 22140581; W2 simulation 22140632; base fits of H and M 22140652. Every array had an add-hours monitor in the small QOS; all completed, with no failures or timeouts.
