# Pre-fit preprocessing v6: correcting the missegregation rate with the level of the neighbour cloud

**Fixed in writing:** 2026-11-20, before any v6 computation.

**User confirmation:** 2026-11-20 around 20:30, "Confirmed, proceed with everything as you recommend".
- two candidate statistics A, B, chosen on the development half;
- the external test uses only the paper's 125 missegregation-rate simulations, with no NK etc. added;
- the Salehi prediction criterion uses (i): β_used for the fit, 5e-5 for the prediction; (ii) is only reported;
- decision thresholds 0.80, 10%, 50%, 5%, 0.95, 1.02, 60%, 0 / −0.1;
- steps 1–3 automatic; if they pass, continue with 4, 5; if not, stop and report.

**Origin:** the lower bound of v5 (`report_v5.md`) almost never triggers. Afterwards it was seen that the level of the neighbour cloud rises monotonically with the true missegregation rate.
- At steady state, neighbour-cloud ratio ≈ missegregation rate × number of routes ÷ fitness disadvantage of the neighbours per division.
- The disadvantage is unknown and can only be calibrated on simulations.
- The rule follows v5: only raise, capped at 10 times.

**Not done:**
- No change to the KARST package, no commit.
- No new simulations generated; the sealed held-out landscapes (seed 98000) are not used.

## Data

**v3 simulations:**
- The non-injected version in the v4 screening record, 1632 datasets.
- Split into development and confirmation halves, following `v4_half()`: calibration on the development half, testing on the confirmation half.

| True multiple | Set it is in |
|---|---|
| ×1 | e1, e2, e3 |
| ×3, ×10, ×30 | C2 |
| ×4, ×16 | e4 |

**Paper simulations:** `PAPER_DATA/raw/ABM_varyMisseg`, 125, used only for testing.
- Missegregation rates 5e-5, 1e-4, 2e-4, 4e-4, 8e-4 (×1, ×2, ×4, ×8, ×16), 25 each; GRF wavelength 0.8.
- Two time points: the windows of the paper benchmark, taking the first and last of the last 2, 4, 8 sampling time points before day 120 (`abm_window_times()`).
- Cell counts:
  - each time point originally has 1000 cells;
  - a 300-cell version is also made: at each time point 300 are drawn without replacement from the 1000, with the seed determined by the dataset name.
- 750 datasets in total.

**Real data:**
- Salehi: 155 trajectories with later passages; the first and the last passage are taken.
- NUMBAT: 41, only the summary is reported.

## Statistics (two candidates)

**Notation:**
- a, b: counts at the two time points;
- number of routes G_k = Σ_i c_i · (1 + 1[c_i ≥ 2]), the same as v5;
- "adjacent" means L1 distance 1.

**A (the user's original idea):**
- High-frequency karyotype: ≥ 10 cells at at least one time point.
- Neighbour cloud: karyotypes adjacent to some high-frequency karyotype that have < 10 cells at both of their own time points, minus those the envelope judges as significantly rising (`env_analyse()`, |z| > 3 and z > 0).
- Computation:
  - the cell count C of the neighbour cloud and the cell count H of the high-frequency karyotypes, both summed over the two time points;
  - Ḡ_A: the G of the high-frequency karyotypes, weighted mean by a + b;
  - **S_A = (C + 0.5) / H / Ḡ_A**.
- Without a high-frequency karyotype, S_A cannot be computed and there is no raise.

**B (no threshold):**
- Within one time point:
  - ordered cell pairs of the same karyotype I_t = Σ_k y_k(y_k − 1);
  - ordered cell pairs of adjacent karyotypes N_t = Σ_k y_k · Σ_{j adjacent to k} y_j.
- Ḡ_B: the G of all karyotypes, weighted mean by a + b.
- **S_B = (N₀ + N₁ + 0.5) / (I₀ + I₁ + 0.5) / Ḡ_B**.

## Step 1: compute the statistics (HPC, no fits)

- One task per v3 folder, the paper simulations split into 5 tasks, Salehi and NUMBAT one each.
- Output: `results/realdata/noise/v6/stats/<tag>.rds`.

## Step 2: development half, choose the statistic and calibrate (local)

**Discrimination threshold:**
- S separates datasets with true ×1 from those with ≥ ×4 at AUC ≥ 0.80.
- Datasets where A cannot be computed are treated as the minimum value.
- If both A and B pass, take the one with the higher AUC (compared after rounding to three decimals); if equal, take A. If neither passes, stop.

**Calibration:**
- Model: log₂ S = c + d · log₂(true multiple) + ε, ordinary least squares, on all datasets in the development half where S can be computed. σ is the residual standard deviation.
- Requires d > 0 and a 95% interval that excludes 0. If not met, stop.
- Invert, taking the conservative 25% prediction quantile: x = (log₂ S − c − 0.674σ) / d.
- **β_used = 5e-5 × 2^min(max(x, 0), log₂ 10)**. When S_A cannot be computed, β_used = 5e-5.
- Report only: the slopes fitted separately on the e sets and on C2.

## Step 3: test (local, no fits)

The confirmation half (v3) and the paper simulations must each pass:

| Criterion | Requirement |
|---|---|
| False raise | among datasets with true ×1, those with β_used ≥ 2 times the default are ≤ 10% |
| Useful | among datasets with true ≥ ×4, those with β_used ≥ 2 times the default are ≥ 50% |
| No overshoot | among all datasets, those with β_used > 2 × the true missegregation rate are ≤ 5% |

- Paper-simulation datasets with true ×2 enter only "no overshoot".
- If it fails, stop and write the report.
- Report only: the fraction raised by set, cell count and interval; the β_used distribution of Salehi and NUMBAT.

## Step 4: fits (HPC)

**Scope:** datasets in the confirmation half and the paper simulations with β_used > 5e-5. Those not raised are the same as T0 and are not fitted.

**Fits:**
- T0: β = 5e-5; T6: β = β_used. Library karst_v4, no `noise_phi` given.
- The v3 settings are the same as v4's `fit_task.R`; the existing v4 T0 fits are used directly.
- Paper simulations:
  - fit settings the same as the paper benchmark: `fit_sparse_joint(x, dt, beta)`, the rest default;
  - truth: the paper's GRF fitness (`grf_fitness()`) times the missegregation retention factor (`sm_self_factor()`);
  - scoring nodes as in v4: `stt_scoring_nodes()`.

**Error ratio:** the mean of RMSE_T6 / RMSE_T0 over distances 0, 1, 2.

**Criteria (computed separately for the two sets):**
- True missegregation rate > ×1 and raised: median error ratio ≤ 0.95, and the 95% upper limit of a cluster bootstrap (2000 resamples) < 1. v3 clusters by landscape, the paper simulations by simulation.
- True ×1 and raised: median error ratio ≤ 1.02. With fewer than 10, only reported.

## Step 5: Salehi (HPC)

**Parallel lineages:**
- Trajectory pairs: same cell line, same treatment at both ends, no shared time points, the same as v4.
- Among trajectory pairs in which at least one has β_used ≥ 2 times, those whose two values differ by at most 2 times are ≥ 60%.
- With fewer than 10 qualifying trajectory pairs, only reported.

**Later-passage prediction:**
- Scope: trajectories with β_used > 5e-5, T0 and T6 fits, settings the same as v4's `salehi_fit.R`.
- Later passages ≤ 60 days, `ll_fc_s_e2`, the support set shared by the three approaches.
- Criterion, using (i) T6 fit with prediction missegregation rate 5e-5, against T0: in the confirmation group (SA906, SA532, SA535, SA1035), gain per cell ≥ 0, and the 95% lower limit of a cluster bootstrap by trajectory > −0.1.
- With fewer than 10 raised trajectories in the confirmation group, only reported.
- (ii) T6 fit, prediction missegregation rate β_used: only reported. The development group is only reported.

For NUMBAT only the summary of β_used is reported.

## Outputs

- Code: `dev/realdata/noise/R/cloud_level.R`, `dev/realdata/noise/v6/`.
- Results: HPC `results/realdata/noise/v6/`, local copy in `dev/results/realdata/noise/v6/`.
- Report: `dev/realdata/noise/report_v6.md`.
