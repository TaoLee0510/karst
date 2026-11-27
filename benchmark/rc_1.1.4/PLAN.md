# KARST 1.1.4 release candidate: missegregation-rate grid and Kester

Approved by the user on 2026-11-26 at about 01:50 EST ("as recommended, confirmed";
options D0a–D5a). All computation runs on the HPC: qos xxlarge, the usual
excluded nodes, one fit per Slurm task, and an add-hours monitor per job.
Development of the patch versions (`patch/`) has stopped.

## 0. Freeze

- **The release candidate:**
  - KARST 1.1.4 committed as 50e5889 on dev and tagged `v1.1.4-rc1`
    (pushed);
  - package files only: DESCRIPTION, NEWS.md, R/alfak-compat.R,
    R/sparse-joint-own.R, two Rd files, and two tests, which pass locally
    together with test-alfak-compat.R.
- **Library:** `benchmark/paper/lib/karst_50e5889` on the HPC clone, built by
  `benchmark/paper/hpc/build_karst.sh` from `git archive` of 50e5889.
  - sha256: 3428bab09a2577599099a50923931bb55321101d4caa2bdd2b59ab7acec4b917
  - lock: `benchmark/paper/hpc/locks/karst_50e5889.lock`
- **Reproduction check (`rc_check.R`; 60 fits, one per task):**
  - 20 KARST fits of the paper round 6ae2a24 (Salehi, both designs,
    minobs 10);
  - 20 of the β grid (`results_6ae2a24/beta_grid`, ×2–×64);
  - 20 own-benchmark held-out datasets (refit as the round does, compared
    with the stored predictions).

  Each is refitted with 50e5889. Pass: the objective, the summary table and
  the predictions are identical to the stored ones, i.e. a largest absolute
  difference of 0. Otherwise stop and ask the user.

## Grid (D1a)

The missegregation rate is 5e-5 × {1, 2, 4, 8, 16, 32, 64}, the same rate
for fitting and forecasting (as `benchmark/paper/docs/beta_grid_salehi_6ae2a24.md`).
×1 is the round itself.

## 1. Paper benchmark (D2a)

- Salehi Fig. 2c / 2d for KARST and ALFA-K at every multiple:
  - from the existing β-grid fits and analyses (1.1.3 = 1.1.4 by default,
    confirmed by the reproduction check);
  - no new fits.
- **Deliverable:** one table per test, over both designs (ends and all) and
  every multiple:
  - Fig. 2c without selection (paired, bootstrap by line);
  - Fig. 2c by the paper rule;
  - Fig. 2d by the paper rule;
  - Fig. 2d without selection (raw and noise held fixed);
  - failure rates.

  All are read from `results_6ae2a24/beta_grid/analysis` and the per-rate
  tables.

## 2. Own benchmark (D3a)

**Runs:**
- Main table (distances 0–2) of the held-out split (500 landscapes; 7,000
  datasets):
  - KARST default (`karst_sparse_joint_balance`) on all 7,000;
  - ALFA-K minobs 10 on the core design (4,000);
- refitted at each multiple ×2–×64 (×1 = round 6ae2a24);
- with the library of 50e5889.

**Code and output:** the benchmark code gets one option, a multiple of the
fitting rate (`--beta-multiple`). Output goes to
`benchmark/results/sparse_two_timepoint_50e5889_x<mm>/`. The data, nodes and
scoring are the round's.

**Estimate:** about 180 CPU-hours per multiple, about 1,100 in total.

**Report:** per multiple and wavelength × distance, the round's scores (RMSE,
Spearman, direction, coverage) for both methods and the zero prediction,
with the change against ×1. The simulations' true rate is 5e-5, so every
other multiple is a misspecified rate.

## 3. Kester (D4a; `dev/KESTER_PREDICTION_PROMPT_CLAUDE.md`)

**Data:** the five `*_karst_counts.tsv` matrices on RED, one trajectory per
culture, fitted separately.
- **Time axis:** one passage = 7 days, because the study split the cultures
  weekly (`dev/kester_sampling/README.md`).
  - KARST keeps its default birth rate of 1/day;
  - ALFA-K gets the same time axis (`dt` = 7).

**Methods:**
- KARST 50e5889 (`fit_sparse_joint`, `min_cn` = 0, as for Salehi);
- ALFA-K (alfakR dev2 f93b2d2, `alfak()` as in `benchmark/paper/R/fit.R`).
- Both at minobs 5, 10 and 20 and the 7 multiples, with n0 / nb as for
  Salehi (1e5 / 1e7).

**Salehi-like (APKS1–3):**
- **Folds:**
  - APKS1 P10/P22 → P25;
  - APKS2 P5/P8 → P14 and P5/P8/P14 → P18;
  - APKS3 P3/P7 → P15, P3/P7/P15 → P19 and P3/P7/P15/P19 → P23;
  - plus the endpoint-only version (first and last training passage) of
    the three folds with ≥ 3 training passages.
- **Forecasts:**
  - as the paper: `alfakR::predict_evo` ABM, 5 replicates, from the last
    training passage to the held-out one;
  - plus the deterministic engine on list E (the training karyotypes and
    their one- and two-step neighbours), from the same start, for a direct
    comparison of the two landscapes.
- **Scores:**
  - Fig. 2c angle of the centroid change, W1 and overlap against the
    held-out passage;
  - per-cell log score against persistence (the last training passage) and
    neutral (flat fitness in the same engine and rate).
  - The predictions are smoothed with weight ε on a uniform over E. ε is the
    share of the last training passage's cells whose karyotype was not seen
    at the earlier training passages, floored at 0.01, so it is set from
    training data only. Held-out cells outside E are counted and reported
    as outside the support.
  - Fig. 2d-style rank correlations over the folds are descriptive only (6
    folds).

**GBM-like (A3, AP1):**
- **Folds:** seeded three-fold cell holdout within each passage
  (`karst_validation_folds()` of kflInfer `Fitting_and_validation.R`, seed
  20261126). The same partitions are used for both methods; training uses the
  remaining cells of both passages.
- **Scores (`karst_score_counts()` of the same file),** on the held-out cells
  at each passage, with the later passage reported separately:
  - KARST's native expected composition (`expected_t*`);
  - both methods propagated with the deterministic engine from the training
    composition of the first passage;
  - persistence and neutral.
- Labelled "within-harvest cell holdout scores".

**Output:** RED `/share/lab_crd/taoli/Project/Kester/results/forecast_50e5889/`:
- the fold manifest;
- training and target cell counts;
- predicted and observed compositions;
- paired scores;
- fit and forecast diagnostics and failures;
- versions and seeds.

Code is in `benchmark/kester/`; the data stay on RED.

## Order and budget (D5a)

- Order: 0 → 1 → 3 → 2, automatic, with a short report after each part.
- Stop if the reproduction check fails.
- Budget about 1,150 CPU-hours; ask the user before 1,400.
