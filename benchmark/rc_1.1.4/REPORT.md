# KARST 1.1.4 rc: missegregation-rate grid on the paper benchmark, the own benchmark and Kester

2026-11-26. Plan: `PLAN.md` (approved 2026-11-26, options D0a–D5a).

**Grid:** the missegregation rate is 5e-5 × {1, 2, 4, 8, 16, 32, 64}. Each
method uses the same rate for fitting and forecasting.

**Release candidate:** KARST 50e5889 (1.1.4 rc, tag `v1.1.4-rc1`). Its
default fits are identical to 1.1.3 (6ae2a24): the reproduction check
refitted 60 fits and all 60 matched exactly.

| Part | Result |
|---|---|
| 0. Freeze | Done: 50e5889, tag `v1.1.4-rc1`, pushed; reproduction check 60/60 identical |
| 1. Paper benchmark (Salehi Fig. 2c / 2d) | `results/paper_grid/paper_grid.md` (from the β-grid study, no new fits) |
| 2. Own benchmark | This report, section 2; tables in `results/own_grid/summary/` |
| 3. Kester | `../kester/REPORT.md` |

## 1. Paper benchmark (summary of `results/paper_grid/paper_grid.md`)

**Fig. 2c without selection (paired median angle):**
- KARST's angle is at or below ALFA-K's at every multiple and design.
- The difference is significant only at ×64 with all time points (8.2°).

**Fig. 2c by the paper rule:**
- No multiple gives a significant difference.
- KARST is closer at ×1 (by 11° with two time points and 14° with all).
- From ×4 on, ALFA-K is mostly closer.

**Fig. 2d by the paper rule (Spearman):**
- KARST rises with the multiple:
  - two time points: 0.22 at ×1 to 0.74 at ×64;
  - all time points: 0.28 at ×1 to 0.80 at ×64.
- ALFA-K stays between 0.27 and 0.71.
- KARST is ahead from ×8 (two time points) or ×32 (all time points), never
  significantly.
- ALFA-K is significantly ahead with all time points at ×1–×4.

**Fig. 2d without selection:** KARST − ALFA-K lies between −0.04 and +0.14,
and every interval includes 0.

## 2. Own benchmark

### Setup

**Fits:**
- The held-out main table (distances 0–2): 500 landscapes, 7,000 datasets.
- KARST default (`karst_sparse_joint_balance`) on all 14 designs.
- ALFA-K minobs 10 on the 8 core designs (4,000 datasets).
- Each refitted at ×2–×64, one fit per Slurm task.
- ×1 is round 6ae2a24 itself.

**Inputs and scoring:**
- The round's datasets, nodes, seeds and scoring, unchanged.
- The simulations' true rate is 5e-5, so ×1 is the correct rate and every
  other multiple is a misspecified one.

**Aggregation (as the round report):**
- Metrics are averaged over designs within a landscape, then over landscapes.
- Error ratio = mean RMSE / mean null RMSE.
- Direction accuracy is computed on one-copy pairs, with ties counted as 1/2
  (`own_direction.R`).
- Bootstrap intervals resample landscapes within wavelength (2,000 draws).

**Check:** the ×1 values reproduce the round's `own_summary.json` exactly for
both methods.

**Two departures from PLAN.md:**
- The rate multiple is applied by a wrapper (`own_task.R`:
  `dataset$beta × multiple`) instead of a new `--beta-multiple` option, so the
  benchmark code is unchanged.
- The fits are in `/share/lab_crd/taoli/Project/karst_rc_1.1.4/results/own_grid/x<mm>/`
  instead of `benchmark/results/sparse_two_timepoint_50e5889_x<mm>/`.

### Completeness

| Multiple | ×1 | ×2 | ×4 | ×8 | ×16 | ×32 | ×64 |
|---|---|---|---|---|---|---|---|
| KARST failed / 7,000 | 0 | 3 | 1 | 1 | 8 | 69 | 156 |
| ALFA-K m10 failed / 4,000 | 145 | 145 | 145 | 145 | 145 | 145 | 145 |

**KARST failures:**
- They grow with the multiple.
- Every one is in the 1,080-day design (`long_500_700_1080d`), mostly at
  λ ≥ 1.2.
- The message is always: "No candidate history (1080, 3000) could be used: the
  posterior mode at the fitted scales could not be found at any."

**ALFA-K failures:** the same datasets fail as in the round, at every
multiple.

Failed fits are scored as the zero prediction (`results/own_grid/summary/karst_failures.tsv`).

### Core designs, all wavelengths (500 landscapes weighted equally)

**Error ratio** (lower is better; zero prediction = 1):

| Multiple | KARST d0 | d1 | d2 | ALFA-K d0 | d1 | d2 |
|---|---|---|---|---|---|---|
| ×1 | **0.553** | **0.440** | **0.406** | 5.889 | 3.732 | 3.130 |
| ×2 | 0.902 | 0.604 | 0.497 | 5.881 | 3.727 | 3.126 |
| ×4 | 1.648 | 1.090 | 0.858 | 5.895 | 3.734 | 3.133 |
| ×8 | 2.647 | 1.710 | 1.319 | 5.927 | 3.755 | 3.151 |
| ×16 | 3.998 | 2.494 | 1.908 | 5.998 | 3.799 | 3.190 |
| ×32 | 4.975 | 3.011 | 2.365 | 6.164 | 3.897 | 3.260 |
| ×64 | 5.453 | 3.183 | 2.641 | 6.516 | 4.103 | 3.408 |

**Direction accuracy of one-copy changes** (higher is better; chance = 0.5):

| Multiple | KARST d0 | d1 | d2 | ALFA-K d0 | d1 | d2 |
|---|---|---|---|---|---|---|
| ×1 | 0.835 | 0.832 | 0.805 | 0.441 | 0.453 | 0.428 |
| ×4 | 0.830 | 0.833 | 0.806 | 0.467 | 0.479 | 0.457 |
| ×8 | 0.829 | 0.833 | 0.807 | 0.485 | 0.498 | 0.480 |
| ×16 | 0.828 | 0.831 | 0.806 | 0.503 | 0.520 | 0.507 |
| ×64 | 0.817 | 0.801 | 0.782 | 0.535 | 0.561 | 0.558 |

**Spearman correlation with the truth:**

| Multiple | KARST d0 | d1 | d2 | ALFA-K d0 | d1 | d2 |
|---|---|---|---|---|---|---|
| ×1 | **0.471** | **0.363** | **0.285** | 0.165 | 0.129 | 0.107 |
| ×8 | 0.412 | 0.295 | 0.231 | 0.194 | 0.151 | 0.123 |
| ×64 | 0.318 | 0.225 | 0.173 | 0.220 | 0.172 | 0.144 |

**KARST's 95% interval coverage:** 0.976 / 0.990 / 0.984 at ×1, and 0.988 to
1.000 at ×2 and above.

**Per wavelength** (`own_grid_summary.md`):
- ×1 is KARST's best error ratio in every wavelength × distance cell, with
  one exception: at λ 0.8, distance 2, ×2 gives 0.565 against 0.570.
- At ×2, KARST is still below 1 everywhere except λ 0.4, distance 0 (1.095).
- From ×4 on, distance 0 is above 1 at every λ.

### What the overstated rate does to KARST: a downward bias, not a reshuffle

Core designs, all wavelengths (`core_scale_diag.tsv`). Bias = mean
(prediction − truth), relative to the reference karyotypes:

| Multiple | KARST bias d0 / d1 / d2 | KARST error sd d0 / d1 / d2 |
|---|---|---|
| ×1 | −0.003 / +0.015 / +0.028 | 0.046 / 0.068 / 0.078 |
| ×2 | −0.060 / −0.072 / −0.071 | 0.053 / 0.071 / 0.080 |
| ×4 | −0.134 / −0.176 / −0.180 | 0.071 / 0.076 / 0.083 |
| ×8 | −0.225 / −0.290 / −0.295 | 0.097 / 0.084 / 0.088 |
| ×32 | −0.434 / −0.513 / −0.481 | 0.163 / 0.099 / 0.105 |

For scale, the true spread (sd) is 0.053 / 0.070 / 0.080, and the truth lies
0.076 / 0.164 / 0.222 below the reference on average.

**How to read this:**
- At an overstated rate, KARST puts every karyotype other than the frequent
  reference ones too low, and the shift grows with the multiple.
- The scatter around that shift grows much less, and the order of the
  karyotypes barely moves:
  - direction accuracy is flat to ×16 and falls by at most 0.03 at ×64;
  - Spearman falls by 0.11–0.15 at ×64.
- So the error ratio and the direction accuracy separate. The error ratio
  rises mainly because of the level shift.
- This is what the model must do to fit the same counts at a higher rate:
  missegregation would then produce the rarer karyotypes faster, so their low
  observed counts have to be explained by strong selection against them.

**ALFA-K barely responds:**
- Its error ratio changes by at most +0.63 (×64, d0; interval includes 0).
- Its error is scatter, not bias: error sd about 0.50 against a true spread
  of 0.05.
- Its direction accuracy rises from 0.43–0.45 at ×1 to 0.54–0.56 at ×64,
  still near chance.

### KARST against ALFA-K at the same multiple (core designs, paired by landscape)

| Multiple | Error ratio, KARST − ALFA-K, d1 [95% CI] | Landscapes where KARST is lower (d0 / d1 / d2) | Direction, KARST − ALFA-K, d1 | Landscapes where KARST is higher |
|---|---|---|---|---|
| ×1 | −3.29 [−4.88, −2.20] | 100% / 100% / 100% | +0.380 | 100% |
| ×2 | −3.12 [−4.79, −2.03] | 97% / 99% / 100% | +0.369 | 100% |
| ×4 | −2.64 [−4.22, −1.57] | 54% / 59% / 73% | +0.354 | 100% |
| ×8 | −2.05 [−3.61, −0.94] | 33% / 34% / 36% | +0.335 | 100% |
| ×16 | −1.31 [−2.95, −0.15] | 19% / 23% / 27% | +0.311 | 100% |
| ×64 | −0.92 [−2.57, +0.34] | 13% / 20% / 21% | +0.241 | 100% |

**Pooled mean against the typical landscape:**
- KARST's pooled error ratio is lower at every multiple.
- From ×8 on, however, ALFA-K has the lower error ratio on most landscapes.
- The pooled mean is carried by ALFA-K's large errors on a minority of
  landscapes (λ 1.6: error ratio 28).

**Direction:** KARST stays ahead on 95–100% of landscapes at every multiple
and distance.

### All 14 designs and the stress designs (KARST)

| Multiple | All designs d0 / d1 / d2 | Stress designs d0 / d1 / d2 |
|---|---|---|
| ×1 | 0.569 / 0.450 / 0.416 | 0.586 / 0.463 / 0.429 |
| ×2 | 0.877 / 0.608 / 0.509 | 0.849 / 0.614 / 0.524 |
| ×8 | 2.423 / 1.670 / 1.314 | 2.184 / 1.621 / 1.307 |
| ×64 | 5.152 / 3.157 / 2.524 | 4.831 / 3.125 / 2.376 |

The paired change against ×1 is positive at every multiple and distance, and
its interval excludes 0 (`paired_bootstrap.tsv`).

## 3. Kester

See `../kester/REPORT.md`.

**Salehi-like folds:**
- KARST forecasts the held-out passage better than ALFA-K on all 6
  all-passage folds at every multiple.
- It beats persistence.
- It ties with neutral at ×1 and beats neutral from ×4 up.

**GBM-like cell holdouts:** KARST is best on both cultures.

**Effect of the multiple:** it changes KARST's scores only slightly; the data
are too small to separate the rates.

## Reading across the three parts

**Own benchmark (true rate known):**
- KARST is most accurate at the true rate.
- Every overstated rate shifts its fitness values down, more so the higher
  the multiple, while their order holds.
- An overstated rate is not a free gain where the fitness values themselves
  are the output.

**Salehi (true rate unknown):**
- Higher multiples raise KARST's Fig. 2d ranking, never significantly.
- This agrees with the earlier varyMisseg study
  (`benchmark/paper/docs/fig2d_misseg_6ae2a24.md`): there, a fit rate below
  the true one reproduces the Salehi reversal.

**Kester:** too small to tell the rates apart.

**What the grid does not do:** it does not pick a rate. The rate still has
to come from outside the data, and 1.1.4 rc keeps 5e-5 as the default.

## Compute

**This plan:** 1,396 CPU-hours (allocated cores × elapsed time; monitors
excluded). That is the 1,400 line; nothing more was submitted.

**Own-benchmark fits:**
- 1,295 CPU-hours: KARST 981, ALFA-K 314.
- PLAN.md estimated about 1,100. KARST's fits at ×2–×4 took longer (median
  39–40 s against 32 s at ×1).

**Scoring and summaries:** 83 CPU-hours. This includes two scoring attempts
that failed on my script errors (a missing function and a relative output
path), about 8 CPU-hours.

## Files

**Local (`benchmark/rc_1.1.4/`):**
- Scripts:
  - `own_task.R` (fit);
  - `own_score.R` and `own_direction.R` (scores);
  - `own_grid_summary.R` and `own_grid_diag.R` (tables);
  - `hpc/` (job scripts).
- Results: `results/own_grid/summary/`.

**HPC (`/share/lab_crd/taoli/Project/karst_rc_1.1.4/`):**
- Fits: `results/own_grid/x02` … `x64`.
- `metrics_heldout_grid.csv` (419,640 rows) and `direction_heldout_grid.csv`
  (231,000 rows).
- The summary folder, logs and `jobs.tsv`.
