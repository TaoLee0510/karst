# Pre-fit preprocessing v4: sensitivity-first noise screening, and two treatments when noise is present

The user confirmed on 2026-11-20: "Confirmed, proceed with everything as you suggest". Four points were confirmed specifically:
- **Noise term of option 2:** added to the variance of the counts (negative binomial), not to the missegregation rate.
- **Clustered anomalies (selection-like):** still go through the normal pipeline.
- **The four numbers used for the verdict:** 1.02, 0.02, 0.95, −0.1.
- **Simulated data:** use the v3 simulations, split into development and confirmation halves; no new ones are generated.

The user's statement of the purpose:
- This step assesses whether the input data are disturbed by strong noise. Detection power comes first; false positives are acceptable.
- When data are judged noisy, help KARST find a suitable treatment, instead of rigidly holding the default missegregation rate.

This file was written at 2026-11-20 16:05 EST, before any v4 computation. It is not changed afterwards; corrections go into the "Corrections" section at the end.

## Step 1: Screening

**Input:** the first and last time points of each dataset, with ρ set to 1 division per day. NUMBAT keeps the 0.141 of N1, report only.

**Method:**
- Envelope, z and Moran's I use `R/envelope.R`; the code is not changed.
- Add `R/envelope_v4.R`, which computes the anomalies separately for each β on the grid (13 values, 5e-5 to 1e-2, i.e. `env_betas`):
  - the expected frequency q(β) and the variance term are the same as in `env_analyse()`;
  - unseen karyotypes likewise get 3/n0 added;
  - the definitions of up, down, z, |z| > 3 and excess are all the same as for the envelope, except that the q of this one β is used, not the maximum and minimum across β.
- "Has anomalies" at a given β means ≥ 3 anomalous karyotypes and excess ≥ 0.05 (as for the envelope).

**Four classes (judged in order):**

| class | condition | destination |
|---|---|---|
| selection-like | envelope verdict selection (has anomalies and Moran p < 0.05) | normal pipeline (T0) |
| noise | envelope verdict sampling or undetermined, i.e. not explained even by the envelope's upper limit, and not clustered or the test cannot be done | T1, T2 |
| needs higher β | envelope verdict no_anomaly, but has anomalies at β = 5e-5 | report T1 only |
| normal | the rest | T0 |

**β̂ (used for T1):**
- Noise class: the β with the fewest anomalous karyotypes; on ties, the smaller one.
- Needs-higher-β class: the smallest β with no anomalies; when every β has anomalies, likewise the one with the fewest anomalous karyotypes.

**φ̂ (used for T2):**
- φ̂ = max(1, mean of the envelope z² of the testable karyotypes).
- z and "testable" are taken from the output of `env_analyse()`; z inside the envelope is 0.

## Step 2: Three fits

| name | settings | used on which datasets |
|---|---|---|
| T0 | β = 5e-5 (current) | noise class, needs-higher-β class |
| T1 (option 1) | β = β̂ | noise class; needs-higher-β class report only |
| T2 (option 2) | β = 5e-5, counts at the second time point negative binomial, variance μ·φ̂ | noise class |

**Other settings:**
- All use two time points (first and last); other arguments as in each set's original fits:
  - group e: birth_rate 1, copy number 1–8;
  - C2: copy number 1 to the simulation's maximum;
  - Salehi: the settings of the paper benchmark.
- For datasets with injection, the fits use the counts after injection, exactly the same as in the screening (saved during screening).

**Library:** `karst_v4`, based on 1.1.3 (6ae2a24), with two patches added:
- `patches/R1_nb1_drift.patch`: checked, can be applied directly to 6ae2a24;
- `patches/V4_noise_phi.patch`: `fit_sparse_joint(noise_phi = φ)`, all states nb_phi = φ − 1; Poisson when φ = 1 or not given.

**Checks (stop if they disagree):**
- Without `noise_phi`, `karst_v4` is bit-identical to `karst_6ae2a24`;
- with `noise_phi = 1`, also bit-identical;
- each checked on at least 2 datasets in group e and in C2.

## Step 3: Simulation evaluation (the v3 simulations)

**Data:**
- e1–e4 and C2, 300 and 1000 cells per time point.
- Three injections: no injection, independent overdispersion with φ = 20, region-wise sampling with σ = 1.0.
  - Seeds as in v3 (`seed_of(dataset_id, variant)`);
  - e2 and e4 had no injection before; it is added this time.

**Split into halves (fixed before any fit):**
- **Group e:**
  - development: nz_L01, L02, L05, L07, L08, L10;
  - confirmation: nz_L03, L04, L06, L09, L11, L12;
  - in each half, 1 or 2 per wavelength class, and 3 each of e4's ×4 (odd) and ×16 (even).
- **C2:**
  - development: nz_c2_L01–L03, L07–L09;
  - confirmation: nz_c2_L04–L06, L10–L12;
  - in each half, 3 each of wavelength 0.8 and 1.2.
- **The development half is used only to get the runs working and fix errors; no parameter is tuned.** The verdict looks only at the confirmation half and is evaluated once.

**Metrics:**
- Truth and scoring nodes as in the autorun phases:
  - group e uses `sm_true_net_rate` and the nodes generated by v3;
  - C2 uses `c2_true_rate` and `stt_scoring_nodes` (seed = dataset seed);
  - scoring uses `stt_score_predictions`.
- **Error ratio:**
  - for each dataset and each distance d0, d1, d2, compute RMSE_T / RMSE_T0, and average the three.
  - This equals the ratio of the "RMSE ratios relative to predicting 0 everywhere", because the denominator is the same within a dataset.
- **False selection (e3):** among d0–d2 nodes, the fraction with |predicted relative fitness| > 1.96 × predicted sd, one value per dataset.

## Step 4: Salehi (and the NUMBAT screening)

**Data:**
- The 155 Salehi trajectories with later passages (`beta/decomp/tasks/fit_ext.tsv`), two time points (ends).
- The 41 NUMBAT datasets are screened only, report only.

**Fits and predictions:**
- Noise-class trajectories are fit with `karst_v4` as T0, T1, T2, written in alfak format as in the paper benchmark (minobs 10), then scored with `forecast_score.R` (unchanged):
  - forecast missegregation rate 5e-5 (×1), nb 1e7;
  - the common support is the union of the three;
  - only later passages ≤ 60 days are used, taking `ll_fc_s_e2`.
- Groups: development group SA609, SA039; confirmation group SA906, SA532, SA535, SA1035.
- Cluster bootstrap by trajectory, 2000 times.

## Verdict fixed in advance (confirmation half; each treatment judged separately)

| criterion | target | requirement |
|---|---|---|
| safety 1 | no injection, e1, e2, e4, C2 flagged as noise | median error ratio T/T0 ≤ 1.02 |
| safety 2 | no injection, e3 flagged as noise | mean false-selection fraction of T exceeds that of T0 by no more than 0.02 |
| useful | with injection (φ = 20 or σ = 1.0), flagged as noise, all five sets | median error ratio ≤ 0.95, and upper limit of the 95% interval from a cluster bootstrap by landscape (2000 times) < 1 |
| real data | Salehi confirmation-group trajectories flagged as noise, each later passage | mean log score gain (T − T0) ≥ 0, and lower limit of the 95% interval > −0.1 |

**Selection:**
- Passing requires all four to pass.
- When both treatments pass, take the one with the smaller "useful" median; when only one passes, use it.
- When neither passes, the screening serves only as a data-quality report, and KARST keeps its defaults.
- The development group (Salehi) only checks direction.

**Report only:**
- Fraction in each class: by set, injection, cell number, bottleneck;
- φ̂ against the injected φ;
- T1 for the needs-higher-β class: settings in C2 whose true β is above the default;
- all results of the development half;
- the NUMBAT classification;
- whether parallel lineages agree.

## Scale and order

1. Screening: locally or on HPC, a few minutes.
2. Patches, build, checks.
3. First a timed trial run on 1 landscape of the development half.
4. The whole development half; fix errors.
5. The confirmation half in one run.
6. Salehi run together with the confirmation half.

Estimate: about 2700 datasets flagged as noise × 3 fits, plus the needs-higher-β class, about 8–9 thousand fits; one fit per task.

## Rules

- Code in `dev/realdata/noise/v4/` and `dev/realdata/noise/R/envelope_v4.R`, results in `results/realdata/noise/v4/`; no commit, the package is not changed.
- HPC: xxlarge, request ≤ 12 hours; the add-hours monitor uses small, 3 days; extra arguments separated by colons.
- Do not use the sealed landscapes (seed 98000).
- After `report_v4.md` is written, stop and wait for the user.

## Stop rules

- The patch checks disagree;
- a script fails for an unknown reason;
- > 5% of the tasks of an array fail, and not because of resources;
- a case the rules do not cover.

## Corrections
- **2026-11-20 16:24 EST: add a robust φ̂ as an extra candidate T2r (user confirmed "add a robust φ̂ as an extra candidate").**
  - **Time:** before any T2 / T2r fit. At that point only the bit-identity checks had been run.
  - **Reason:** after the screening, the originally defined φ̂ (mean of z² outside the envelope) was found to be dominated by a few karyotypes with very large deviations:
    - with injection φ = 20, median about 100, maximum about 680,000;
    - without injection, median 26–56;
    - Salehi noise class, median 47, maximum 3735.
    - Such a φ̂ would make the second time point almost ineffective.
  - **Settings of T2r:** β = 5e-5, `noise_phi` = φ̂_r; the rest the same as T2.
  - **Definition of φ̂_r:**
    - at β̂ (noise class: the β with the fewest anomalous karyotypes), compute the full signed z for every testable karyotype:
      - z = (y1 − n1 q̃) / sqrt(n1 q̃ (1 − q̃) + n1² v);
      - q̃ = q(β̂), with 3/n0 added for karyotypes not seen at the first time point; the variance term as in `env_by_beta()`;
      - not set to 0 inside the envelope.
    - "Testable": y1 ≥ 5, or seen at the first time point and n1 q ≥ 5.
    - φ̂_r = max(1, median(z²) / 0.4549), where 0.4549 is the median of the chi-squared distribution with 1 degree of freedom.
    - When fewer than 5 karyotypes are testable, φ̂_r = 1 (T2r is then the same as T0).
  - **Verdict:** T2r uses the same four criteria. Of the three candidates (T1, T2, T2r), among those passing all four, take the one with the smallest "useful" median.
  - **Risk:** the candidates grow from two to three, which slightly increases the chance of passing by chance. Noted in the report.
  - **Scale:** the simulations add about 2739 fits, Salehi adds 67.
