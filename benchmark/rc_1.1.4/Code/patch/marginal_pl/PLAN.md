# Per-chromosome marginal pseudo-likelihood: plan and pre-registered criteria

Approved by the user on 2026-11-26 at about 00:25 EST ("proceed as recommended"; D1a,
D2a, D3a, D4a, D5a). Code here; KARST is not modified; nothing is committed;
HPC results in `/share/lab_crd/taoli/Project/karst_patch/marginal_pl/`.

## Why

On Salehi, single karyotypes carry almost no count-change information (SA609:
no karyotype has ≥ 3 cells at both time points), but every chromosome's
copy-number distribution uses all cells. On SA609, 11 of 22 chromosomes shift
significantly between the two training time points, and 67% of the cells are
outside the domain, so they never enter KARST's likelihood
(`../forecast_rate/diag/info_check.R`).

## Method (`R/mpl.R`)

**The pseudo-likelihood:**

PL(δ) = Σᵢ Σ_c y1[i, c] log p̂1[i, c](δ)

- y1[i, c] is the number of cells at the last time point with c copies of
  chromosome i.
- p̂1 is the same distribution after every karyotype observed at the first
  time point is grown for dt at net rate (2P_ii − 1)·exp(log birth).
- p̂1 is mixed with a uniform over the copy states seen at t1 with weight
  1e-3. Those states cannot arise without missegregation (D3a).

**Update:**
- Maximise −½ δ′ Σ⁻¹ δ + w·PL(δ), where δ is the change of the trend's
  dosage, curvature and shell coefficients, and Σ is the fit's covariance of
  those coefficients.
- The intercept is fixed: a common shift of fitness leaves p̂1 unchanged.
- **D1a:** in the propagation and in the predictions, karyotypes in the
  domain, and unseen karyotypes, get log birth + Δtrend(k). Karyotypes in O
  keep KARST's common level. The reference (majors weighted by their cells)
  is recomputed, and the sd columns are unchanged.
- **D1b (report only):** O also get Δtrend.
- **D2a:** w ∈ {1/22, 1/8, 1/3, 1}, chosen on the development half.

## Data

- **Own benchmark:**
  - development: dev split (168 datasets, refit with KARST 1.1.3 as the
    round);
  - **D4a, confirmation:** a random half of the held-out landscapes (250 of
    500, `sample()` with seed 20261126; all their datasets, 3,500), refit,
    with the identity check against the stored predictions. The chat plan said
    "175 of 350"; the held-out split has 500 landscapes of 14 datasets each
    (corrected before any result).
- **Paper simulations:** main ABM (fig2d_truth step 2) and varyMisseg
  (fig2d_misseg, β fixed and β true) at Salehi cell counts, from the round's
  stored fits. Odd replicates form the development half and even replicates
  the confirmation half.

**Accuracy sets (truth: the true net growth rate):**
- own: the round's scored nodes;
- paper: O, up to 3,000 sampled unseen out-of-domain karyotypes within two
  steps of the observed ones, and the observed in-domain non-majors.

**Metrics:**
- `sp_ou`: Spearman correlation over O ∪ unseen out-of-domain karyotypes (the
  K1 and gate metric);
- `sp_un`: the same over unseen out-of-domain karyotypes only;
- `sp_in`, `rmse_in`: over observed in-domain non-majors, the RMSE after
  centring the prediction and the truth on their own means over the set.

## Step 0: information (development half; no update)

PL per t1 cell, for the truth, KARST's current prediction and neutral (all
equal), for the karyotypes of t0:
- **S0a:** median of PL(truth) − PL(neutral) > 0 in both families.
- **S0b:** median of PL(truth) − PL(KARST) > 0 in both families, and at least
  20% of the median of PL(truth) − PL(neutral).

If either fails, stop. Salehi's PL(KARST) − PL(neutral) is reported,
descriptively.

## Step 1: development

**Gain:** per dataset, `sp_ou`(update) − `sp_ou`(current).
- Clusters are weighted equally, as in `../outside_rank` step 1: own by
  landscape_id, paper by (source, sim_id).
- Only datasets with ≥ 5 karyotypes in O ∪ unseen enter.

**Choice of w:**
- the w with the largest smaller-family gain;
- among those within 0.005 of it, the smallest w.

**Gate:**
- gain ≥ 0.05 in both families;
- cluster-weighted mean `rmse_in` ratio (update / current) ≤ 1.02 in both
  families.

If the gate is not passed, stop.

## Step 2: confirmation

- **K1:** `sp_ou` gain ≥ 0.05 with a 95% interval above 0, in both families
  (2,000 whole-cluster draws, `R/boot.R`).
- **K2 (own held-out half, landscape-weighted as in `../outside_rank`):**
  - d0, d1 and d2 RMSE ratios ≤ 1.02;
  - d0 Spearman not lower;
  - d0 95% coverage lower by ≤ 0.02.
- **K3 (paper even half):** deterministic forecasts as `../outside_rank`
  step 2 (aligned arms), against current:
  - median W1 ratio ≤ 1.01;
  - median near out-of-domain error ratio ≤ 1.01;
  - median angle change ≤ +1°.

  All three hold over all inputs and over the inputs whose forecast changed.
  The paired comparison with ALFA-K (`../forecast_rate` alfak_g1, where
  available) is reported.
- **Pass:** all three.

## Step 3: Salehi (descriptive, no verdict)

On the 34 windows and 41 daughters of `../outside_rank` step 3, with the
update at the chosen w:
- angle, W1, overlap, near out-of-domain direction and new-karyotype error;
- compared with current KARST and ALFA-K.

## Running

- D5a: steps 0–3 run automatically and stop at a stop branch; then report.
- **HPC:**
  - own refits: one fit per array task;
  - stored fits: one analysis job;
  - qos xxlarge, the excluded-node list, and an add-hours monitor per job.
- **Budget:** about 130 CPU-hours; ask the user before 200.
- At most 3 heavy local processes.
