# B: how KARST should use an externally given missegregation rate

Approved by the user on 2026-11-25 at about 19:15 EST ("proceed as recommended, confirmed";
options D3a, D4a, D5a of the chat plan). Code here; KARST is not modified;
nothing is committed; HPC results in
`/share/lab_crd/taoli/Project/karst_patch/forecast_rate/`.

## Scope

- This study does **not** estimate β from the data (that line stays closed,
  2026-11-20). It asks what to do with a rate given from outside:
  1. use it for fitting, for forecasting, or both;
  2. how much it buys when it is right;
  3. how much it costs when it is wrong.
- **Mixture arm (D3a):** a forecast mixture over a fixed range of rates,
  which uses no information from the data. If it passes, it is only a
  candidate option; packaging needs the user's separate approval.

## Background (rechecked before planning)

- Every attempt to estimate β from the data failed (alfak2 β auto; E1–E3;
  noise v1–v7; the dispersion gate).
- With the right rate, fitted-rate accuracy improves modestly on varyMisseg
  (1 − r ratio 0.72–0.80 at ×8/×16) and v7 (0.86, 0.71).
- Overshooting costs a lot (×2 on ×1 data: 1.09 on varyMisseg, 1.63 on v7);
  the e sets lose.
- On Salehi, changing only the forecast rate raised the later-passage log
  score by +0.45 per cell, while changing only the fitted β gave +0.03.
- The effect of the forecast rate alone on forecasts (angle, W1, new
  karyotypes) has never been measured on simulations with known truth.

## Data

**Simulations:**
- varyMisseg: 125 simulations (true rate 5e-5 × {1, 2, 4, 8, 16}, 25 each)
  × 3 designs (w2 / w4 / w8 ends) at Salehi cell counts.
  - KARST 1.1.3 fits at f = ×1 come from the round (fig2d_misseg beta_fixed);
    f = ×2, ×4, ×8 and ×16 come from the dispersion-gate study
    (`validation/sparse_joint/dispersion_gate/results/step1/vm`, read only).
  - Checked on one input: same counts as ×1, and the fit at the true rate is
    bit-identical to fig2d_misseg beta_true.
- Main ABM (true rate ×1): replicates 1–40 of each wavelength × 3 designs (480
  inputs), round fits (fig2d_truth step 2).
- **Halves:** odd replicates form the development half and even replicates
  the confirmation half (as in `../outside_rank`).

**Forecast (as `../outside_domain`):**
- start from the fitted parent;
- list E (the training karyotypes and their one- and two-step neighbours);
- the deterministic engine with a 0.1-day step, to the second sampling time
  after the parent;
- the observed value is the simulation's 1,000-cell sample.

## Arms (`scripts/sim_task.R`)

- **KARST native division rate** of a fit at f:
  `fitness_mean / (2 P_ii(f) − 1)` from `predict()` on E. The engine is
  then run at the forecast rate g (multiples of 5e-5).
- `f1_g{1,2,4,8,16}`: fit at ×1, forecast at g.
- `f{m}_g{m}`, m = 2–16 (varyMisseg only): fit and forecast at the same rate.
- `f1a_g*`: the same as `f1_g*`, level-aligned to the truth (report only).
- `mix_1_{4,8,16}`: equal-weight mixture of the `f1_g` forecasts with g from
  ×1 up to ×4, ×8 or ×16.
- `truth_g1`, `truth_gtrue`: the true landscape at ×1 and at the true rate.

**Metrics:**
- angle of the centroid change;
- W1 of the full distribution (floor 1e-5);
- exact-karyotype overlap;
- change of the near out-of-domain share (≤ 4 steps from a major): its
  error and its direction;
- the share of new karyotypes (not seen in training): error of the forecast
  share;
- the daughter's share on list E.

## Questions and decision rules (fixed before any result)

Paired per input; medians; 95% intervals by resampling whole simulations
(2,000 draws).

**Q1, where to use the rate.** On varyMisseg inputs with true rate ≥ ×4,
using all replicates:
- FO = `f1_g{true}` against base = `f1_g1` (median W1 ratio m_FO);
- BOTH = `f{true}_g{true}` against FO (m_B), and against base.

| Condition | Recommendation |
|---|---|
| m_FO ≤ 0.95 and m_B ≥ 0.98 | forecast only |
| BOTH against base ≤ 0.95 and m_B < 0.98 | fit and forecast |
| neither FO nor BOTH ≤ 0.95 against base | a correct rate does not improve W1 |

The angle and the new-karyotype error are reported for the same comparisons.

**Q2, the gain when right:** for each true rate, the ratios of the
recommended usage against base (W1, angle, new-karyotype error, halo error),
and the share of inputs where the near out-of-domain direction is right.

**Q3, the cost when wrong:**
- for each true rate t, `f1_g` against `f1_g{t}`, and `f{g}_g{g}` against
  `f{t}_g{t}`, for g = t/4, t/2, 2t, 4t on the grid;
- on ×1 data (main ABM and varyMisseg 5e-5), `f1_g{2,4,8,16}` against
  `f1_g1`.

**Q4, the mixture:**
- **On the development half,** choose the range among ×1–×4, ×1–×8 and
  ×1–×16 with the smallest median W1 ratio against base at true rate ≥ ×4,
  subject to a median W1 ratio ≤ 1.02 at ×1. ×1 is judged separately on main
  ABM and on varyMisseg 5e-5, and both must hold. If no range qualifies,
  Q4 fails.
- **On the confirmation half, the chosen range passes** if:
  - the median W1 ratio is ≤ 1.02 at ×1, on main ABM and on varyMisseg 5e-5
    separately;
  - the median W1 ratio is ≤ 0.95 at ≥ ×4, with the 95% interval below 1;
  - the median angle change at ×1 is ≤ +1°.

## Salehi (D4a, descriptive, no verdict; `scripts/salehi_task.R`)

- **Data:**
  - codex's 35 windows (shared full parent, list E, the engine of
    `../outside_rank` step 3);
  - `f1_g` for g = ×1 to ×64 from the round fits;
  - `f{g}_g{g}` from the β-grid fits (×2 to ×64);
  - the three mixtures, and ALFA-K at ×1.
- **Check:** the round fit at ×1 must reproduce codex's karst_native map
  (largest relative difference ≤ 1e-6); windows that do not are left out.
- **Report:**
  - median angle, W1, overlap, near out-of-domain direction and
    new-karyotype error, by family;
  - wins against ALFA-K.
- The best rate is labelled in-sample and optimistic. No rate is chosen per
  line, because that would be estimating it.

## Running

- Runs automatically to the end (D5a) and reports.
- **HPC:**
  - one input or window per array task;
  - qos xxlarge, the excluded-node list, and an add-hours monitor per job.
- Shared budget with A′: about 200 CPU-hours; ask the user before the total
  goes over 300.
