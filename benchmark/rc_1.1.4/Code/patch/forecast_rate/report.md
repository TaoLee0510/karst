# B: using an externally given missegregation rate: results

2026-11-25. Plan and pre-registered decision rules: `PLAN.md`. KARST was not
modified; β is not estimated from the data anywhere here.

## Verdict

- **Q1, by the pre-registered rule: a correct rate does not improve W1.**
  - Forecast only: 3.3% worse.
  - Fit and forecast: 2.3% better.

  Neither reaches the 0.95 threshold.
- **Q4: the mixture failed** on the development half. It costs 3–37% in W1 on
  ×1 data, against an allowed 2%.
- **What a correct rate does change:**
  - the share of new karyotypes and the near out-of-domain direction improve
    strongly with forecast-only;
  - on Salehi, a high rate used for both fitting and forecasting improves
    the angle strongly (descriptive, in-sample).

## Simulations

837 inputs (18 varyMisseg inputs have no ×1 fit, as in the round); all fits
match their rebuilt inputs. Results: `results/final/`, `results/final_log.txt`.

**Q1 (varyMisseg, true rate ≥ ×4; 207 inputs, 75 simulations).** Median
paired ratios, 95% interval by simulation:

| Comparison | W1 | Angle change | New-karyotype error | Near out-of-domain error |
|---|---|---|---|---|
| Forecast only (f ×1, g true) vs base (×1, ×1) | 1.033 (1.026, 1.041) | +0.1° | 0.42 | 0.95 |
| Fit and forecast (f = g = true) vs forecast only | 0.953 (0.937, 0.964) | −1.3° | 1.42 | 1.04 |
| Fit and forecast vs base | 0.977 (0.961, 0.990) | −0.7° | 0.59 | 0.96 |

**Q2 (by true rate; direction = share of inputs predicting a rise where the
near out-of-domain share rose):**

| True rate | Base direction | Forecast-only direction | Fit-and-forecast direction | True landscape at ×1: W1 / angle |
|---|---|---|---|---|
| ×2 | 65% | 87% | 65% | 0.82 / −13° |
| ×4 | 40% | 93% | 46% | 0.88 / −15° |
| ×8 | 56% | 99% | 53% | 0.93 / −33° |
| ×16 | 29% | 34% | 49% | 0.90 / −50° |

- The true landscape forecast at ×1 is far better than KARST at any rate. So
  W1 and the angle are set mainly by fitness, not by the rate.
- At the true rate, the true landscape improves W1 only a little further
  (×16: 0.887 against 0.899), but it fixes the direction (×16: 91% against
  54%).
- With β fixed at 5e-5, the observed share of new karyotypes reaches 0.29 at
  ×4 and 0.77 at ×16, while KARST forecasts 0.05–0.07. At the true rate, the
  true landscape itself forecasts only 0.10 and 0.39: offspring outside list
  E are dropped.

**Q3 (misspecification, median W1 ratio against the true rate):**
- On ×1 data with forecast only:
  - main ABM: ×2 1.05, ×4 1.20, ×8 1.61, ×16 2.40;
  - varyMisseg ×1: 1.02 / 1.13 / 1.45 / 2.11.
- Fitting and forecasting at too high a rate costs less: on varyMisseg ×1,
  1.02 / 1.07 / 1.12 / 1.16.
- Rates that are too low do not hurt W1 (0.96–0.98 with forecast only).
  W1 favours low forecast rates; the new-karyotype error favours the true
  one.

**Q4 (mixtures from f = ×1, development half):**

| Range | W1 gain at ≥ ×4 | Cost on main ×1 | Cost on varyMisseg ×1 |
|---|---|---|---|
| ×1–×4 | 1.00 | 1.07 | 1.03 |
| ×1–×8 | 1.01 | 1.17 | 1.09 |
| ×1–×16 | 1.03 | 1.37 | 1.22 |

No range qualifies, so the confirmation half was not judged.

## Salehi (descriptive; in-sample)

- **Coverage:** 34 windows. SA000_X6 has no ends fit and is left out, as in
  `../outside_rank`; 41 daughters.
- **Check:** in every window, the ×1 fit reproduces codex's karst_native map
  exactly.
- **Rerun:** one window (SA609_X7) overflowed under a high fitted rate and
  was rerun with rescaling added to this study's engine copy (shares
  unchanged; `R/engine.R`).

**Medians over the 41 daughters:**

| Arm | Angle | W1 | Overlap | Direction right | New-karyotype error |
|---|---|---|---|---|---|
| Current (fit ×1, forecast ×1) | 88.6° | 1.46 | 0.249 | 10% | 0.45 |
| Forecast only ×8 | 89.0° | 1.61 | 0.228 | 87% | 0.16 |
| Forecast only ×16 | 89.6° | 1.65 | 0.172 | 90% | 0.17 |
| Fit and forecast ×8 | 80.8° | 1.47 | 0.287 | 23% | 0.23 |
| Fit and forecast ×32 | 74.2° | 1.40 | 0.311 | 29% | 0.11 |
| Fit and forecast ×64 | **52.8°** | **1.40** | 0.296 | 35% | 0.08 |
| Mixture ×1–×16 | 87.2° | 1.54 | 0.231 | 84% | 0.23 |
| ALFA-K (×1) | 68.2° | 1.67 | 0.239 | 48% | 0.41 |

- Fit and forecast at ×64 against current:
  - the angle is smaller in 28 of 41 daughters;
  - W1 is smaller in only 18.
- **By family** (median angle, current → fit and forecast ×64; ALFA-K):

  | Family | Current | ×64 | ALFA-K |
  |---|---|---|---|
  | SA609 | 59.5° | 48.9° | 68.2° |
  | other | 88.6° | 69.9° | 68.3° |
  | SA906 | 94.9° | 72.3° | 47.9° |

  SA906's W1 worsens (1.18 → 1.37).
- **Angle wins against ALFA-K:**
  - current: 22;
  - fit and forecast ×4 and ×64: 24;
  - forecast only ×16: 18.

## Reading

1. **If a rate is given from outside, use it for both fitting and
   forecasting, not only for forecasting.**
   - In the simulations, forecast-only makes W1 worse (1.033), while both
     improve it slightly (0.977).
   - On Salehi, forecast-only spreads cells into new karyotypes and worsens
     the angle and W1, while both improve the angle (to 52.8° at ×64).
2. **The gain from a correct rate is small on the pre-registered metric**
   (W1 2.3% in the simulations). The large gains are on the new-karyotype
   share and the near out-of-domain direction, which W1 does not reward.
3. **Overshooting costs more than undershooting.** That is the same
   asymmetry as before, and the reason a rate should come from outside
   rather than be guessed.
4. **A fixed mixture over rates is not a safe default:** it costs 3–37% in
   W1 on ×1 data.
5. **Fitness is the main lever:** the true landscape at ×1 is 7–20% better in
   W1 and 7–50° better in angle than KARST at ×1.

## Cost and files

- **CPU:** about 100 CPU-hours (simulations 89.5, Salehi 9.8, pilots
  about 0.5).
- **Total with A′:** about 411 CPU-hours, past the 300 line; the user was
  told and released the Salehi part.
- **Code:**
  - `R/fr.R`, `R/engine.R` (with the rescaling), `R/common.R`;
  - `scripts/sim_task.R`, `scripts/salehi_task.R`, `scripts/analyse.R`;
  - `hpc/`.
- **Results:**
  - `results/final/`: sim_rows, q1–q4, salehi_rows, verdict;
  - `results/out/`: per-input outputs;
  - on the HPC: `/share/lab_crd/taoli/Project/karst_patch/forecast_rate/`.

## Addendum (2026-11-26): ALFA-K arms in the simulations (descriptive)

The user asked for this after the report. `scripts/alfak_task.R` builds the
ALFA-K rate map on list E:
- ALFA-K's minobs-10 landscape on its own list (always inside E);
- its stable Kriging prediction elsewhere, as `benchmark/paper/R/align.R` and
  codex's Salehi maps do.

The engine, parent, daughter and scores are the same as for the KARST arms.

**Coverage and arms:** 838 inputs. The 17 misseg inputs without fits are the
same as for KARST, and 8 main-ABM inputs at λ 0.2 have an undefined ALFA-K
angle because the centroid does not move. The arms are:
- `alfak_g1`: the β-fixed fit, forecast at ×1;
- `alfak_gtrue`: the same fit, forecast at the true rate;
- `alfakT_gtrue`: the fit at the true β, forecast at the true rate.

**Cost:** 36.6 CPU-hours. Results: `results/final/karst_vs_alfak.tsv`.

**KARST against ALFA-K.** The W1 column is the median paired ratio
KARST / ALFA-K, with a 95% interval by simulation:

| True rate | Both default: W1 | KARST W1 smaller | Angle change | Both fitted and forecast at the true rate: W1 | Angle change |
|---|---|---|---|---|---|
| main ABM (×1) | 0.718 (0.685, 0.768) | 80% | −4.9° | — | — |
| varyMisseg ×1 | 0.662 (0.541, 0.823) | 88% | −11.8° | — | — |
| ×2 | 0.706 (0.638, 0.807) | 88% | −14.1° | 0.708 (0.651, 0.810) | −14.0° |
| ×4 | 0.801 (0.716, 0.890) | 91% | −6.1° | 0.819 (0.762, 0.873) | −7.1° |
| ×8 | 0.900 (0.817, 0.926) | 88% | −4.7° | 0.915 (0.864, 0.945) | −5.9° |
| ×16 | 0.946 (0.896, 0.990) | 72% | −0.5° | 0.929 (0.859, 0.971) | −4.7° |

- KARST's forecasts are closer than ALFA-K's at every true rate.
- The lead shrinks as the population gets more dispersed. At ×16 with both
  at their defaults, KARST's angle is smaller in only half of the inputs.
- With the true rate used for both fitting and forecasting, KARST's lead at
  ×16 widens in angle (−4.7°) and stays at 7% in W1.
- KARST at the true rate against ALFA-K at its default: ×16 0.872, ×8 0.879.
- The true landscape against ALFA-K: 0.52–0.82.

**ALFA-K predicts a rise of the near out-of-domain share more often than
KARST.** Share of inputs where the share rose and the arm predicts a rise:

| Arm | Main ABM | varyMisseg ×1 | ×4 | ×16 |
|---|---|---|---|---|
| ALFA-K | 89% | 89% | 84% | 54% |
| KARST | 76% | 75% | 40% | 29% |

- ALFA-K's Kriging extrapolates a smooth surface, which keeps out-of-domain
  karyotypes near the observed levels.
- KARST's trend declines with distance from the majors.
- This does not translate into a better W1 or angle for ALFA-K.

**Difference from the Fig. 2d result:** in fig2d_misseg, with β fixed,
KARST's Fig. 2d (rank of the fitness change) fell below ALFA-K's from 2e-4
(×4) onward. On the full-distribution forecast (W1, angle), KARST stays ahead
at every rate here. The two tests measure different things.
