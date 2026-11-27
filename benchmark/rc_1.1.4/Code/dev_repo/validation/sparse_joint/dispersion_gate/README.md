# Dispersion gate for the missegregation rate

Plan and pre-registered criteria: `PLAN.md` (approved 2026-11-23, D1a–D5a).
Status: **stopped after step 1, which failed** (2026-11-23 19:45 EST).
Step 2 was not run. KARST keeps β = 5e-5.

## Step 0: dispersion across sources (passed, `once` only)

`step0_dispersion.R`, 17,600 existing datasets, no fits;
results in `results/step0/` (`s0.tsv`, `cells.tsv`, `step0.pdf`).

| statistic | n | τ | ×16 above τ | ×4 above τ |
| --- | --- | --- | --- | --- |
| once | 100 | 0.602 | 62% | 29% |
| once | 300 | 0.402 | 67% | 32% |
| new | 100, 300 | 1 | 0% | 0% |

- `new` fails: in the e2 sweeps and with long intervals, t1 shares no
  karyotype with t0 even at ×1.
- `once` passes, but only through smooth landscapes:
  - at λ 0.2 and 0.4, ×16 is almost never above τ (0% outside e3);
  - τ is set by the ×1 datasets at λ 1.6 (own benchmark, e3).
- The source matters: at ×1, the median of `once` (n = 100) is 0.12 in the
  paper's ABM at λ 1.6 and 0.41 in the own benchmark at λ 1.6.

## Step 1: calibration (failed)

Scripts:
- `step1_fit_v7.R` and `step1_pick_v7.R`: HPC fits;
- `step1_analyse.R`: the analysis.

Data and fits:
- v7 dev half at 200/200 and 500/700 cells (2,400 datasets), KARST 1.1.3 at
  ×1, ×2, ×4, ×8, ×16. The ×1, ×4 and ×16 fits are bit-identical to
  `dev/improve`'s.
- ABM_varyMisseg (125 sims × 3 designs at the Salehi cell counts), KARST at
  the same multiples.
- 188 v7 datasets end in a domain error at every rate.

Results: `results/step1/analysis/` (`step1.txt`, `k_choice.tsv`,
`closer.tsv`, `pick.tsv`, `gate_cells.tsv`, `step1.pdf`).

**The gate itself behaves as planned.**
- Raised at ×1: at most 4.2% of any development cell (2 datasets, e3 λ 1.6),
  and 0% in varyMisseg.
- Raised at ×16: all datasets at λ ≥ 0.8; varyMisseg 85–89% at ×8/×16.

**No k helps on the raised datasets.** The ratio is the error at ×k over the
error at ×1 (v7: RMSE at distance 1; varyMisseg: 1 − r):

| k | H (×1, n = 2) | G (v7, truth ≥ ×4, n = 582) | v7 ×4 | v7 ×16 | varyMisseg G |
| --- | --- | --- | --- | --- | --- |
| 2 | 2.05 | 0.957 | 1.06 | 0.93 | 0.92 |
| 4 | 4.02 | 1.103 | 1.50 | 0.90 | 0.87 |
| 8 | 6.05 | 1.364 | 2.20 | 0.91 | 0.78 |

- The criterion needs H ≤ 1.25 and G ≤ 0.90. G fails for every k, so the
  verdict does not rest on H, which has only two datasets.
- On the raised datasets, the higher rate beats ×1 in 54% (×2), 47% (×4)
  and 42% (×8) of v7 cases.

**By source** (raised, truth ≥ ×4, k = 4; post hoc):
- C2 (passaging): better, 0.65–0.85;
- e1: mixed, 0.83–2.4;
- e2: mostly worse, 0.75–1.9;
- e3 (neutral): worse, 1.9–4.2, except λ 1.6 at ×16 (0.85);
- varyMisseg: better, 0.85–0.91.

A post hoc look at the correlation error (1 − r at distance 1) on the same
raised datasets:

| | ×2 | ×4 | ×8 | ×16 |
| --- | --- | --- | --- | --- |
| C2 | 0.95 | 0.79 | 0.66 | 0.60 |
| e1 | 0.98 | 0.97 | 0.94 | 0.93 |
| e2 | 1.05 | 1.12 | 1.06 | 1.14 |
| e3 | 0.99 | 0.99 | 0.98 | 0.99 |

So much of the RMSE loss in the e sets comes from the scale of the fitted
rates, not their order. This does not change the verdict: the
pre-registered metric is the RMSE.

**Is a rate closer to the truth better?** (median ratio against ×1, all
development datasets; `closer.tsv`)

| truth | ×2 | ×4 | ×8 | ×16 |
| --- | --- | --- | --- | --- |
| v7 ×1 | 1.63 | 2.94 | 5.05 | 6.98 |
| v7 ×4 | 0.86 | 0.86 | 1.47 | 2.71 |
| v7 ×16 | 0.93 | 0.83 | 0.75 | 0.71 |
| varyMisseg ×1 | 1.09 | 1.33 | 1.77 | 2.62 |
| varyMisseg ×4 | 0.96 | 0.96 | 1.06 | 1.29 |
| varyMisseg ×8 | 0.91 | 0.85 | 0.80 | 0.86 |
| varyMisseg ×16 | 0.95 | 0.92 | 0.77 | 0.72 |

- On average, yes: the best rate is at or near the truth.
- The gain from the true rate is modest: 14% at ×4 and 29% at ×16 in v7.
- Overshooting costs far more: ×2 on ×1 data is already 1.6 times worse in v7.
- On the dispersed datasets the gate picks, the average gain does not hold
  in the e sets.

**Pick step:** none of the objective, posterior objective and cell-half
hold-out picks the more accurate fit better than chance (46–53% right, 70%
needed). On varyMisseg (reported only), KARST's CV score picks right in 70%
and the objective in 73%.

**Cost:** 474 CPU-h (v7 fits 369, varyMisseg 71, pick 33), over the
150–250 estimated for step 1 because each task carried about a minute of
image and scoring overhead. Total within the 650 CPU-h budget.

## Decision

The run stopped as the plan requires. β stays 5e-5; no package change.
What to do next is the user's decision.
