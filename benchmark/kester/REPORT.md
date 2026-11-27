# Kester organoids: KARST 1.1.4 rc against ALFA-K by missegregation multiple

2026-11-26. Plan: `benchmark/rc_1.1.4/PLAN.md`, part 3 (options D1a, D4a).
Prompt: `dev/KESTER_PREDICTION_PROMPT_CLAUDE.md`.

Tables: `results/summary/` (local) and
`/share/lab_crd/taoli/Project/Kester/results/forecast_50e5889/` (RED). The RED
folder holds:
- inputs, the fold manifest and the task table;
- every fit and ABM forecast;
- scores and compositions;
- versions and seeds.

## Setup

**Versions:** KARST 50e5889 (1.1.4 rc, `lib/karst_50e5889`); alfakR 0.0.1.0
(dev2 f93b2d2 image).

**Time axis:** one passage = 7 days, from the weekly passaging. Both methods
get the same time axis (passages rebased to the first training passage,
`dt` = 7). ALFA-K is fitted as in the paper benchmark (`passage_times = NULL`),
not with `dt` = 1 and passage times as the prompt suggests; the time axis is
the same.

**Rates:** missegregation 5e-5 × {1, 2, 4, 8, 16, 32, 64}, the same for
fitting and forecasting. n0 / nb are 1e5 / 1e7 (Salehi's); minobs 5, 10 and
20.

**Fits:**
- 476 fit tasks, 14 CPU-hours.
- The only failures are minobs 20 (no karyotype with ≥ 20 cells) in four
  GBM-like folds, for both methods.
- 378 ABM forecasts.

**Folds (`fold_manifest.tsv`):**
- Salehi-like: APKS1–3, 6 folds on all training passages, plus 3 on the
  endpoints only;
- GBM-like: A3 and AP1, seeded three-fold cell holdout within each passage
  (seed 20261126), with kflInfer's partitions and scoring;
- held-out target cells outside list E: 0 in five folds, 0.8% in one.

**Data are tiny:**
- 23–541 cells per passage;
- 3–20 distinct training karyotypes per fold;
- AP1 has 10 held-out cells per fold.

All results are descriptive.

## Salehi-like folds (APKS1–3)

**Per-cell log score on the held-out passage** (engine on list E, the same
start and horizon for every arm). Means over the 6 all-passage folds:

| Multiple | KARST | ALFA-K m10 | Neutral | Persistence |
|---|---|---|---|---|
| ×1 | −1.900 | −2.963 | −1.888 | −2.111 |
| ×4 | −1.885 | −2.962 | −2.140 | −2.111 |
| ×8 | −1.880 | −3.001 | −2.503 | −2.111 |
| ×16 | −1.887 | −3.129 | −3.106 | −2.111 |
| ×64 | −2.043 | −3.661 | −4.804 | −2.111 |

**Paired differences per cell** (folds where the first arm is better, of 6):

| Comparison | ×1 | ×8 | ×64 |
|---|---|---|---|
| KARST − ALFA-K m10 | +1.06 (6/6) | +1.12 (6/6) | +1.62 (6/6) |
| KARST − neutral | −0.01 (4/6) | +0.62 (6/6) | +2.76 (6/6) |
| KARST − persistence | +0.21 (5/6) | +0.23 (6/6) | +0.07 (2/6) |

- KARST beats persistence. It ties with neutral at ×1, as on Salehi.
- It beats neutral at higher multiples, because a flat landscape spreads into
  new karyotypes faster than the observed composition moves.
- KARST's own score changes little with the multiple (best at ×8).

**Fig. 2c-style angle** (mean over folds):

| Design | Forecast | KARST | ALFA-K |
|---|---|---|---|
| All passages | engine | 60°–66° | 91°–106° |
| All passages | paper ABM, own landscape | 62°–65° | 77°–104° |
| Endpoints only | paper ABM | 78° (×1) to 101° (×64) | m5: 52°–74° |

- On the all-passage folds, KARST's angle is much smaller than ALFA-K's at
  every multiple, with both forecast types.
- On the endpoint folds, ALFA-K m5 has the smaller angle.

**All passages against the endpoints** (3 folds, same targets):
- KARST's log score is the same within 0.02 at every multiple.
- ALFA-K's endpoint fits score about 3 nats per cell worse than its
  all-passage fits.

## GBM-like folds (A3, AP1): within-harvest cell holdout scores

Negative log likelihood per held-out cell at the later passage (lower is
better), mean of 3 folds; kflInfer's scoring (half-count pseudocount; all
held-out cells scored):

| Culture | Multiple | KARST native expected | KARST engine | ALFA-K m10 engine | Neutral | Persistence |
|---|---|---|---|---|---|---|
| A3 | ×1 | **1.163** | 1.174 | 1.213 | 1.432 | 1.223 |
| A3 | ×8 | 1.210 | 1.224 | 1.341 | 2.608 | 1.223 |
| AP1 | ×1 | 0.783 | 0.842 | 5.92 | 2.133 | 2.396 |
| AP1 | ×8 | **0.706** | 0.722 | 5.91 | 3.453 | 2.396 |

- On both cultures, KARST is best.
- A3 has 3–4 karyotypes, so all arms are close.

## Why ALFA-K's engine scores are low on the smallest inputs

ALFA-K's fitted rates on these tiny inputs are unstable in scale (×1,
`diag_alfak.R`):

| Input | ALFA-K rate range (per day) |
|---|---|
| AP1 | 15 to 100 |
| APKS2 endpoints | −4 to 20 |
| APKS3 (all and endpoints) | 0.01 to 0.2 |

KARST's rates lie between 0.57 and 1.03 per day on all of them.

With such spreads, the deterministic engine moves almost all cells onto the
karyotype with the highest rate, often a Kriged neighbour. The log score then
collapses (AP1, the endpoint folds). ALFA-K's paper-style ABM forecast on its
own landscape keeps a good angle on the endpoint folds.

These ALFA-K numbers reflect its fits on very few cells, not a comparison of
landscape shapes at equal scale.

## Reading

- On the Kester trajectories with several passages, and on the two-passage
  cell holdouts, KARST forecasts the held-out cells better than ALFA-K and
  persistence, and at least as well as neutral.
- The missegregation multiple changes KARST's scores only slightly (best
  around ×8). The data are too small to separate the rates.
- With 6 + 3 folds and 2 cultures, nothing here is a test with a verdict.
