# D4 held-out evaluation: F2 did not pass

Written 2026-11-20 01:17 EST. Script `d4.R judge`, output `dev/results/realdata/CE/d4_judge.log`. Evaluated only this once.

**Procedure:** Following user decision 1 (PLAN.md, 01:15), the development threshold was relaxed to 0.945. On re-selection, F2 qualified at quantile 0.979: development GRF coverage 0.947 / 0.961 / 0.948 / 0.954, Salehi development-group parallel lineages 0.066. F1 and F3 still did not qualify. The decision criteria on the held-out data were not changed.

## Conclusion: did not pass (D4-1, D4-2, D4-4 not met)

### Coverage on simulations

Each category is "first averaged within the category for each dataset, then averaged over datasets".

| Data | Form | Observed < 5 | Observed ≥ 5 | Unobserved d1 | Unobserved d2 | Fraction of intervals excluding 0 | Of which sign correct | Median width |
|---|---|---|---|---|---|---|---|---|
| Held-out GRF | F2 | **0.896** | 0.955 | 0.914 | 0.932 | 0.41 | 0.999 | 0.259 |
| Held-out GRF | Current | 0.969 | 0.977 | 0.990 | 0.985 | 0 | — | 0.251 |
| RMF | F2 | 0.965 | 0.958 | 0.904 | **0.813** | 0.98 | 0.998 | 0.258 |
| RMF | Current | 0.960 | 0.972 | 1.000 | 1.000 | 0 | — | 0.463 |
| NK (report only) | F2 | 0.577 | 0.893 | 0.519 | 0.661 | 0.98 | 0.999 | 0.333 |
| NK (report only) | Current | 0.704 | 0.899 | 0.785 | 0.755 | 0 | — | 0.709 |
| e1 | F2 | 0.943 | **0.834** | **0.883** | **0.861** | 0.52 | 0.965 | 0.283 |
| e1 | Current | 0.987 | 0.965 | 0.992 | 0.986 | 0 | — | 0.293 |
| e2 | F2 | **0.800** | **0.451** | **0.810** | **0.877** | 0.33 | 0.989 | 0.255 |
| e2 | Current | 0.803 | 0.425 | 0.871 | 0.910 | 0 | — | 0.261 |
| e3 (neutral) | F2 | 0.901 | **0.831** | **0.695** | **0.543** | 0.35 | **0.590** | 0.270 |
| e3 (neutral) | Current | 1.000 | 1.000 | 1.000 | 1.000 | 0 | — | 0.233 |

### Real data

| Check | F2 | Current | Criterion |
|---|---|---|---|
| Salehi confirmation-group parallel lineages \|z\| > 1.96 | **0.135** | 0.035 | ≤ 0.10 |
| GBM / Ovarian agreement between the two halves (10+ cells) | 1 / 0.985 | 1 / 1 | ≥ 0.90 |
| Intervals excluding 0 among karyotypes with ≥ 10 cells (Salehi confirmation group / GBM) | 0.55 / 0.62 | 0 / 0 | median ≥ 0.20 |

### Verdict on each criterion

- **D4-1 did not pass:** held-out GRF, RMF and e1–e3 fall below 0.93 / 0.90 in several places; NK is much lower than with the current constant.
- **D4-2 did not pass:** parallel lineages 0.135 > 0.10.
- **D4-3 passed.**
- **D4-4 did not pass:** the informative fraction meets the criterion, but the sign-correct rate is 0.590 on e3 and 0.965 on e1, below 0.95.

## What we saw

- **The most important case is the neutral data e3.** On e3 the true fitness differences are 0, yet F2 gives 35% of the karyotypes an interval excluding 0, and all of these are false positives. The current form is so wide that it never excludes 0, which is exactly what guards against this.
  - That is, intervals calibrated on GRF and binned by |z| carry over the underestimate of the model sd on other data as well.
  - This agrees with the phase 3 conclusion: the interval shortfall comes from model misspecification, the scale differs between data types, and a constant from one data type cannot be used to calibrate another.
- **On e2 (clone replacement), the current form also reaches only 0.425 on karyotypes observed ≥ 5 times,** comparable to F2. e2 is already the weak spot of the current intervals.
- **On real data, F2 gives more than half of the karyotypes an interval excluding 0.** But the agreement between parallel lineages gets worse with it (0.035 → 0.135), which shows that on real data these intervals are also too narrow.
