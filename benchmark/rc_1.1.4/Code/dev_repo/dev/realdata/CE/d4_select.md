# D4 selection stage: no qualifying form, D4 failed

Written 2026-11-20 00:47 EST. Script `d4.R select`, output `dev/results/realdata/CE/d4_select.log`.

**Data used:** development data only, i.e. the 96 datasets of the dev split of the development GRF (cd342db base fits, 361 141 karyotypes), the Salehi development group and Ovarian. Nothing was computed on the held-out data.

| Form | Quantile | Coverage: observed < 5 | Coverage: observed ≥ 5 | Coverage: unobserved d1 | Coverage: unobserved d2 | Fraction excluding 0 in development simulations | Of which sign correct | Parallel lineages \|z\|>1.96 (Salehi development) | Excluding 0 among karyotypes with ≥10 cells (Salehi development) | Same as left (Ovarian) | Qualified |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Current constants | — | 0.982 | 0.972 | 0.991 | 0.984 | 0 | — | 0.023 | 0 | 0 | Not a candidate |
| F1 refit (floor 0.0069, deficit 0.33) | — | 0.958 | 0.950 | 0.960 | **0.939** | 0.002 | 1.000 | 0.060 | 0.31 | 0.38 | No (d2 < 0.95) |
| F2 two sides calibrated separately | 0.98 (upper limit) | **0.949** | 0.967 | 0.950 | 0.956 | 0.306 | 0.999 | 0.063 | 0.43 | 0.55 | No (observed < 5 short by 0.001) |
| F3 symmetric conformal | 0.95 | 0.950 | 0.951 | 0.950 | 0.950 | 0.133 | 1.000 | **0.081** | 0.33 | 0.38 | No (lineage > 0.07) |

Per the PLAN and clarification C1, only forms that meet "development GRF coverage ≥ 0.95 in every class, and Salehi development group parallel lineages ≤ 0.07" are candidates. None of the three meets this, so D4 does not enter the held-out evaluation and is judged failed.

## What was seen (does not change the verdict)

- **The current intervals are wide because of the deficit.** After refitting the deficit from 0.6 to 0.33, the intervals can exclude 0: in the Salehi development group 31% of the intervals of karyotypes with ≥10 cells exclude 0, and in the development simulations the signs of such intervals are all correct. F1 falls somewhat short on coverage only for unobserved karyotypes at distance 2 (0.939).
- **F2 is the most informative, short by only 0.001.**
  - The interval is narrow on the side towards 0 and wide on the side away from 0, consistent with an error structure of "estimates slightly shrunk".
  - In the Salehi development group 43%, and in Ovarian 55%, of the intervals of karyotypes with ≥10 cells exclude 0; in the development simulations the sign accuracy of such intervals is 0.999.
  - The shortfall is in the "observed, < 5 cells, |z| ≤ 1" bin: the quantiles towards 0 / away from 0 are 5.7 / 7.7 times u, so the error tails are very heavy.
- All three new forms are at 0.06–0.08 on the Salehi development group parallel lineages, higher than the 0.023 of the current constants, but far below the 0.31 of the model sd.

## For the user to decide

1. Accept that D4 failed.
2. Relax the threshold and judge once more (for example development coverage ≥ 0.945, or a quantile upper limit of 0.985). This is a change made after seeing the development results, so you need to decide. After the change, the held-out data are still evaluated only once.
