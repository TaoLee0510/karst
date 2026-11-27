# Pre-fit preprocessing v3: retest on new simulated data against the corrected criteria

**Basis:**
- The plan is in `PLAN_v3.md` (fixed 2026-11-20 13:06; the seed base was corrected at 13:07, before any data were generated).
- Method and decision thresholds as in v2; `R/envelope.R` was not changed.
- Numbers come from `v3/v3_eval.R` (criteria) and `v3/v3_extra.R` (report-only breakdowns). Results are in `results/realdata/noise/v3/evaluation.rds`.
- 8160 analyses in total, all successful:
  - 2112 datasets of the e set and 864 of C2 on new landscapes;
  - plus the various injections.

## Conclusion: not passed; stopped as the plan specifies

| Criterion | T = 300 | T = 500 | Requirement |
|---|---|---|---|
| 1. e3, bottleneck 20000, no injection: called sampling-dominated | 6.3% (9/144) passed | 9.4% (9/96) passed | ≤ 10% |
| 2. e3 drift order (bottleneck 200 > 2000 > 20000, gap ≥ 0.30) | 49% > 33% > 6%, passed | 61% > 45% > 9%, passed | in order and gap ≥ 0.30 |
| 3. Simulations with selection (e1, e2, e4, C2): called sampling-dominated | 10.4% (210/2016), **not passed** | 14.6% (196/1344), **not passed** | ≤ 10% |
| 4. Power (e3, bottleneck 20000, φ = 20) | 45% (65/144), **not passed** | 63% (60/96), **not passed** | ≥ 70% |

- Neither threshold passed, so no threshold was enabled.
- Under the rule you set, the preprocessing is not enabled; KARST takes the normal path, with β fixed at the default 5e-5.
- As planned, the Salehi and NUMBAT data are not used, and no P5 plan is written.

**What decides the conclusion is criterion 4 (power), not the marginal gap on criterion 3.**
- Criterion 3 at T = 300 is 10.4%, 95% interval 9.1–11.8%, only slightly off the 10% threshold.
- Criterion 4 is far off: after injecting strong overdispersion with φ = 20, only 10% are called sampling-dominated at 300 cells.

## Why it did not pass: power comes only with more cells, but then real drift is detected as well

### 1. Power rises with the number of cells per time point

Fraction called sampling-dominated for e3, bottleneck 20000, injected φ = 20:

| Cells per time point | 200 | 300 | 500 | 1000 |
|---|---|---|---|---|
| Called sampling-dominated | 0% | 10% | 46% | 79% |

Only at 1000 cells does it exceed 70%.

### 2. At 1000 cells, passage drift at bottleneck 20000 is also detected

Fraction called sampling-dominated for e3, bottleneck 20000, no injection, by the interval between the two time points:

| Interval | 18 days | 42 days | 180 days | 420 days |
|---|---|---|---|---|
| 300 cells | 0% | 0% | 0% | 0% |
| 500 cells | 0% | 0% | 0% | 8% |
| 1000 cells | 0% | 0% | 42% | 25% (another 42% undecidable) |

- At short intervals there are no false calls.
- At intervals of 180 days or more, passage drift accumulates and at 1000 cells is already large enough to be detected.
- So these are not "noise-free data called wrongly": the method detects real drift, and drift is exactly one of the two kinds of noise you defined.
- But criterion 1 treated bottleneck 20000 in advance as "no drift", so by the criterion it counts as a failure.
- With the threshold set at 1000 (report only, not an option in the plan):
  - criterion 4 is 79%, passed;
  - criterion 1 is 17% (8/48), not passed;
  - criterion 3 is 19% (129/672), not passed.

### 3. In simulations with selection, sampling-dominated calls come mainly from strong-drift settings

Fraction called sampling-dominated, no injection:

| Set | Drift | 300 cells | 500 cells | 1000 cells |
|---|---|---|---|---|
| C2 | strong (bottleneck 300) | 1% | 17% | 44% |
| C2 | weak (bottleneck 3000) | 1% | 3% | 11% |
| e1 | strong (bottleneck 200) | 8% | 21% | 29% |
| e1 | weak (bottleneck 2000, 20000) | 2% | 9% | 10% |
| e2 | weak (bottleneck 2000, 20000) | 4% | 8% | 3% |
| e4 | strong (bottleneck 200) | 2% | 17% | 31% |
| e4 | weak (bottleneck 2000, 20000) | 0% | 3% | 8% |

- Looking only at the weak-drift settings, the values under the three thresholds are 5.3%, 7.2% and 8.6%, all within 10%.
- With strong drift (only 200 or 300 cells kept at each passage), selection and drift are both present, and the method reports "sampling-dominated".
  - This also counts as detecting real noise, but the prior criterion 3 does not distinguish strong from weak drift, so by the criterion it counts as a failure.
- In C2, sampling-dominated calls concentrate in settings with a low missegregation rate plus strong drift:
  - β = 1.5e-4 or 5e-4, bottleneck 300, 1000 cells: 65–67%;
  - at β = 1.5e-3, all settings are 0%, almost all called "no clear anomaly", because the true missegregation rate itself falls inside the envelope.

### Taken together

- With 300 or 500 cells per time point, power is insufficient.
- At 1000 cells, power is sufficient, but at long intervals passage drift of normal strength is also reported.
- This test uses only two time points and cannot separate "sequencing sampling noise" from "passage drift". Both are noise, and they look the same at two time points.
- So no threshold can satisfy all four criteria at once.

## Report-only parts

### Region-wise sampling is still often called selection

Cell count ≥ 300, fractions called selection and sampling-dominated:

| Set | Injection | Bottleneck | Called selection | Called sampling-dominated |
|---|---|---|---|---|
| e3 | σ = 0.5 | 200 / 2000 / 20000 | 2% / 9% / 11% | 41% / 31% / 10% |
| e3 | σ = 1.0 | 200 / 2000 / 20000 | 7% / 20% / 25% | 22% / 16% / 9% |
| e1 | σ = 1.0 | 200 / 2000 / 20000 | 6% / 18% / 26% | 9% / 3% / 0% |
| C2 | σ = 1.0 | 300 / 3000 | 23% / 16% | 12% / 4% |

- As in v2: the weaker the drift, the more easily regional noise is taken for selection.
- When the clone families in one region are over- or under-sampled, it looks "patchy", just like real selection.

### The ability to recognise real selection is also low

Fraction called selection with no injection and 1000 cells:
- e1: 17%;
- e2 (clonal replacement): 51%;
- e4: 17%;
- C2: 0–15%.

The neutral e3 also has 5% called selection at 1000 cells.

### With 200 cells per time point, almost everything is called "no clear anomaly" or "undecidable"

| Set | Injection | No clear anomaly | Undecidable | Selection | Sampling-dominated |
|---|---|---|---|---|---|
| e1 | none | 66% | 34% | 0% | 0% |
| e2 | none | 38% | 62% | 0% | 1% |
| e3 | none | 69% | 31% | 0% | 0% |
| e4 | none | 78% | 22% | 0% | 0% |
| e3 | φ = 20 | 24% | 76% | 0% | 0% |

This matches the rule "if undecidable, take the default path", but it also shows that the preprocessing has no effect at this cell count.

### Other

- φ = 5 (variance inflated 5-fold, weaker than φ = 20):
  - at bottleneck 20000, the fractions called sampling-dominated at 300 / 500 / 1000 cells are 15% / 52% / 67%;
  - at 300 and 500 cells this is actually slightly higher than for φ = 20 (with φ = 20 more fall into "undecidable"). I did not investigate the reason.
- Full table split by setting: the last table of `Rscript dev/realdata/noise/v3/v3_eval.R`.

## Relation to v2

The v3 results agree with v2, with the same two limitations:
- power is sufficient only at about 1000 cells;
- regional noise cannot be separated.

v3 fixed the problem with v2's criterion 1 by taking bottleneck 20000 out on its own. But at 1000 cells and long intervals, the drift at bottleneck 20000 is also detectable, so the same problem appears in a different place.

## Next steps (waiting for your decision; I have not done these)

The plan stops here. Directions that could be considered:

1. **Accept the status quo:** the preprocessing does not go into KARST; two-time-point data all take the normal path, with β fixed at the default. This is consistent with the rule you set.
2. **Change the goal:** no longer distinguish "drift" from "sampling"; rename "sampling-dominated" to "noise-dominated", and use it only to alert the user.
   - criteria 1 and 3 would have to be changed to look only at settings where drift is negligible;
   - this needs a new plan and new simulated data; this round's results cannot be re-judged after the fact.
3. **Use more than two time points:** drift accumulates over time while sampling error does not, so with three or more time points the two can in principle be separated. This is an extension built on top of the two-time-point design, not a replacement for it.

## Files

- Code: `dev/realdata/noise/v3/` (`sim_eset.R`, `sim_c2.R`, `v3_task.R`, `v3_eval.R`, `v3_extra.R`).
- HPC results: `results/realdata/noise/v3/{sims,analysis}`.
- Local copy: `dev/results/realdata/noise/v3/analysis`, `evaluation.rds`.
- Jobs:

  | Content | Job id | Tasks |
  |---|---|---|
  | truth | 22169368 | 22 |
  | e-set datasets | 22169370 | 48 |
  | C2 | 22169372 | 72 |
  | analysis | 22169472 | 170 |

  All COMPLETED, no failures.

## Addendum (2026-11-20 15:43 EST): the user pointed out that noise does not distinguish sampling from drift in the first place

User: "Within the noise you can't tell whether it is sampling or drift."

- Under the user's original definition, noise includes both drift and sampling. What the preprocessing has to separate is "noise" from "selection", not "sampling" from "drift".
- So "sampling-dominated" in the calls should be read as "noise-dominated".
- Criteria 1 and 3 counted the detection of real drift as an error; this is a problem with how I set the criteria, the same error as v2's criterion 1, just in a different place.

**Re-examining the v3 data under this reading (post hoc analysis, not a verdict, because the criteria were changed after seeing the results):**

| Question | Result |
|---|---|
| Data with truly almost no noise (e3, bottleneck 20000, interval 18 or 42 days) | 0% called noise-dominated at every cell count |
| Simulations with selection and weak drift (bottleneck ≥ 2000) | 5.3% / 7.2% / 8.6% called noise-dominated at thresholds 300 / 500 / 1000 |
| Simulations with selection and strong drift (bottleneck 200 / 300) | being called noise-dominated no longer counts as wrong (both are present) |
| Power (φ = 20) | 45% / 63% / 79% at thresholds 300 / 500 / 1000, independent of the reading |
| Region-wise sampling noise called selection | at most 26%, independent of the reading |

**The low power is due to the structure test, not the anomaly test:**
- At φ = 20 and 300 cells, 34 of 48 datasets are "undecidable".
- Of these, 32 are because there are fewer than 10 pairs of adjacent anomalous karyotypes, so the Moran test cannot be done.
- These datasets have a median of 8 anomalous karyotypes, which shows that the anomaly test under the envelope has already seen the anomaly; what is missing is the step that separates noise from selection.

**Meaning for real data:**
- Under this reading, all four could be satisfied only at threshold 1000.
- With ≥ 1000 cells at both time points: 0 of 277 Salehi trajectories, 3 of 41 NUMBAT datasets.
- Threshold 500: Salehi 33 trajectories, NUMBAT 8 datasets. Threshold 300: Salehi 123 trajectories, NUMBAT 18 datasets.
- So even with the corrected reading, almost all existing real data would fall into "undecidable, take the normal path".
