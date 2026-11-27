# Why the forecast is too concentrated: plan and decision rules set in advance

The user replied "3" on 2026-11-20 (after 08:38), i.e. item 3 of section 5 of `../decomp/report.md`: "find out why the forecast is too concentrated". This file was written on 2026-11-20 08:52 EST, before any new forecast. It is not changed afterwards; corrections go into the "Corrections" section at the end, with time and reason.

## Background

Conclusion of the E4 split:
- With the fit β kept at 5e-5 and the forecast missegregation rate raised to ×16–×128, the log score on later passages improves on every cell line (confirmation group +0.45 per cell).
- That is, starting from the 5e-5 landscape and forecasting with 5e-5 is too concentrated for the later passages: the karyotypes that actually appear get too low a forecast probability.

**How the forecast is made (`forecast_score.R`, same as the paper benchmark):**
- alfakR's ABM: 1e5 initial cells, allocated by the raw frequencies of the parent passage sample;
- when the cell number reaches 1e7, 1% is kept at random;
- in each division, each chromosome copy missegregates with probability p, so one division can change several chromosomes at once;
- a daughter cell whose karyotype is not in the landscape is discarded;
- average over 5 replicates.
- Scoring is done only on the landscape (i.e. the common support); cells outside the landscape are not scored.

## Four candidate causes

| ID | Cause | Mechanism |
|---|---|---|
| H1 | Incomplete sampling of the parent | Only 49–400 cells were measured in the parent passage; rare karyotypes present in the population but not measured start from 0 in the forecast |
| H2 | Passage bottleneck (drift) | The ABM bottleneck is 1e5 cells, almost no drift; real PDX passages may be much smaller |
| H3 | Multi-chromosome changes | New karyotypes in later passages are far from the parent and cannot be reached in 15–60 days at a missegregation rate of 5e-5 |
| H4 | Karyotypes outside the landscape | Later-passage cells fall outside the landscape, cannot be forecast and are not scored; in the ABM, daughter cells that leave the landscape are discarded |

## Method

**Subjects:**
- The 155 Salehi trajectories with later passages (137 actually usable), only later passages ≤ 60 days.
- The fit is fixed to the ×1 (5e-5) alfak results (minobs 10), same as in the split study.

**Script:** a new `disp_forecast.R`; the forecast and scoring code is copied from `forecast_score.R`, with only three additions:
- the initial frequencies can be replaced (I1);
- the ABM bottleneck can be replaced (I2);
- it additionally saves, on each later passage, the forecast probability, karyotype class and distance of the observed karyotypes.

**Check:** without changes (base in the table below), `ll_fc_s_e2` at F = ×1 and ×16 must be identical to the results of the split study (same seeds). If different, stop.

**9 forecasts per trajectory:**

| Variant | Forecast missegregation rate F | Content |
|---|---|---|
| base | ×1, ×16 | Same as before |
| I1 smooth | ×1, ×4, ×16 | Parent frequencies plus the unmeasured part (see below) |
| I2 drift1e3 | ×1, ×16 | ABM starts with 1e3 cells, keeps 10% at 1e4 (bottleneck about 1e3) |
| I2 drift1e4 | ×1, ×16 | ABM starts with 1e4 cells, keeps 10% at 1e5 (bottleneck about 1e4) |

**Initial frequencies of I1:** x0 = (1 − u) × parent frequency + u × g.
- **u:** the Good–Turing estimate of the unseen part of the parent sample, i.e. the number of karyotypes seen only once ÷ the number of parent cells, capped at 0.5.
- **g:** distributed over the karyotypes in the landscape that the parent sample does not have.
  - Allocated by their summed counts in the earlier passages of the same lineage.
  - When there are no such karyotypes, shared equally among the one-step neighbours, in the landscape, of the parent sample karyotypes.

## Diagnostics (description only, no decision)

**D1 Which class of cells the gain comes from.** Cells of a later passage that are in the landscape fall into four classes:

| Class | Definition |
|---|---|
| A | Karyotypes present in the parent sample |
| B | Not in the parent sample, but seen in an earlier passage of the same lineage |
| C | Never seen before, one step from some parent sample karyotype (copy-number L1 distance = 1) |
| D | Never seen before, ≥ 2 steps from all parent sample karyotypes |

Reported for each class: fraction of cells, log score per cell at F = ×1 and ×16, and the contribution to the total gain (fraction × difference).

**D2 Parent sampling (H1):**
- Compare two things, grouped by the day of the later passage:
  - the fraction of cells in the later passage with karyotypes absent from the parent sample (including cells outside the landscape);
  - the Good–Turing unseen part u of the parent sample.
- If the two are close at 15 days, the new karyotypes are mainly ones not measured in the parent; only the part of the fraction that grows with time is newly generated or expanded.

**D3 Distance spectrum (H3, H4):**
- The distance of later-passage cells from the nearest parent sample karyotype, in four bins 0, 1, 2, ≥ 3, given as fractions, including cells outside the landscape.
- Compared with the forecast probability in each bin (F = ×1, ×16, inside the landscape).
- The fraction of cells outside the landscape is also reported.

**D4 Drift (H2):**
- Take the karyotypes of class A with forecast expectation ≥ 5 cells, compute the Pearson dispersion φ = Σ(y − nq)²/(nq) ÷ (number of terms − 1).
- Reported by day; with sampling noise only, φ ≈ 1.
- A rough estimate of the effective bottleneck N_e per passage from φ − 1 ≈ n × number of passages / N_e. Error in the selection estimate also raises φ, so this is only a lower bound.

## Decision rules for the interventions (set in advance)

**Gain:** the difference in `ll_fc_s_e2` relative to base F = ×1, averaged over the later passages; 95% interval by cluster bootstrap over trajectories, 2000 times. Groups as before: development group SA609, SA039, confirmation group SA906, SA532, SA535, SA1035.

**An intervention "explains a substantial part" if both hold:**
- at F = ×1, the confirmation group gain ≥ 0.5 × 0.45 (the leave-one-out gain of changing only F in the split study), and the interval excludes 0;
- development group gain ≥ 0.

The two bottlenecks of I2 are decided separately.

**Report only:**
- whether the gain still rises when F = ×4, ×16 is added to the intervention, i.e. whether the optimum of F moves towards ×1;
- the gain decomposition of base F = ×16 (D1).

## Rules

- Code in `dev/realdata/beta/disp/`, not committed; the package is not changed, cd342db is used.
- HPC: xxlarge, one forecast per task; the add-hours monitor in small QOS, 3 days; extra arguments separated by colons.
- After the runs finish and `report.md` is written, stop and wait for the user to decide.

## Stop rules

- The base check does not agree;
- a build or script fails for an unknown reason;
- > 5% of the tasks of an array fail, and it is neither a resource problem nor an error of the package itself;
- a situation the rules do not cover.

## Corrections
