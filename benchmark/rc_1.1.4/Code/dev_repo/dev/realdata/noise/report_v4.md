# Pre-fit preprocessing v4: sensitivity-first noise screening, and handling when there is noise

**Basis:**
- Plan: `PLAN_v4.md`, fixed 2026-11-20 16:05; the 16:24 correction added the robust-φ̂ candidate T2r, before any T2 or T2r fit.
- The numbers come from `v4/v4_eval.R`; results are in `results/realdata/noise/v4/evaluation.rds` (confirmation half) and `evaluation_dev.rds` (development half).
- When the evaluation script first ran, the pattern for reading the Salehi prediction files accepted only lowercase `.rds` and missed the `.Rds` files. It was corrected and rerun. This was a reading error; neither the criteria nor the data were touched.

## Conclusion: none of the three treatments passed; by the plan, the screening serves only as a data-quality report, and KARST keeps its defaults

The verdict looks at the simulations of the confirmation half and at the Salehi confirmation group:

| Criterion | T1 (missegregation rate set from the data) | T2 (default missegregation rate + noise term, original φ̂) | T2r (default missegregation rate + noise term, robust φ̂) |
|---|---|---|---|
| Safety 1: selection simulations wrongly flagged as noise, median error ratio ≤ 1.02 | 1.00, passed | 1.19, **failed** | 1.10, **failed** |
| Safety 2: e3 wrongly flagged as noise, increase in false-selection fraction ≤ 0.02 | +0.0002, passed | −0.0002, passed | −0.0002, passed |
| Useful: data with injection, median error ratio ≤ 0.95 and upper interval limit < 1 | 1.000 [1.000, 1.045], **failed** | 0.774 [0.730, 0.821], passed | 0.796 [0.757, 0.853], passed |
| Real data: Salehi confirmation group, gain ≥ 0 and lower interval limit > −0.1 | −0.146 [−0.270, −0.041], **failed** | −0.028 [−0.205, +0.103], **failed** | −0.015 [−0.182, +0.102], **failed** |

**Sample sizes:**
- Safety 1: 238 datasets; Safety 2: 76; Useful: 1122;
- Salehi confirmation group: 42 trajectories, 72 later passages.

**Development half:** same conclusion.
- T1: 1.00, +0.0003, 1.000;
- T2: 1.20, −0.0001, 0.749;
- T2r: 1.14, −0.0001, 0.791.

**Fit failures:** 16 datasets, failing in all four fits, all with "at some time point no cell is in the model domain". This is a limitation of the data themselves, unrelated to the treatment, about 0.5%.

## 1. Option 1 (missegregation rate set from the data): no benefit, clearly worse on real data

- In the noise class, 49% of β̂ are exactly the default 5e-5, and then T1 is identical to T0.
- The other 51% of β̂ are higher, and for these datasets the median error ratio is **1.46**, clearly worse.
- The per-cell log score of the Salehi confirmation group drops by 0.146, with an interval excluding 0. The drop comes mainly from SA1035 (−0.42) and SA535 (−0.28).
- This agrees with earlier results: raising only the missegregation rate used for fitting does not improve prediction (`beta/decomp/report.md`, confirmation group +0.03, interval contains 0).
- The reason also matches the expectation in PLAN: data judged as noise are, by definition, not explained even at β = 1e-2. Raising β does not remove these anomalies, and instead explains real fitness differences as missegregation.

## 2. Option 2 (default missegregation rate + negative binomial noise term): clearly better when there really is noise, but two other kinds of data get worse

**Where it improves:** median error ratio, confirmation half, data flagged as noise:

| Data | T2 | T2r |
|---|---|---|
| Neutral e3 (with passage drift), no injection | 0.52–0.72 | 0.62–0.71 |
| e1 (selection + drift), no injection | 0.82–1.12 | 0.68–1.00 |
| Independent noise of φ = 20 added, each set | 0.54–0.92 | 0.64–0.93 |

**Where it gets worse:** C2 and e4, medians without injection:

| Group | T2 | T2r |
|---|---|---|
| C2 | 1.55 | 1.44 |
| e4 | 1.16 | 1.12 |
| e2 | 1.07 | 1.02 |

- With region-wise sampling noise added, C2 also gets worse: T2r is 1.14–1.48; e4 at 1000 cells is 1.11.
- Safety 1 failing is caused mainly by C2 and e4.

**What C2 and e4 have in common: the true missegregation rate is higher than the default.** C2 is 3–30 times higher, e4 4 or 16 times higher.

My conjecture (not verified separately):
- in these two groups, part of the "anomaly" comes from the higher missegregation rate, that is, more new karyotypes are produced;
- when the fit still uses the default missegregation rate, the negative binomial term flattens these new karyotypes as noise, and with them weakens the real fitness differences.

**Salehi: neither noise term gives a visible gain, and the intervals are wide.**

| Cell line | T2 | T2r |
|---|---|---|
| SA906 | +0.10 | +0.06 |
| SA532 | +0.08 | +0.08 |
| SA609 | −0.05 | −0.02 |
| SA1035 | −0.10 | −0.04 |
| SA535 | −0.28 | −0.19 |

**The robust φ̂ is better than the original definition, but not enough to pass.** It lowers Safety 1 from 1.19 to 1.10, the Salehi interval is also slightly narrower, and "Useful" is slightly weaker (0.796 vs 0.774).

## 3. The screening itself: on Salehi, the calls between parallel lineages are close to random

- Trajectory pairs of the same cell line, with the same treatment at both ends and non-overlapping time points: 1232 pairs in total. The agreement rate of the call "noise or not" is only **0.511**.
- If two trajectories were each called noise with probability 43%, independently, the expected agreement rate would be about 0.51.
- That is, on the current Salehi data, this call is almost not reproducible between parallel lineages.

## Report-only part

### Fractions of each class

Confirmation and development halves together; see the first table output by `v4_eval.R`.

| Data | Normal | Needs higher β | Noise | Selection-like |
|---|---|---|---|---|
| Salehi (155 trajectories) | 68 | 20 | 67 (43%) | 0 |
| NUMBAT (41) | 17 | 1 | 9 | 14 |
| Simulation (4896) | — | — | 2739 | — |

### φ̂ against injection

Data flagged as noise, quantiles 25% / 50% / 75%:

| Injection | Cells per time point | Original φ̂ | Robust φ̂_r |
|---|---|---|---|
| None | 300 | 13 / 56 / 237 | 7.8 / 16.6 / 38.1 |
| None | 1000 | 10 / 26 / 119 | 6.7 / 10.9 / 18.4 |
| φ = 20 | 300 | 62 / 92 / 168 | 18.4 / 31.4 / 67.0 |
| φ = 20 | 1000 | 69 / 103 / 165 | 18.3 / 35.4 / 57.6 |
| Region σ = 1 | 300 | 16 / 39 / 175 | 7.6 / 11.8 / 26.2 |
| Region σ = 1 | 1000 | 14 / 41 / 352 | 7.1 / 11.1 / 18.4 |

- The scale of the robust φ̂_r is more reasonable: in data with φ = 20 injected, both time points have noise added, and there is also drift of their own, so φ̂_r is about 31–35.
- But without injection, φ̂_r is also 11–17. These data have drift of their own, or selection.

### The class that needs a higher β, with T1

Median error ratio:
- e1: 1.59 (32);
- e3: 1.45 (19);
- C2 β 1.5e-4: 1.46 (6);
- e2: 1.02 (9);
- e4 ×4, C2 β 5e-4: 1.00;
- e4 ×16: 1.34 (1);
- only the 1 dataset of C2 β 1.5e-3 is 0.61.

So even for this class, raising the missegregation rate from the data mostly makes things worse.

## Next steps (waiting for your decision; I have not done them)

The plan stops here. Possible directions:

1. **Accept the current state:** preprocessing serves only as a data-quality report, and KARST keeps the default missegregation rate without a noise term.
2. **Add the noise term only where it really helps:**
   - Option 2 improves clearly on neutral and selection data with drift (error down by about 20–50%), and gets worse only on C2 and e4, where the true missegregation rate is high.
   - If a signal can be found that, with two time points, distinguishes "missegregation rate too high" from "drift / sampling noise", the noise term could be added only for the latter. For example, whether new karyotypes are concentrated on the one-step neighbours of existing karyotypes.
   - This needs a new plan; the verdicts cannot be changed after the fact with these data.
3. **Improve the estimate of the robust φ̂:**
   - At present φ̂_r is 11–17 even without injection, partly from selection.
   - One could consider estimating φ only from karyotypes far from existing clones, which are unlikely to be under selection.
   - This also needs a new plan.

## Files

- Code:
  - `dev/realdata/noise/R/envelope_v4.R`;
  - `dev/realdata/noise/v4/`: `screen_task.R`, `screen_real.R`, `fit_task.R`, `make_fit_tasks.R`, `check_identity.R`, `salehi_fit.R`, `salehi_horizon.R`, `v4_eval.R`;
  - patch `dev/patches/V4_noise_phi.patch`, stacked on `R1_nb1_drift.patch`, base 6ae2a24.
- Library: HPC `karst_dev/lib/karst_v4`. With `noise_phi` not given or given as 1, it is bit-identical to 1.1.3; checked.
- Results: HPC `results/realdata/noise/v4/{screen,fits,salehi}`, local copy in `dev/results/realdata/noise/v4/`.
- Jobs:

  | Content | Job ID |
  |---|---|
  | Screening | 22188519, 22188522 |
  | Build | 22189041 |
  | Bit-identity check | 22202681, 22202683 |
  | Trial run | 22205524 |
  | Development half | 22224242 (5452) |
  | Confirmation half | 22225935 (5930) |
  | Salehi fits | 22206030 |
  | Salehi prediction | 22206035 |

  All COMPLETED.
