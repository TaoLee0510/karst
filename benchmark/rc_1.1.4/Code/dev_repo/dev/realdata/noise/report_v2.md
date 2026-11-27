# Pre-fit preprocessing v2 (anomalous fluctuations and their structure): P4 results

**Basis:**
- The plan is in `PLAN_v2.md`, fixed on 2026-11-20 at 12:30, with no corrections.
- The numbers come from `v2_eval.R`, with results in `results/realdata/noise/v2_evaluation.rds`.
- 3648 analyses in total, all successful.

## Conclusion: P4 did not pass; stopped as the plan says, not applied to Salehi

| Criterion | Result | Passed |
|---|---|---|
| 1. e3 (neutral, no injection) judged sampling-dominated ≤ 10% | 20.5% | No |
| 2. e1, e2, e4, C2 with selection judged sampling-dominated ≤ 10% | 5.9% | Yes |
| 3. e3 + independent overdispersion φ = 20 judged sampling-dominated ≥ 70% | 34.7% | No |

**But on a closer look, the reasons for the failure are more specific than "the method does not work":**

1. **I set criterion 1 wrongly: e3 is not free of noise.**
   - e3 is neutral, but it has passage bottlenecks (200, 2000, 20000 cells), i.e. it has drift, and drift is exactly one of the kinds of noise the user defined.
   - With 1000 cells per time point, the fraction judged sampling-dominated varies with the strength of drift:

     | Bottleneck | 200 | 2000 | 20000 |
     |---|---|---|---|
     | Judged sampling-dominated | 58% | 50% | 12% |

   - The method in fact detected real drift; at bottleneck 20000 (almost no drift) it is 12%, close to the 10% threshold.
   - The verdict still counts as "not passed" under the criteria set in advance; this breakdown only explains the reason for the failure.
2. **With only 200 cells per time point, there is almost no power.**
   - In e3 the median number of testable karyotypes is only 6, with 4 neighbouring karyotype pairs and 1 anomalous karyotype.
   - The decisions almost all fall into "no clear anomaly" or "cannot tell", which go through the default pipeline.
   - Criterion 3 failed mainly because of this: at φ = 20, 0% are judged sampling-dominated with 200 cells, and 52–85% with 1000 cells.
   - This in itself fits your decision rule (if it cannot tell, use the default), but it shows that this preprocessing is useful only when there are many cells per time point (about 1000).
3. **Sampling by region is indeed often misjudged as selection; this is the limitation that was to be measured.**
   - For e3 with σ = 1 regional noise and 1000 cells, the fraction judged "selection":

     | Bottleneck | 200 | 2000 | 20000 |
     |---|---|---|---|
     | Judged selection | 15% | 33% | 56% |

   - The weaker the drift, the more easily regional noise is taken for selection. That is, when a clone family in one region is over- or under-sampled, it looks "patchy" just like real selection, and the structure test cannot separate them.
4. **Real selection is rarely recognised as "patchy" in these simulations.**
   - The simulations have low missegregation rates, and selected new clones rarely have patches of close relatives.
   - e1 with 1000 cells is judged "selection" in only 15–17%; e2 (clone replacement) in 54%.
   - C2 at β_sim = 1.5e-3 is 100% judged "no clear anomaly": the true missegregation rate itself is inside the envelope, and the envelope is very wide, so the anomalies cannot be detected.
   - So the criterion "selection produces patches of new karyotypes" is essentially unusable in simulations with low missegregation rates. At high missegregation rates, on the other hand, the envelope is so wide that even selection cannot be seen.

## 1. Decision fractions per set

All cell counts pooled. "No anomaly" means anomalous karyotypes < 3 or excess fraction < 0.05.

| Set | Injection | n | No anomaly | Cannot tell | Selection | Sampling |
|---|---|---|---|---|---|---|
| e1 | none | 288 | 0.54 | 0.27 | 0.08 | 0.11 |
| e2 | none | 192 | 0.25 | 0.47 | 0.28 | 0.01 |
| e3 | none | 288 | 0.55 | 0.21 | 0.04 | 0.21 |
| e4 | none | 288 | 0.62 | 0.20 | 0.11 | 0.07 |
| C2 β 1.5e-4 (bottleneck 300 / 3000) | none | 48 / 48 | 0.54 / 0.79 | 0.35 / 0.15 | 0 / 0.04 | 0.10 / 0.02 |
| C2 β 5e-4 (bottleneck 300 / 3000) | none | 48 / 48 | 0.60 / 0.98 | 0.33 / 0.02 | 0.02 / 0 | 0.04 / 0 |
| C2 β 1.5e-3 (both bottlenecks) | none | 96 | 1.00 | 0 | 0 | 0 |
| e3 | overdispersion φ = 5 | 288 | 0.24 | 0.42 | 0.02 | 0.32 |
| e3 | overdispersion φ = 20 | 288 | 0.14 | 0.50 | 0.01 | 0.35 |
| e3 | region σ = 0.5 | 288 | 0.40 | 0.32 | 0.05 | 0.22 |
| e3 | region σ = 1.0 | 288 | 0.32 | 0.38 | 0.18 | 0.13 |
| e1 | overdispersion φ = 20 | 288 | 0.32 | 0.45 | 0.07 | 0.17 |
| e1 | region σ = 1.0 | 288 | 0.42 | 0.34 | 0.14 | 0.09 |

After adding noise to C2, most are judged "cannot tell": there are too few testable neighbouring karyotype pairs; at β 1.5e-3 the median is only 0–1 pairs.

## 2. A property of the envelope

- **At the β = 1×10⁻² end, the lower edge is almost 0.** With many divisions, almost no karyotype can stay unchanged, so any fall is inside the envelope.
  - Example: with an interval of 18 days and 2ρΔt = 36, q_min is about 1e-7.
- **Therefore anomalies can essentially be seen only among "rises";** only when the interval is very short can a "fall" be judged anomalous.
- This is a direct consequence of "letting missegregation explain the change as far as possible".

## 3. Summary with respect to the user's idea

- **"Whether the anomalies are patchy" can distinguish drift (independent fluctuations) from patchy change, but only if there are about 1000 cells per time point.**
  - In neutral data, the fraction judged sampling-dominated rises with the true strength of drift (58% → 12%).
  - After adding independent overdispersion, 52–85% are judged sampling-dominated with 1000 cells.
- **"Patchy" cannot distinguish real selection from sampling by region.** When drift is weak, 56% of the regional noise is judged selection.
- **At low missegregation rates, real selection is not necessarily patchy.** It often shows up as the expansion of a single clone.
- **At high missegregation rates, the envelope is too wide, and neither selection nor noise looks anomalous.**

## 4. What could be done next (for you to decide)

1. **Per the decision rule, with two time points the default is the normal pipeline** (β fixed at the default). P4 did not pass, so the preprocessing is not used for automatic routing for now.
2. **To continue in this direction, two corrections could be made (a new plan is needed):**
   - change criterion 1 to "neutral data without drift (bottleneck 20000) judged sampling-dominated ≤ 10%", and add another one: the stronger the drift, the more are judged sampling-dominated;
   - enable the preprocessing only on data with ≥ 500 cells per time point.
   - These are changes made after seeing the results; they must be validated again on new simulated data (new landscapes), and passing cannot be claimed on the same data.
3. **That sampling by region cannot be separated from selection is a fundamental limitation of this approach.** To distinguish them would require multiple samples from the same time point (multi-region sampling), or three or more time points: selection keeps going in one direction, while regional sampling errors do not.
4. **P5 (externally estimated φ + negative binomial) is not done for now**, because it depends on P4 passing first.
