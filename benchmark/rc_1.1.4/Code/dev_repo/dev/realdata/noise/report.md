# Pre-fit noise assessment: results of N1–N3

**Basis:**
- The plan is in `PLAN.md` (fixed 2026-11-20 10:52); confirmations and refined definitions were added at 11:40, all before any r̂ was computed.
- β_hi = 1×10⁻², confirmed by the user; the rationale is in `n0_beta_hi.md`.
- Numbers come from `n23_eval.R`; results are in `results/realdata/noise/evaluation.rds`.

**Notation:**
- r̂ and the truth are written as multiples of 5e-5.
- **Definition of r̂:** r̂ = −ln[(1 − n_new)/(1 − u)] ÷ (C̄ · D).
  - n_new: fraction of cells at the later time point whose karyotype is not in the earlier time point's sample;
  - u: the Good–Turing unseen fraction of the earlier time point's sample;
  - C̄: mean copy number;
  - D: number of divisions.

## Conclusion: neither N2 nor N3 passed; N4 is not entered

| test | criterion | result | passed |
|---|---|---|---|
| N2 accurate | r̂ within a factor of 2 of the truth ≥ 70% | 37.7% | no |
| N2 monotone | median rises with the truth | e1 ×1 < e4 ×4 < e4 ×16; C2 ×3 < ×10 < ×30 | yes |
| N2 no false flags | flagged as noise-dominated in simulations ≤ 10% | 13.4% | no |
| N2 drift | N̂_e of e3 within a factor of 3 of the true bottleneck ≥ 70% | 37.8% | no |
| N3 prediction | confirmation-group gain ≥ 0.227, interval excludes 0 | +0.227 [0.114, 0.323] | yes (just over the line) |
| N3 development group | gain ≥ 0 | +0.94 | yes |
| N3 parallel lineages | r̂ differing within a factor of 2 ≥ 70% | 14.6% (3703 pairs) | no |

**In one sentence:** the "apparent diversification speed" measured directly from two time points measures the missegregation rate accurately in the neutral case. But under selection it mainly reflects selection, not missegregation or noise, so it cannot be used to separate noise from missegregation.

**What was seen in detail:**
1. **Accurate when neutral, strongly biased under selection.**
   - Neutral set e3 (truth ×1): 76% within a factor of 2, median 0.92 times.
   - Under selection, the bias goes in two opposite directions:
     - e1 (also truth ×1, with selection): median only ×0.04, r̂ = 0 in 42% of datasets. Negative selection cleared out the new karyotypes, so it looks as if there were no diversification.
     - e2 (clonal replacement): 66% saturated (the later time point is almost all new karyotypes), r̂ = ∞, flagged as "noise-dominated". In fact it is a selective sweep.
   - So r̂ mixes missegregation and selection, and the two cannot be separated with two time points.
2. **On high-missegregation simulations, r̂ roughly follows the truth, but with much scatter.**
   - C2: ×3 → median 1.25, ×10 → 8.4, ×30 → 48; 42–63% within a factor of 2.
   - This is comparable to A1 of E1 in the β study (38–60%), not better.
3. **The simulations flagged as noise-dominated all come from saturation**: e2, C2 ×30 (19%), e4 ×16 (21%). The cause is not noise but new karyotypes expanding to become the bulk.
4. **The drift lower bound is accurate only for small bottlenecks.**
   - In e3 with true bottleneck 200, the median of N̂_e / truth is 0.79, with 63% within a factor of 3;
   - with bottlenecks 2000 and 20000, the median is only 0.28 and 0.08, with only 39% and 13% within a factor of 3.
   - With large bottlenecks drift is very weak; F_c comes mainly from new karyotypes produced by missegregation and from sampling details, and N̂_e is greatly underestimated.
5. **On Salehi, r̂ as the forecast missegregation rate is indeed useful, but unstable.**
   - Using the r̂ computed from each trajectory's own two input time points as the forecast missegregation rate, without any later passage: development group +0.94, confirmation group +0.23. The leave-one-out gains when only the forecast missegregation rate is changed are +1.03 and +0.45.
   - At a division rate of 0.141/day, the confirmation group gives +0.29.
   - But between parallel lineages of the same cell line and the same treatment, only 15% of r̂ differ within a factor of 2. On single trajectories it is very unstable.
6. **With β_hi = 1×10⁻², no Salehi or NUMBAT data are flagged as noise-dominated.**
   - Median r̂ per Salehi cell line: SA039 ×0.07, SA532 ×1.3, SA906 ×1.6, SA1035 ×2.0, SA609 ×14, SA535 ×15.
   - All far below ×200.

**So, regarding the user's idea:**
- **The step "measure dispersion with a missegregation rate close to the physiological upper limit" can be done.** The apparent diversification of Salehi is all within the physiological range (at a division rate of 1/day).
- **The step "classify anything above the upper limit as noise" finds no use in the data:**
  - no real data exceed the upper limit;
  - the simulations that exceed it are in fact all selective sweeps, not noise.
- **"An additive noise rate ν" is still a reasonable form, but it cannot be determined from the data this way.** With two time points, ν is mixed with selection.

## 1. N2: r̂ on simulations (two time points)

| set (truth) | datasets | within 2× | median r̂ | quartiles | r̂ = 0 | saturated | flagged as noise | median u | median n_new |
|---|---|---|---|---|---|---|---|---|---|
| e1 (×1) | 288 | 0.17 | 0.04 | 0 – 0.33 | 0.42 | 0.01 | 0.01 | 0.08 | 0.09 |
| e2 (×1) | 192 | 0.06 | ∞ | 3.2 – ∞ | 0.03 | 0.66 | 0.66 | 0.09 | 1.00 |
| e3 (×1) | 288 | 0.76 | 0.92 | 0.67 – 1.31 | 0.06 | 0.00 | 0.00 | 0.07 | 0.30 |
| e4 (×4) | 144 | 0.22 | 0.65 | 0.01 – 2.9 | 0.24 | 0.02 | 0.02 | 0.36 | 0.38 |
| e4 (×16) | 144 | 0.26 | 14.8 | 2.6 – 43 | 0.04 | 0.21 | 0.21 | 0.50 | 0.93 |
| C2 (×3) | 96 | 0.42 | 1.25 | 0.68 – 1.9 | 0.04 | 0 | 0 | 0.19 | 0.28 |
| C2 (×10) | 96 | 0.63 | 8.4 | 4.8 – 12.5 | 0.01 | 0 | 0 | 0.50 | 0.72 |
| C2 (×30) | 96 | 0.61 | 48 | 32 – 92 | 0 | 0.19 | 0.19 | 0.50 | 0.99 |

- **Difference between 200 and 1000 cells per time point:** e3 is accurate in both (medians 1.16 and 0.82); e1 is low in both (0 and 0.07). The amount of sampling is not the main problem.
- **Segments between adjacent time points (report only):** e3 median 0.92, within a factor of 2; e1 0.31; e4 6.1. With short intervals, the underestimate of e1 is reduced.

## 2. N3: Salehi

**Using r̂ as the forecast missegregation rate** (fit at ×1, taking the nearest level on the forecast grid already in the decomposition study), per-cell log score gain relative to ×1:

| division rate | group | later passages | trajectories | gain [95% interval] |
|---|---|---|---|---|
| 1/day | development | 136 | 51 | +0.94 [0.79, 1.11] |
| 1/day | confirmation | 165 | 85 | +0.23 [0.11, 0.32] |
| 0.141/day | development | 136 | 51 | +1.01 [0.84, 1.18] |
| 0.141/day | confirmation | 165 | 85 | +0.29 [0.16, 0.43] |

- **The threshold is passed just barely:** the threshold is half of 0.4531, i.e. 0.2266; the confirmation group is 0.2270, only 0.0004 higher.
- **Distribution of the choices (1/day):**
  - SA609 mostly at ×16–×256; SA535 at ×8–×32;
  - SA906 and SA1035 at ×1/4–×4; SA532 scattered.
- **Medians per cell line (1/day):**

  | cell line | r̂ | u | n_new | one step | two steps | ≥ 3 steps | N̂_e (per passage) |
  |---|---|---|---|---|---|---|---|
  | SA039 | ×0.07 | 0.02 | 0.03 | 0.02 | 0.00 | 0.00 | 31 |
  | SA1035 | ×2.0 | 0.10 | 0.20 | 0.14 | 0.03 | 0.01 | 31 |
  | SA532 | ×1.3 | 0.41 | 0.49 | 0.27 | 0.09 | 0.10 | 25 |
  | SA535 | ×15 | 0.50 | 0.85 | 0.21 | 0.19 | 0.43 | 18 |
  | SA609 | ×14 | 0.50 | 0.90 | 0.25 | 0.20 | 0.34 | 40 |
  | SA906 | ×1.6 | 0.11 | 0.40 | 0.15 | 0.04 | 0.04 | 161 |

- For SA535 and SA609 the unseen fraction reaches the upper limit 0.5; 85–90% of cells at the later time point are new karyotypes, and more than one third of these are more than 3 steps from the earlier time point. This agrees with the conclusion of the disp study.
- N̂_e is only 18–161 cells per passage. The e3 results show that this quantity is badly underestimated with large bottlenecks, so it can only serve as a lower bound, not as the true bottleneck.

## 3. NUMBAT (report only, 0.141/day)

| cancer type | median r̂ | u | n_new |
|---|---|---|---|
| GBM | ×0.18 | 0.01 | 0.06 |
| HNSC | ×0.68 | 0.11 | 0.18 |
| Ovarian | ×6.4 | 0.21 | 0.48 |
| Melanoma | ×59 | 0.21 | 0.50 |

None are flagged as noise-dominated.

## 4. What could be done next (for you to decide)

1. **N4 is not done.** By the plan, N4 is entered only if both N2 and N3 pass.
2. **r̂ as an empirical value for the forecast missegregation rate can be tested separately.** Without using later passages, it obtains half to nine tenths of the leave-one-out gain from changing only the forecast missegregation rate. But it is very unstable on single trajectories. A pooled estimate per cell line (all trajectories of the same cell line together) may be more stable; this needs new criteria fixed in advance.
3. **To really separate drift, sampling and selection, information beyond two time points is needed:**
   - three or more time points (the variance of drift grows linearly with time, while the effect of selection keeps pointing in one direction);
   - or parallel lineages (same starting point, each drifting on its own).
   - Salehi has parallel lineages, which can serve as the basis of the next design.
