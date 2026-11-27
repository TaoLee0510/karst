# Pre-fit preprocessing v2: anomalous fluctuations and their structure

The user confirmed on 2026-11-20 (after 12:10):
- the null hypothesis uses the envelope over β ∈ [5e-5, 1e-2];
- selection-dominated cases also go through the normal pipeline for now;
- the validation adds structured noise from "sampling by region", to measure how large this limitation is;
- the direction of the noise pipeline (P5) is "externally estimated φ + negative binomial fit".

This file was written on 2026-11-20 at 12:30 EST, before any computation. It is not changed afterwards; anything that needs correcting goes into the "Corrections" section at the end.

## The user's idea

With only two time points, do preprocessing first.
- If missegregation cannot be separated from other factors: go through the normal pipeline, with β fixed at the default.
- If they can be separated: go through the noise pipeline.

**What to look at:** at the second time point (especially when the interval is short), karyotype fractions that rise or fall far beyond the expectation under "neutral + missegregation only". The cause is either a strong fitness change or sampling error.

**Criteria for telling them apart:**
- If, after assuming an extremely high missegregation rate, many karyotypes still fluctuate anomalously, it is most likely a sampling problem.
- Real selection produces a batch of new karyotypes with similar fitness that are contiguous in karyotype space; the anomalies come in patches.
- Sampling errors are isolated.

## P1 Anomalies (per dataset, two time points t0, t1)

**Expectation under neutral + missegregation:**
- **Single division:** the per-chromosome daughter-cell transition probability is the same as in KARST (`chromosome_transition_probability`): of c copies, m ~ Binomial(c, β) missegregate, and each missegregated copy enters this daughter cell with probability 1/2.
- **Number of divisions in the interval:** under exponential growth, the number of divisions in the ancestral lineage of a random cell over time t is approximately Poisson(2ρt). So the transition over the interval Δt is

  P = Σ_m Pois(m; 2ρΔt) · K^m

  Here it is computed independently per chromosome type and the results are multiplied (ignoring the correlation from different chromosomes sharing division times); the copy-number cap is the observed maximum + 4.
- **Division rate ρ:** the simulations use 1/day (same as the fit); Salehi uses 1/day, with 0.141/day also reported; NUMBAT uses 0.141/day.
- **β grid:** 13 log-spaced values of β on [5e-5, 1e-2].
- **Expected fractions:** for each β, q_k(β) = Σ_j p0_j · P(j → k), where p0 is the frequency at time point t0; the variance of q due to sampling at t0 is (Σ_j p0_j P_jk² − q_k²) / n0.
- **Envelope:** q_max_k = max_β q_k(β), q_min_k = min_β q_k(β).
  - For karyotypes absent from time point t0, q_max gets an extra 3/n0, i.e. the 95% upper bound for a zero count, representing the possibility that they were present at t0 but not measured.

**Standardised residual z:**
- **Rise:** when x1_k ≥ 5 and x1_k > n1·q_max, z = (x1 − n1·q_max) / √(n1·q_max(1 − q_max) + n1²·Var(q) at the β where q_max is attained).
- **Fall:** when k appears at t0, n1·q_min ≥ 5 and x1_k < n1·q_min, compare with q_min in the same way.
- **Otherwise:** z = 0 (inside the envelope).

**Anomalies and excess:**
- Anomaly: |z| > 3.
- Excess fraction A = [Σ_rises (x1 − n1·q_max) + Σ_falls (n1·q_min − x1)] ÷ n1.

## P2 Structure (whether neighbouring karyotypes move in the same direction)

- **Nodes T:** testable karyotypes, i.e. x1 ≥ 5, or present at t0 with n1·q_min ≥ 5.
- **Edges:** pairs of karyotypes in T at copy-number L1 distance 1.
- **Statistic:** Moran's I = (|T| / |edges|) · Σ_edges z_a·z_b ÷ Σ_T z².
- **Null distribution:** z permuted at random within T 999 times, p = (1 + number of times the permuted I ≥ the observed I) / 1000. Fixed seed.
- **Report only:** among anomalous karyotypes, the fraction that have another anomalous karyotype in the same direction within 2 steps, and its permutation expectation.

## P3 Decision (thresholds set in advance)

| Outcome | Condition | Pipeline |
|---|---|---|
| No clear anomaly | anomalous karyotypes < 3, or A < 0.05 | normal pipeline, default β |
| Cannot tell | edges < 10; or neither of the two conditions below holds | normal pipeline, default β |
| Selection-dominated | p < 0.05 | normal pipeline (for now) |
| Sampling-dominated | p ≥ 0.20 and anomalous karyotypes ≥ 5 | noise pipeline (P5) |

## P4 Validation on simulations

**Base data:** the two-time-point inputs of the β study (first and last time points of the e sets, C2), 1344 datasets in total.

**Injected noise:** injected once per dataset, with the seed determined by the dataset name. The two time points are injected independently, then resampled multinomially at the original cell counts.
- **Independent overdispersion:** p' ~ Dirichlet(α·p_obs), α = (n − 1)/(φ − 1) − 1, φ = 5, 20.
- **Sampling by region:**
  - each observed karyotype is assigned to the nearest, by L1 distance, of the 5 karyotypes most abundant over t0 and t1 combined, forming "clone families";
  - each family is multiplied by a common factor exp(N(0, σ²)), σ = 0.5, 1.0, then renormalised.
- **What gets injected:**
  - e3 gets all four;
  - e1 and C2 get φ = 20 and σ = 1.0, report only.

**To pass, all must hold:**
1. e3 (neutral, no injection): ≤ 10% judged sampling-dominated.
2. e1, e2, e4, C2 with selection (no injection), pooled: ≤ 10% judged sampling-dominated.
3. e3 + independent overdispersion φ = 20: ≥ 70% judged sampling-dominated.

**Report only:**
- the results for φ = 5;
- under sampling by region, the fractions judged selection, sampling, and cannot tell (this is the limitation the user wants measured);
- the results for e1 and C2 with injected noise;
- results split by cells per time point (200 / 1000), by interval, and by set.

## Afterwards

- **If P4 passes:** apply to Salehi (277 trajectories; mainly 1/day, 0.141/day also reported) and NUMBAT, and report:
  - the distribution of decisions per cell line and per cancer type;
  - whether the decisions agree between parallel lineages;
  - which datasets are sampling-dominated.

  Then write a separate P5 plan: estimate φ from the residuals of the unstructured part, then fix φ and fit with the negative binomial.
- **If P4 fails:** write the report and stop.

## Rules

- Code in `dev/realdata/noise/`, not committed; the package is not changed, and the data are only read.
- HPC: xxlarge, add-hours monitor with small QOS, 3 days; extra arguments separated by colons.
- Stop after writing `report_v2.md`, and wait for the user.

## Limitations stated in advance

- **Sampling by region may be misjudged as selection.** A region is often a whole clone family, whose sampling error is also patchy in karyotype space. P4 will measure the fraction of such cases judged "selection".
- **The simulations have low missegregation rates, so the new karyotypes created by selection rarely have patches of close relatives.** In the simulations β is mostly 5e-5; for selectively expanded new clones, their closely related karyotypes rarely reach 5 cells, so the anomalies may look isolated and be judged "sampling". This is exactly what criterion 2 tests. C2 (β_sim 1.5e-4–1.5e-3) is closer to the high-CIN case the user has in mind and is reported separately.
- **Other effects:**
  - The division rate affects the envelope width: overestimating ρ makes the envelope wider, so more are judged "no clear anomaly", which is conservative.
  - Computing the transitions independently per chromosome is an approximation.

## Corrections
