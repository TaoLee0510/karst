# Pre-fit noise assessment: apparent diversification speed and an additive noise rate

The user confirmed on 2026-11-20 (after 10:41) that the two plans are done at the same time. Specific decision: there is no ready literature value for the physiological upper bound β_hi; N0 sets it from the literature, it is shown to the user first, and N1–N3 continue only after confirmation.

This file was written 2026-11-20 10:52 EST, before any analysis. It is not changed afterwards; anything that needs correcting goes into the "Corrections" section at the end, with the time and the reason.

## The user's idea (original meaning)

Add a preprocessing step before KARST processes the data: judge the degree of dispersion from the given data, combining the intervals between time points and the degree of change between karyotypes, to measure drift and selection.

- The missegregation rate has a physiological upper bound: the number of chromosomes is fixed, so cells cannot have a very high missegregation rate.
- Use a missegregation rate close to the physiological upper bound to assess the dispersion of the data.
- If even at such a high missegregation rate the differences in karyotype dispersion between the finite time points still cannot be explained, the data are very likely affected by large noise. The noise has two aspects: drift and sampling.
- In that case, add a parameter that is additive with the missegregation rate and acts jointly with it on the karyotype trajectories.

## Choice of approach

Do not estimate it by likelihood. The β studies have already shown that the two-time-point likelihood is almost flat in β and cannot pin down β.

Instead, measure the "apparent diversification speed" r̂ (per chromosome copy per division) directly from the data, and compare it with the physiological upper bound β_hi:
- r̂ ≤ β_hi: can be explained by missegregation within the physiological range;
- r̂ > β_hi: flagged as "noise-dominated".

The additive noise rate is ν̂ = r̂ − β, where β is the missegregation rate the model uses, default 5e-5.

## N0 set the constants (stop after this step, for the user to confirm)

- **β_hi:** look up measured missegregation rates of human cells in the literature, including:
  - near-diploid stable cell lines, CIN cancer cell lines, organoids;
  - tumours and PDX, inferred from single-cell karyotype diversity;
  - the upper bound at which proliferation can still be sustained.
- Units unified to "per chromosome copy per division"; state the original meaning, the conversion and the source of every number.
- **Model-internal upper bound (as a reference):** ln2 ÷ number of chromosome copies. Above it, a karyotype can no longer maintain itself by faithful division. About 1.2e-2 at 60 copies; the point where the alfak format conversion fails at ×256 is close to this.
- **Result:** written into `n0_beta_hi.md`; stop and wait for the user to confirm β_hi.

## N1 assessment quantities (input data only, no fit needed)

For each interval between adjacent time points (t0, t1):

- **Sampling:** Good–Turing unseen fraction of the t0 sample, u = number of karyotypes seen only once ÷ number of t0 cells, capped at 0.5.
- **Diversification:**
  - n_new = the fraction of t1 cells whose karyotype is not in the t0 sample;
  - reported separately by the distance to the nearest karyotype in the t0 sample: 1, 2, ≥ 3.
- **Number of divisions:** D = interval × division rate.
  - simulations use the true value (1 for the e group and C2);
  - Salehi uses KARST's setting in the paper benchmark (1/day), and the sensitivity at 0.141/day is also reported;
  - NUMBAT uses 0.141/day.
- **Apparent speed:** under the neutral missegregation model, the probability that a lineage's karyotype stays unchanged over D divisions is about exp(−C̄ r D), where C̄ is the mean total copy number of the t0 cells. After removing the sampling part:

  r̂ = −ln[(1 − n_new) / (1 − u)] / (C̄ · D)

  r̂ = 0 when n_new ≤ u.
- **Noise flag:** when r̂ > β_hi, recorded as noise-dominated, ν̂ = r̂ − β.
- **Lower bound on drift:** use the classical temporal method (the F_c of Waples 1989):
  - take the karyotypes present at both time points; after F_c removes the variance of the two samplings, the effective bottleneck per passage is N̂_e = number of passages ÷ F_c;
  - with two time points, drift and selection cannot be separated, so this is only a lower bound.

## N2 validation on simulations (the existing e1–e4 and C2, no new simulations)

It passes only if all of the following are met:
- **Accurate:** datasets with r̂ within a factor of 2 of the true missegregation rate ≥ 70%. Reported separately for each cell of the e group and each C2 setting; the verdict uses everything pooled.
- **Monotone:** the median of r̂ rises monotonically with the true missegregation rate (e1 ×1 < e4 ×4 < e4 ×16; C2 ×3 < ×10 < ×30).
- **No false flags:** the true values in the simulations are all below β_hi; datasets flagged as noise-dominated ≤ 10%.
- **Drift:** on the neutral set e3, N̂_e within a factor of 3 of the true bottleneck ≥ 70%.

Reported separately: the difference between 200 and 1000 parental cells.

## N3 test on Salehi (no new fits needed, and no new predictions either)

- Each trajectory uses the r̂ computed from its own two input time points as the prediction missegregation rate, and takes the nearest level by log distance on the split study's existing prediction grid (fit ×1, prediction ×1/4–×256).
- **To pass, all of the following must hold:**
  - the gain of the confirmation group relative to ×1 is ≥ 0.23 (half of the leave-one-out gain of 0.45 obtained when only the prediction missegregation rate is changed), and its interval excludes 0 (cluster bootstrap by trajectory, 2000 times);
  - the development-group gain is ≥ 0;
  - parallel lineages (as defined in B) whose r̂ differ by no more than a factor of 2 ≥ 70%.
- This step does not tune parameters on later passages. Only a pass shows that r̂ measures something real.
- For NUMBAT and GBM, only r̂, u and N̂_e are reported.

## N4 (only after both N2 and N3 pass; a separate plan is written for the user to confirm)

Put ν into the model, in two ways:
- **Control:** ν is fully added to β and used together with it everywhere. This is equivalent to raising β, and existing results suggest it does not work.
- **Candidate:** ν is added only on the side of diffusion to neighbours (including the prediction), and does not enter the outflow term, which would in turn change the fitness estimates. The basis is the split study: raising the diffusion of the prediction helps, raising the fit β does not.

## Rules

- Code in `dev/realdata/noise/`, not committed; the package is not changed. N1–N3 only read existing data and results.
- After N0, stop and wait for the user; after N3, stop, write the report and wait for the user.

## Limitations stated in advance

- **β and ν are the same thing in the transition term.** Separating them relies entirely on the assumption that "β does not exceed the physiological upper bound". So ν̂ means "apparent diversification beyond the upper bound", not a direct measurement of drift.
- **Selection makes r̂ deviate from the true missegregation rate:** negative selection makes new karyotypes disappear, so r̂ is too low; positive selection makes new clones expand, so r̂ is too high. N2 tests this on the simulations with selection.
- **D depends on the division rate:** the division rate of Salehi is unknown and r̂ is inversely proportional to it, so both 1/day and 0.141/day are reported.
- **The sampling correction has a large effect:** the median unseen fraction of the Salehi parental samples is 0.37, so the sampling correction has a large effect on r̂.

## Corrections
- **2026-11-20 11:40 EST:** the user confirmed β_hi = 1×10⁻² ("Confirm β_hi = 1×10⁻², start N1–N3"). The basis is in `n0_beta_hi.md`.
- **2026-11-20 11:40 EST (refined definitions, before computing any r̂):**
  1. **Input:** N2 uses the same two-time-point input as the β studies (the first and last time points for the e group, the two time points of C2); N3 uses the first and last passages of the Salehi trajectories. Results for each segment between adjacent time points are also reported.
  2. **Saturation:** when n_new = 1, r̂ = ∞; this counts as noise-dominated and also as not within a factor of 2.
  3. **Drift F_c:**
     - **Classes:** each karyotype with ≥ 5 cells in t0 is its own class, and the remaining karyotypes in t0 are merged into one class; in t1 only cells whose karyotype appeared in t0 are counted.
     - **Formula:** F_c = Σ(x − y)² / ((x + y)/2 − xy) ÷ (number of classes − 1) (Nei–Tajima).
     - **Effective bottleneck:** cells are sampled in the haploid way, so N̂_e (per passage) = P ÷ (F_c − 1/S0 − 1/S1); when the denominator is ≤ 0, N̂_e = ∞.
     - **P, S0, S1:** P is the number of passages between the two time points; S0 and S1 are the numbers of cells used in the calculation.
  4. **Number of passages and true bottleneck:**
     - e group: use the phase 0 passage replay (`phase0/process_noise/replay`) to get the number of passages in the interval; the true bottleneck is the harmonic mean of the numbers of surviving cells after each passage.
     - C2: one passage every 15 days; the bottleneck is the set value B.
     - Salehi: the number of passages is the difference between the labels of the two time points.
  5. **Agreement of parallel lineages:** both having r̂ = 0 counts as agreement; only one being 0 counts as disagreement.
  6. **N3 grid:** the N3 grid includes ×256. With fit ×1, all trajectories have a prediction at ×256.
