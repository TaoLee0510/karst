# Pre-fit preprocessing v5: frequent karyotypes only, a lower bound on the missegregation rate from neighbour growth

**Basis:**
- Plan `PLAN_v5.md`, fixed 2026-11-20 before any v5 computation. Confirmed by the user at about 19:55.
- Numbers come from `v5/step12_eval.R`; the output is in `dev/results/realdata/noise/v5/step12_eval.log`.
- Section 3 is a post-hoc description, from `posthoc.log` and `posthoc2.log`, not a decision.
- Computation: job 22236370, all 32 tasks COMPLETED. 5093 records in total (simulations 4896, Salehi 155, NUMBAT 42), all succeeded.

## Conclusion: step 2 failed on the development half; stopped according to the plan, steps 3 and 4 not done

**Decision (development half, simulations without injection, 816 datasets):**

| Version | Valid: L exceeds the true missegregation rate (all / worst set) | Useful: fraction triggered when the true missegregation rate is ≥ 4 times | Fraction triggered among all datasets |
|---|---|---|---|
| a | 0 / 0, pass | **0.3% (1 of 336), fail** (required ≥ 20%) | 0.12% |
| b | 0 / 0, pass | **0.3%, fail** | 0.12% |

- The lower bound is never wrong: no dataset has L above the true missegregation rate. But it almost never triggers, so it is of no use.
- Neither version passed; the confirmation half is no longer used for a decision.
- KARST keeps the default missegregation rate 5e-5.

## 1. Step 1: premise checks (report only)

### 1.1 Neutral e3: fraction of "significant changes" among frequent karyotypes

A significant change is one that no missegregation rate in [5e-5, 1e-2] can explain, with |z| > 3.

| Bottleneck (cells) | Cells per time point | Frequent karyotypes | Called significant change | Datasets with at least one |
|---|---|---|---|---|
| 200 | 300 | 399 | 62% | 94% |
| 200 | 1000 | 1035 | 59% | 100% |
| 2000 | 300 | 223 | 29% | 50% |
| 2000 | 1000 | 1018 | 52% | 90% |
| 20000 | 300 | 170 | **2.9%** | 8% |
| 20000 | 1000 | 676 | 16% | 54% |

- **Total 44.7%** (3521 frequent karyotypes, 288 datasets).
- By interval: 18 days 25%, 42 days 32%, 180 days 48%, 420 days 62%.
- **"Noise is a small perturbation" (≤ 5%) holds only at bottleneck 20000 with 300 cells per time point.** At bottlenecks 200 and 2000, 30–60% of the changes of neutral frequent karyotypes cannot be explained even by the upper-bound missegregation rate. These changes come from passage drift.
- One reservation: ρ = 1 is higher than the true division rate in the simulations (see 1.4). This makes the envelope's expected retention of existing karyotypes too low, and at long intervals may over-report some rises. I did not quantify this term separately.

### 1.2 Simulations with selection: whether the direction of significant changes agrees with the true fitness

| Set | Called significant | Agree |
|---|---|---|
| e2 | 1457 | 99.3% |
| e1 | 1271 | 72.3% |
| e4 ×4 / ×16 | 805 / 598 | 72.0% / 71.9% |
| C2 1.5e-4, bottleneck 300 / 3000 | 471 / 176 | 68.8% / 93.2% |
| C2 5e-4, bottleneck 300 / 3000 | 417 / 39 | 88.2% / 100% |
| C2 1.5e-3, bottleneck 300 / 3000 | 46 / 1 | 100% / 100% |

- **Total 81.8%** (5281), below the expected 90%.
- **By bottleneck (e sets): 61% at 200, 92% at 2000, 98% at 20000.**
- So "large change = selection" largely holds at bottleneck ≥ 2000. Under strong drift about 4 in 10 of the directions are wrong.
- **The envelope can only detect rises.** At the upper-bound missegregation rate, almost all existing karyotypes would leave, the expected retention is close to 0, and a fall can never be detected. Everything called significant is a rise; there is not a single fall in the simulations.

### 1.3 Salehi parallel lineages

- Of 1232 trajectory pairs, 729 pairs have frequent karyotypes in common.
- **According to the plan (envelope test):**
  - 533 shared karyotypes (302 pairs) are called significant in both, with 100% agreement in direction.
  - **This number carries no information.** The envelope only detects rises, so the direction necessarily agrees.
  - Of the shared karyotypes called significant, 62% are significant in only one of the trajectories.
- **Post hoc (two-proportion test |z| > 3, not a decision):**
  - 934 (510 pairs) are significant in both, with **76.7%** agreement in direction;
  - by cell line: SA609 0.96 (164), SA906 0.76 (598), SA532 0.62 (163);
  - significant in only one: 50%.
  - The agreement rate lies between 50% and the 75% written in the plan, towards the upper end, with large differences between cell lines.

### 1.4 Upper bound on the division rate

- True net growth rate of the reference karyotypes in the simulations: median 0.887, maximum 0.944, all below ρ = 1.
- So the inflow term in the lower bound is overestimated, and the lower bound is conservative. This is consistent with "never exceeds the true missegregation rate".

## 2. Step 2: why it almost never triggers

Statistics by set, both halves combined, report only. Ratio r = number of cells in the neighbour cloud ÷ number of cells of the reference karyotype; β_pt is the point value of the lower bound.

| Set | True multiple | No reference karyotype | β_pt > 0 | Triggered | Median β_pt (× default) |
|---|---|---|---|---|---|
| e1 | ×1 | 0.7% | 27% | 0 | −0.08 |
| e2 | ×1 | 0 | 53% | 0 | 0.00 |
| e3 | ×1 | 9% | 52% | 0 | 0.01 |
| C2 1.5e-4, bottleneck 300 / 3000 | ×3 | 0 / 0 | 25% / 39% | 0 / 0 | −0.30 / −0.24 |
| e4 ×4 | ×4 | 5% | 23% | 0.7% | −0.30 |
| C2 5e-4, bottleneck 300 / 3000 | ×10 | 6% / 52% | 34% / 31% | 2.1% / 0 | −0.84 / −0.90 |
| e4 ×16 | ×16 | 43% | 29% | 0.7% | −0.53 |
| C2 1.5e-3, bottleneck 300 / 3000 | ×30 | 65% / 99% | 44% / 0 | 0 / 0 | −0.57 / −28 |

**Two reasons:**

1. **The neighbour cloud usually does not grow relative to the reference karyotype; it stays flat or shrinks.**
   - The lower bound requires the neighbour cloud to grow clearly between the two time points. When the population is close to steady state the growth is 0; when the reference karyotype itself rises under selection, the neighbour cloud cannot keep up and the ratio actually falls.
   - Of the 5 datasets where version a triggered (4 simulated, 1 Salehi), in 3 the cell number of the reference karyotype was falling, just not reaching "significant fall" (e.g. 24 → 15), so the ratio rose. In the other 2 the neighbour cloud really grew: C2 1 → 18 cells, SA1035 3 → 36.
2. **The higher the true missegregation rate, the fewer frequent karyotypes are available.**
   - Fraction without a reference karyotype: C2 ×30 65–99%, e4 ×16 43%, C2 ×10 bottleneck 3000 52%.
   - The datasets that most need an upward adjustment are exactly the ones with no lower bound to compute.

**Versions with injected noise (report only):** dm20 and reg10 also never exceed the true missegregation rate, triggering 0–0.4%.

**Real data (report only):**
- Salehi: 1 of 155 triggered (SA1035, confirmation group, L ≈ 2.4 times). 16 / 67 SA609 trajectories have no reference karyotype.
- NUMBAT: 0 of 41 triggered, 5 have no reference karyotype.

## 3. Post-hoc description (not a decision, not verified)

**The level of the neighbour cloud rises monotonically with the true missegregation rate.** What is looked at here is the ratio itself, not its growth.

| Set | True multiple | Fraction of cells in frequent karyotypes (median) | r at the second time point (median) |
|---|---|---|---|
| e1 / e2 / e3 | ×1 | 0.82 / 0.76 / 0.74 | 0.12 / 0.19 / 0.23 |
| C2 1.5e-4, bottleneck 300 / 3000 | ×3 | 0.69 / 0.46 | 0.28 / 0.57 |
| e4 ×4 | ×4 | 0.44 | 0.34 |
| C2 5e-4, bottleneck 300 / 3000 | ×10 | 0.21 / 0.01 | 0.76 / 1.22 |
| e4 ×16 | ×16 | 0.02 | 0.53 |
| C2 1.5e-3, bottleneck 300 / 3000 | ×30 | 0.00 / 0.00 | 1.30 / 0.80 (only 34 and 1 datasets have a reference karyotype) |

**On Salehi, ordered by the same two quantities:**

| Cell line | Fraction of cells in frequent karyotypes | r at the second time point |
|---|---|---|
| SA906 | 0.82 | 0.16 |
| SA1035 | 0.74 | 0.32 |
| SA532 | 0.47 | 0.43 |
| SA535 | 0.26 | 0.57 |
| SA609 | 0.07 | 1.18 |

- This order is roughly the same as the direction selected with later passages in the β study: SA906 leans towards ×1/2–×1, the others towards ×8–×64 (`beta/report.md`).
- **But this cannot be used directly as an estimate of the missegregation rate.** At steady state r ≈ missegregation rate × number of routes ÷ fitness disadvantage of the neighbours, and the size of the disadvantage is unknown. The distribution of the disadvantage in the simulations is set by the same kind of GRF landscape, so it looks monotone on the simulations, but this need not hold for another kind of landscape. r also differs with the bottleneck (C2 ×3: 0.28 vs 0.57).
- Using it would need a new plan, and a test on different types of landscape of whether it transfers.

## Next steps (waiting for your decision; I have not done these)

1. **Accept the current state:** the missegregation rate fixed at 5e-5. Add "fraction of cells in frequent karyotypes" and "neighbour-cloud ratio r" to the data-quality report, as descriptions.
2. **New plan: correct the missegregation rate with the level (not the growth) of the neighbour cloud.**
   - The unknown fitness disadvantage has to be handled, e.g. by calibrating on several kinds of landscape (different wavelengths, non-GRF) and testing on landscapes not used in the calibration.
   - The rule of only adjusting upwards, capped at 10 times, can then be kept.
   - Risk: the calibration may only be valid for the type of landscape used in the simulations.
3. **The finding of step 1 itself:** "large change = selection" is largely reliable only when the passage bottleneck is ≥ 2000 cells. I do not know the passage bottleneck of Salehi; this needs your judgement.

## Files

- Code:
  - `dev/realdata/noise/R/anchor.R`;
  - `dev/realdata/noise/v5/step12_task.R`, `step12_eval.R`;
  - task table `v5/tasks/step12.tsv`.
  - The post-hoc description scripts `v5_posthoc.R` and `v5_posthoc2.R` are in the results directory and were not put in the repository. The by-set section of `posthoc.log` was computed wrongly because of repeated dataset names (different C2 settings, e1/e3 with the same names) and has been removed; for by-set numbers `posthoc2.log` is authoritative.
- Results:
  - HPC `results/realdata/noise/v5/anchor/`;
  - local `dev/results/realdata/noise/v5/`: `anchor/`, `step12_evaluation.rds`, `step12_eval.log`, `posthoc.log`, `posthoc2.log`.
- Job: 22236370 (32 tasks), library karst_v4 (used only to read the truth, no fitting).
