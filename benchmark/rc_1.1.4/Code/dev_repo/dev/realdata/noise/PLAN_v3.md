# Pre-fit preprocessing v3: corrected criteria, re-validated on new simulated data

The user confirmed on 2026-11-20 (after 13:03): "Confirmed, proceed entirely as you suggest". This includes:
- 12 new landscapes each for the e group and for C2;
- testing the two enabling thresholds 300 and 500 at the same time, and taking the smaller one that passes;
- not using the sealed held-out landscapes (seed 98000).

This file was written 2026-11-20 13:06 EST, before any new data were generated. It is not changed afterwards; anything that needs correcting goes into the "Corrections" section at the end.

## Unchanged parts

- The method is exactly as in `PLAN_v2.md` and `R/envelope.R`, with no code changes: envelope β ∈ [5e-5, 1e-2], z, Moran's I, decision thresholds (anomalous karyotypes ≥ 3, excess ≥ 0.05, edges ≥ 10, p < 0.05 for selection, p ≥ 0.20 and anomalous ≥ 5 for sampling).
- The input takes the first and last time points of each dataset, ρ = 1/day.
- The injected noise is the same as in v2: independent overdispersion φ = 5, 20; region-wise sampling σ = 0.5, 1.0.

## New simulated data

The generation code is not changed (`validation/sparse_joint/multi_timepoint`, `series_methods`, `dev/realdata/CE/c2_sim.R`); only the landscapes, seeds, number of cells per time point and output directories change.

**e group:**
- 12 new landscapes: `nz_L01`–`nz_L12`, seed base 101000 (checked as never used); wavelength classes as before (3 each of 0.4, 0.8, 1.2, 1.6).
- The design of the four sets is as before:
  - e1: selection, passage bottleneck 200 / 2000 / 20000, interval 6 / 60 days, 4 or 8 time points;
  - e2: clone replacement, starting from day 240 or 480, interval 120 days, bottleneck 20000 / 2000;
  - e3: neutral, design as in e1;
  - e4: true missegregation rate 4 times the default (odd landscapes) or 16 times (even).
- 200, 300, 500 and 1000 cells drawn per time point, sampled from the same simulated population. Previously only 200 and 1000.
- In total e1 576, e2 384, e3 576, e4 576 datasets.

**C2:**
- 12 new landscapes: `nz_c2_L01`–`nz_c2_L12`, seed base 102000; 6 each of wavelengths 0.8 and 1.2, as before.
- The 6 settings as before (β_sim 1.5e-4, 5e-4, 1.5e-3 × bottleneck 300, 3000).
- 300, 500 and 1000 cells drawn at each of the two time points (previously 300 and 450), intervals 15 / 30 / 45 / 60 days.
- 864 datasets in total.

**Trial run first:** first run all steps through on 1 landscape and time them, then submit everything.

## Criteria (only datasets in which both time points reach threshold T; T = 300 and 500 judged separately)

| No. | Criterion |
|---|---|
| 1 | Neutral data without drift (e3, bottleneck 20000, no injection): judged sampling-dominated ≤ 10% |
| 2 | Ordering by drift (e3, no injection): the fraction judged sampling-dominated is bottleneck 200 > 2000 > 20000, and the difference between 200 and 20000 is ≥ 0.30 |
| 3 | Simulations with selection (e1, e2, e4, C2, no injection, pooled): judged sampling-dominated ≤ 10% |
| 4 | Power (e3, bottleneck 20000, with independent overdispersion φ = 20 added): judged sampling-dominated ≥ 70% |

**How the threshold is chosen:**
- if all four are met at T = 300, the enabling threshold is set to 300;
- otherwise, if all are met at T = 500, it is set to 500;
- if neither, it does not pass.

**Report only:**
- the fractions judged selection and sampling under each bottleneck for region-wise sampling (σ = 0.5, 1.0);
- the distribution of verdicts at 200 cells;
- the results for φ = 5;
- the results for e1 and C2 with injected noise;
- the results for each C2 setting separately.

**What noise is injected into:** all four kinds on e3; φ = 20 and σ = 1.0 on e1 and C2.

## Afterwards

- **Pass:** apply it to the Salehi and NUMBAT datasets in which both time points reach the threshold, and report the verdict distribution and whether parallel lineages agree; then write a separate P5 plan (external estimate of φ + negative binomial).
- **Not pass:** write the report and stop.

## Rules

- Code and results in `dev/realdata/noise/v3/` and `results/realdata/noise/v3/`, not committed; the package is not changed.
- HPC: xxlarge, add-hours monitor with the small QOS, 3 days; extra arguments separated by colons.
- After writing `report_v3.md`, stop and wait for the user.

## Stop rules

- A generation or analysis script fails, and the cause is unclear;
- > 5% of the tasks in an array fail, and it is not a resource problem;
- a situation the rules do not cover.

## Corrections
- **2026-11-20 13:07 EST (before any data were generated): the seed bases are changed to 111000 (e group) and 112000 (C2).**
  - Reason: the 101000 and 102000 written above are not unused; `validation/sparse_joint/cloud_level/setup.R` uses 101000–103000 (and also 104000). My earlier check command failed on a zsh glob error and was never actually executed, yet I treated it as checked.
  - After re-checking, the seed bases already used in the repository are 61000, 63000, 64000, 71000, 82000, 93000–99000, 100000, 101000–104000; 110000–129000 do not appear.
  - The derived seeds (landscape_seed × 1000 + …, then × 10 + …) are all within integer range, and do not overlap with those derived from other bases.
