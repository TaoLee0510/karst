# The 35 parental training windows of Fig. 2c: fit domain and ploidy count census

2026-11-24, benchmark round `6ae2a24`. Taking as the entry point the 44 daughters frozen in the original full-parent version of Fig. 2c, and counting each parental training window only once, there are 35 windows in total. **In 16 windows the fraction of out-of-domain cells rises, and in 16 windows the high-copy fraction within the out-of-domain group rises; both rise together in 11.** These 11 windows include 6 of the SA906 family, 1 of the SA609 family and 4 of the other lineages. All conclusions describe the change between the two training time points. [Mapping of the 44 daughters to training windows](../results_6ae2a24/analysis/fig2c_two_timepoints/shared_full_parent/daughter_angles.tsv) · [Per-window data table](fig2c_ploidy_census_6ae2a24.tsv)

## Grouped results

| Family | Windows | Corresponding daughters | Out-of-domain of all ↑ | High-copy of out-of-domain ↑ | Both together ↑ | Among windows with out-of-domain rising, both together ↑ |
|---|---:|---:|---:|---:|---:|---:|
| SA906 family | 10 | 10 | 7 | 7 | 6 | 6/7 |
| SA609 family | 11 | 17 | 1 | 2 | 1 | 1/1 |
| Other lineages | 14 | 17 | 8 | 7 | 4 | 4/8 |
| Total | 35 | 44 | 16 | 16 | 11 | 11/16 |

- In the SA609 family, the only window where both rise is **SA609 X3 → SA609R2 X6**. The fraction of out-of-domain cells goes from 76.02% to 76.67%, and high-copy of out-of-domain from 2.31% to 5.49%. Family membership is determined by the source in the metadata; the SA000 branch is also assigned to SA609.
- Of the 4 windows in the other lineages where both rise, SA532 contributes 1 and SA535 contributes 3. The 3 windows of SA1035 have zero out-of-domain high-copy counts at both end points; their within-group high-copy fractions are all recorded as unchanged.
- The two rising sets both contain exactly 16 windows. Besides the 11 in the intersection, another 5 have only the out-of-domain fraction rising (of which the high-copy fraction within the out-of-domain group falls 3 times and is unchanged 2 times), and another 5 have only the high-copy fraction within the out-of-domain group rising. All 35 windows have out-of-domain cells at both end points, so all denominators are defined.

## Per-window cell counts and fractions

Each window uses `summary$karyotype` of its archived KARST fit object to determine the direct fit domain, and both ends share this domain. High copy is defined as a sum of the integer copy numbers of the 22 autosomes ≥66; low copy means a sum <66. The sums of the four classes of raw cell counts recover the total cell counts of the two time points respectively. Out-of-domain of all = (out-of-domain low + out-of-domain high)/all; high-copy of out-of-domain = out-of-domain high/(out-of-domain low + out-of-domain high). A rise is judged by cross-multiplying the raw integer numerators and denominators, and unchanged is listed separately. The percentages below are shown to two decimal places. “SA535-T, SA535-U” abbreviate `SA535_CISPLATIN_CombinedT` and `SA535_CISPLATIN_CombinedU` in the panel respectively.

### SA906 family

| Parent ID | Training start → training end | Start in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | End in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | Out-of-domain of all, start→end | High-copy of out-of-domain, start→end |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 | SA906a X10 → SA906a X15 | 777 | 25 | 12 | 9 | 750 | 6 | 5 | 6 | 2.55% → 1.43% | 42.86% → 54.55% |
| 16 | SA906a X10 → SA906a X25 | 779 | 25 | 10 | 9 | 196 | 2 | 5 | 9 | 2.31% → 6.60% | 47.37% → 64.29% |
| 9 | SA906a X15 → SA906a X30 | 750 | 0 | 5 | 12 | 538 | 45 | 10 | 97 | 2.22% → 15.51% | 70.59% → 90.65% |
| 11 | SA906a X30 → SA906a X40 | 538 | 51 | 10 | 91 | 253 | 97 | 6 | 40 | 14.64% → 11.62% | 90.10% → 86.96% |
| 12 | SA906a X30 → SA906a X50 | 538 | 51 | 10 | 91 | 134 | 199 | 7 | 49 | 14.64% → 14.40% | 90.10% → 87.50% |
| 15 | SA906a X10 → SA906b X20 | 779 | 25 | 10 | 9 | 857 | 16 | 23 | 21 | 2.31% → 4.80% | 47.37% → 47.73% |
| 14 | SA906b X20 → SA906b X30 | 861 | 0 | 19 | 37 | 878 | 0 | 30 | 41 | 6.11% → 7.48% | 66.07% → 57.75% |
| 8 | SA906b X20 → SA906b X35 | 861 | 0 | 19 | 37 | 1084 | 0 | 26 | 88 | 6.11% → 9.52% | 66.07% → 77.19% |
| 7 | SA906a X10 → SA906b X40 | 779 | 25 | 10 | 9 | 499 | 2 | 28 | 34 | 2.31% → 11.01% | 47.37% → 54.84% |
| 17 | SA906b X30 → SA906b X45 | 875 | 0 | 33 | 41 | 539 | 0 | 14 | 49 | 7.80% → 10.47% | 55.41% → 77.78% |

### SA609 family

| Parent ID | Training start → training end | Start in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | End in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | Out-of-domain of all, start→end | High-copy of out-of-domain, start→end |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 66 | SA609 X1 → SA000 X5 | 30 | 0 | 323 | 36 | 82 | 0 | 373 | 16 | 92.29% → 82.59% | 10.03% → 4.11% |
| 74 | SA609 X1 → SA000 X6 | 30 | 0 | 323 | 36 | 80 | 0 | 367 | 1 | 92.29% → 82.14% | 10.03% → 0.27% |
| 59 | SA609 X1 → SA609 X4 | 30 | 0 | 323 | 36 | 124 | 0 | 165 | 5 | 92.29% → 57.82% | 10.03% → 2.94% |
| 63 | SA609 X1 → SA609 X5 | 30 | 0 | 323 | 36 | 134 | 0 | 177 | 2 | 92.29% → 57.19% | 10.03% → 1.12% |
| 69 | SA609 X1 → SA609 X6 | 31 | 0 | 322 | 36 | 181 | 0 | 176 | 2 | 92.03% → 49.58% | 10.06% → 1.12% |
| 76 | SA609 X5 → SA609 X7 | 140 | 0 | 171 | 2 | 209 | 0 | 103 | 1 | 55.27% → 33.23% | 1.16% → 0.96% |
| 87 | SA609 X1 → SA609 X8 | 30 | 0 | 323 | 36 | 33 | 0 | 32 | 1 | 92.29% → 50.00% | 10.03% → 3.03% |
| 88 | SA609 X4 → SA609 X9 | 125 | 0 | 164 | 5 | 366 | 0 | 136 | 3 | 57.48% → 27.52% | 2.96% → 2.16% |
| 65 | SA609 X1 → SA609R2 X5 | 30 | 0 | 323 | 36 | 69 | 0 | 365 | 29 | 92.29% → 85.10% | 10.03% → 7.36% |
| 70 | SA609 X3 → SA609R2 X6 | 41 | 0 | 127 | 3 | 122 | 0 | 379 | 22 | 76.02% → 76.67% | 2.31% → 5.49% |
| 75 | SA609 X1 → SA609UnBU X6 | 30 | 0 | 323 | 36 | 154 | 0 | 160 | 20 | 92.29% → 53.89% | 10.03% → 11.11% |

### Other lineages

| Parent ID | Training start → training end | Start in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | End in-domain low | In-domain high | Out-of-domain low | Out-of-domain high | Out-of-domain of all, start→end | High-copy of out-of-domain, start→end |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 23 | SA1035T X5 → SA1035T X6 | 728 | 0 | 71 | 0 | 373 | 0 | 37 | 0 | 8.89% → 9.02% | 0.00% → 0.00% |
| 25 | SA1035U X4 → SA1035T X7 | 515 | 0 | 14 | 0 | 78 | 0 | 7 | 0 | 2.65% → 8.24% | 0.00% → 0.00% |
| 27 | SA1035U X6 → SA1035U X7 | 308 | 0 | 31 | 0 | 99 | 0 | 1 | 0 | 9.14% → 1.00% | 0.00% → 0.00% |
| 32 | SA532 X1 → SA532 X2 | 92 | 0 | 62 | 11 | 94 | 0 | 71 | 6 | 44.24% → 45.03% | 15.07% → 7.79% |
| 34 | SA532 X1 → SA532 X3 | 93 | 0 | 61 | 11 | 494 | 0 | 74 | 8 | 43.64% → 14.24% | 15.28% → 9.76% |
| 35 | SA532 X1 → SA532 X4 | 93 | 0 | 61 | 11 | 356 | 0 | 46 | 15 | 43.64% → 14.63% | 15.28% → 24.59% |
| 36 | SA532 X2 → SA532 X5 | 90 | 0 | 75 | 6 | 192 | 0 | 37 | 3 | 47.37% → 17.24% | 7.41% → 7.50% |
| 37 | SA532 X5 → SA532 X7 | 190 | 0 | 39 | 3 | 86 | 0 | 38 | 1 | 18.10% → 31.20% | 7.14% → 2.56% |
| 38 | SA532 X3 → SA532 X8 | 481 | 0 | 87 | 8 | 205 | 0 | 38 | 7 | 16.49% → 18.00% | 8.42% → 15.56% |
| 45 | SA535-T X6 → SA535-T X7 | 143 | 21 | 300 | 8 | 106 | 4 | 198 | 15 | 65.25% → 65.94% | 2.60% → 7.04% |
| 48 | SA535-U X5 → SA535-T X8 | 201 | 0 | 41 | 0 | 103 | 4 | 396 | 11 | 16.94% → 79.18% | 0.00% → 2.70% |
| 52 | SA535-T X6 → SA535-T X9 | 142 | 21 | 301 | 8 | 59 | 8 | 81 | 2 | 65.47% → 55.33% | 2.59% → 2.41% |
| 47 | SA535-U X5 → SA535-U X7 | 201 | 0 | 41 | 0 | 357 | 18 | 69 | 2 | 16.94% → 15.92% | 0.00% → 2.82% |
| 50 | SA535-U X5 → SA535-U X8 | 199 | 0 | 43 | 0 | 158 | 28 | 107 | 12 | 17.77% → 39.02% | 0.00% → 10.08% |

## Why the counts of SA906a X25 differ from the earlier singleton analysis

The training window of the fixed panel that predicts SA906a X30 is **X10 → X25**. In this domain, X25 has 196 low-copy and 2 high-copy cells inside the domain, and 5 low-copy and 9 high-copy cells outside the domain. The number of out-of-domain high-copy cells is 9 at both ends; the number of out-of-domain cells falls from 19 to 14, so the high-copy fraction within the out-of-domain group rises from 47.37% to 64.29%. The rise describes a share; the absolute count did not increase.

The earlier study devoted to the 11 singletons used a different fit, **X15 → X25**. It put all 11 high-copy cells of X25 outside the domain, and high-copy of out-of-domain fell from 12/16 to 11/16. In both analyses the raw total of X25 is 212 and the high-copy total is 11. With a different training start and fit domain, the resulting allocation inside and outside the domain, and the direction of change, also differ. These 35 windows count only the X10 → X25 selected by the frozen panel. [X15 → X25 singleton analysis](sa906a_alfak_karst_high_copy_6ae2a24.md) · [X15 → X25 raw numerators and denominators](karst_ploidy_prior_review_plan.md)

## Review and interpretation

- The training start, end, family and archived fit of each parent in the fixed panel are carried over; the actual domain is read with the corresponding `ends/karst/minobs_*/karst_fit.Rds`. Only the first and last columns of each training file are taken; the complete raw counts are checked row by row against the archive, and then cross-validated against the total cell count and in-domain cell count of each time point in the fit diagnostics. All 35 windows are kept, including the 6 windows the original fit marked as not converged; domain membership and raw counts can both be checked directly.
- The per-window data table contains the full time-point names, the training interval, the mapping of the 44 daughters to parents, the four classes of counts, unrounded fractions and percentage-point changes. Approximate 95% sampling intervals of the changes are attached as extra columns, combined from Wilson score intervals of two independent proportions (Newcombe method), conditional on the current fit domain. The intervals reflect only cell sampling, and do not cover clonal correlation, passage bottlenecks or uncertainty in the choice of domain; the grouped summaries are not filtered by interval.
- The 11 windows where both rise can be used to check ploidy candidate models; some windows share training time points, and a few changes are very small. The direction of the counts by itself is not enough to identify a ploidy-specific fitness advantage.
- Running `Rscript --vanilla benchmark/paper/scripts/fig2c_ploidy_census.R` from the root directory of the KARST repository rebuilds the per-window data table. The environment variable `BENCH_RESULTS` points to the round directory containing the archived fits, and `PAPER_DATA` points to the ALFA-K data directory. This time only a count census of existing inputs and fit objects was done. [Recomputation script](../scripts/fig2c_ploidy_census.R) · [Per-window data table](fig2c_ploidy_census_6ae2a24.tsv)
