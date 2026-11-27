# O (overdispersion) paused under the stopping rule

Written on 2026-11-20 at 00:41 EST. None of O's decision fits have been submitted yet.

## What happened

The method for O fixed in the PLAN is: from the base fit take the states with expected count ≥ 5, compute the Pearson dispersion D₅ at the second time point, then take φ̂ = max(0, D₅ − 1).

After the patch `dev/patches/O_nb1_est.patch` was written and checked:
- with `overdispersion` not given or given as 0, it is bit-for-bit identical to cd342db;
- trying `"estimate"` on one e1 dataset (bottleneck 200, every 6 days, 1000 cells, first and last of 8 time points): only 3 states have an expected count ≥ 5, D₅ = 0.07, φ̂ = 0.

The reason is that each major has its own free fitness parameter, the counts at the two time points are fitted almost completely by them, and the residuals are absorbed. Under this rule φ̂ is 0 on almost all datasets, and O is the same as base. This is a case of "some criterion is unreasonable", so the work stopped.

## Where B's dispersion comes from

On B's base fits, the Pearson statistic is split by karyotype class. The table shows dataset medians.

| Group | Total dispersion D | Share from majors | Share from non-majors with expectation 1–5 | Share from non-majors with expectation < 1 that did appear | D of majors | D for expectation 1–5 | D for expectation < 1, appeared |
|---|---|---|---|---|---|---|---|
| NUMBAT pooled | 1.16 | 0.27 | 0.05 | 0.53 | 0.52 | 0.51 | 6.0 |
| Salehi confirmation | 4.18 | 0.03 | 0.02 | **0.95** | 0.50 | 0.48 | **6.5** |
| Salehi development | 1.65 | 0.06 | 0.07 | **0.87** | 0.55 | 0.45 | 2.3 |

- The "overdispersion" flagged by L6 comes almost entirely from rare karyotypes whose expected count is below 1 but which did appear at the second time point, not from drift in abundant clones.
- B1's ratio of zero-count neighbours is 0.83: fewer kinds of rare karyotypes appear than Poisson expects, but when they appear they have more cells than expected. That is, rare karyotypes appear in clusters, like individual small clones.
- In KARST such a rare karyotype gets its fitness pushed up, which may be a false selection signal.

## Possible approaches (for the user to decide)

1. **O′: use NB only for non-major states.** φ is estimated from the states with "expectation < 5", and majors still use Poisson. This targets exactly the clustered rare karyotypes. The method and the criteria of O would have to be rewritten.
2. **Do not do O.** R5 of phase 3 (random rare-karyotype cloud) has already failed on simulations, and it too targeted rare karyotypes.
