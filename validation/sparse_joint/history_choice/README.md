# History choice: which objective picks the candidate history

Study of 2026-11-19, asked for by the user before the selection rule of
karst 1.1.1 is changed.

## Question

A fit with the balance prior fits each candidate history and keeps the one
with the lowest objective. Up to karst 1.1.1 that objective is nlminb's: the
Laplace objective at the inner mode the scale optimisation ended with.

The posterior the fit reports sits at another mode, the one the last inner
search found:
- usually from TMB's best mode so far;
- from the initial latent values when that search failed (a candidate
  "rescued" by the second search of 1.1.1).

The objective at that mode can differ from nlminb's. Three rules are compared
on the same candidate fits:

| Rule | Value compared |
| --- | --- |
| current | nlminb's objective for every candidate (karst 1.1.1) |
| A | the objective at the posterior mode for a rescued candidate, nlminb's for the others |
| B | the objective at the posterior mode for every candidate |
| oracle | the candidate closest to the truth; for scale only |

## Design

Each candidate is fitted once. `capture.R` wraps `fit_sparse_joint()` and
`sparse_joint_inner_mode()` in the karst namespace of the R process, so that
the real fit keeps every nested single-history fit and the objective at its
posterior mode. The package is not changed. Each task checks that the kept
candidate's record is the fit returned. Every rule is then applied after the
fact (`analyse.R`).

### Development benchmark (two time points)

- The 168 scoreboard datasets (96 core, 72 stress), fitted as
  `karst_sparse_joint_balance`, the package default (`dev_task.R`).
- Accuracy uses the benchmark's own scoring at distances 0, 1 and 2: RMSE,
  Spearman correlation and 95% coverage against the true growth rates.

### ALFA-K paper benchmark, ABM and missegregation simulations

Covers the anchored, pairwise and state-space methods (`abm_task.R`). The
tasks of round 50415cf's tables:
- the 71 ABM and missegregation fit tasks that failed in round 2c60af1 with
  `karst_inner_failure`. Rescued candidates occur only in these;
- 306 others, 34 drawn at random (seed 1) from each design × method stratum,
  from every memory tier.

Accuracy is the Pearson and Spearman correlation of KARST's relative fitness
with the true fitness of the simulation's Gaussian random field, over the
domain states, the majors, and shells 1 and 2.

## Criteria, agreed before the results

On the fits whose choice rule B changes, B is adopted if:
- **Development benchmark:**
  - the median RMSE ratio B / current is at most 1.00 over distances 0–2 together, and at most 1.02 at each distance;
  - the median 95% coverage is at most 0.01 lower.
- **ABM, random set:**
  - the median difference of Spearman correlations (B − current, all domain states) is at least 0;
  - B is not significantly worse by a sign test (one-sided p ≥ 0.05).

If B is clearly worse, rule A is taken. In between, the numbers go to the user.

## Running it (HPC)

```
R=validation/sparse_joint/history_choice
Rscript $R/make_tasks.R KARST_DIR 50415cf benchmark/paper/build/tmp_inner_failures.tsv $R/results
$R/hpc_submit.sh            # one array per memory tier, add-hours monitors
Rscript $R/analyse.R $R/results
```

- Tasks run in the paper benchmark's container with the karst library of its
  lock (`hpc_task.sh`).
- Outputs stay on HPC in `results/` (not tracked); `report.md` is copied
  here when the study is done.

## Results

Run on HPC on 2026-11-19 with karst 1.1.1 (50415cf); full tables in
[report.md](report.md).

- **Every fit ran.** 168 of 168 development datasets and 377 of 377 ABM tasks
  were fitted. In every one, the kept candidate's record was the fit returned
  and the current rule reproduced the fit's own choice.
- **The posterior-mode value often differs from nlminb's, but rarely changes
  the choice.** It differs by more than 1e-3 at 22.5% of development
  candidates and 19.4% of ABM candidates. A candidate rescued by the second
  search sits at a worse mode (median +8.95 on the ABM). Even so:

  | Fits | Rule A changes | Rule B changes |
  | --- | --- | --- |
  | development (168) | 2 | 5 |
  | ABM random (306) | 0 | 3 |
  | ABM failed in 2c60af1 (71) | 3 | 6 |

  A and B never choose a rescued candidate. The current rule chose one in 2
  development datasets and 3 ABM tasks.
- **Where the choice changes, the effect is small and mixed:**
  - development: A's 2 datasets are better (mean RMSE ratio 0.947 and
    0.999); of B's 5, these 2, one more better (0.994) and 2 equal;
  - ABM random set: B's 3 tasks are better (Spearman +0.002 to +0.010);
  - ABM failed set: A is better at 1 of 3 and worse at 2 (−0.015, −0.016);
    B is better at 1 of 6 and worse at 5 (median −0.012, sign test
    p = 0.11).

  Pooled over each set, the three rules agree to the third decimal.
- **B meets all four agreed criteria:**
  - development: median RMSE ratio 1.000 overall, at most 1.000 at each
    distance; coverage unchanged;
  - ABM random set: median Spearman difference +0.005, 3 better and 0 worse.

  The failed set was not part of the criteria; there B leans worse, without
  significance.
- **The objective is not where the accuracy is lost.** The candidate closest
  to the truth differs from the current choice at 115 of 168 development
  datasets (median RMSE ratio 0.98–0.99) and 206 of 306 random ABM tasks
  (median Spearman +0.008). That gap
  is far larger than any difference between the rules, and none of the rules
  closes it.

Decision (2026-11-19): rule B, in karst 1.1.2. Candidate histories and the
two trends of a cohort are compared on `diagnostics$posterior_objective`.
These scripts were written for karst 1.1.1, whose own choice is the current
rule. Under 1.1.2 the check that the current rule reproduces the fit's own
choice would fail where B differs.
