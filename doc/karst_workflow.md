# KARST workflow, step by step

This is how KARST fits a landscape, from a two-time-point count matrix
to a fitted landscape and, optionally, to the output files of alfak. It is
Section 4 of [`karst_method.tex`](karst_method.tex) in plain text; the model
and its equations are Section 2 there, and the design decisions and benchmark
results are in
[`validation/sparse_joint/README.md`](../validation/sparse_joint/README.md).

## Entry points

All three run the same model, the sparse joint count model.

| Entry point | Use | Returns |
| --- | --- | --- |
| `fit_karst(counts, dt, beta)` | Main entry point | An `karst_sparse_joint_fit` |
| `fit_sparse_joint(counts, dt, beta, ...)` | The same model with its arguments given directly | The same |
| `alfak_compat(yi, outdir, ...)` | The arguments of `alfakR::alfak()`; writes alfak's files | The cross-validation R-squared, invisibly, and the files in `outdir` |

`fit_karst()` passes `counts`, `dt` and `beta` (the last two only when
given, so that a prepared `karst_data` object keeps its own rate) and the
list `sparse_joint_control` on to `fit_sparse_joint()`; giving an argument
both ways, or any other argument, is an error.

## Part 1. Fitting (`fit_sparse_joint`)

1. **Prepare the input.** A two-column matrix of integer counts, one column
   per time point, with karyotype labels such as `"2.2.3.2"` as row names:
   unsigned integers separated by single dots. Labels are rewritten without
   leading zeros (`"2.02.3.2"` is `"2.2.3.2"`); two rows that name the same
   karyotype, or a malformed label (`"2.2.3.2."`, `" 2.2.3.2"`), are an error.
   Rows with no cells are dropped. `dt` is the time between the two time points
   (days for the GBM inputs) and `beta` the missegregation probability per
   chromosome copy and division, which must be positive. Copy numbers must lie
   between `min_cn` (default 1) and `max_cn` (default one above the largest
   observed), and both time points need cells.

2. **Find the major karyotypes and the domain.**
   - A karyotype is *major* if it holds at least 1% of the cells at either time
     point and at least 3 cells over both. Without any, the most abundant one is used.
   - The *domain* is the majors and every one-step neighbour of a major (one
     copy more or fewer of one chromosome), observed or not. An unobserved
     neighbour enters the likelihood because its absence next to an abundant
     clone is evidence that it is less fit.
   - Observed karyotypes more than one step from every major are left out of
     the likelihood together with their cells. They would enter only because
     they were observed, while their unobserved siblings would not, which
     biases their fitness upwards. They are predicted later at the level of
     the observed non-major karyotypes.
   - If the domain would exceed `max_states` (5000), the share required of a
     major is raised until it fits.

3. **Build the division kernel.** For every pair of neighbouring domain
   states, the probability that a daughter of one has the other's karyotype,
   from the exact single-chromosome missegregation kernel (the same as
   alfak's `pij`), and for every state the probability `P_ii` that a
   daughter keeps the parent's karyotype.

4. **Build the trend design.** The centre is the count-weighted mean
   karyotype. Each state's covariates are an intercept, each chromosome's
   deviation from the centre (dosage), its square (curvature), and the
   distance to the nearest major raised to the power 0.85 (the shell term).

5. **Set the prior on initial abundances.**
   - Seen at the first time point: centred on its count, essentially free.
   - Seen only at the second: centred on its frequency there.
   - Never seen: a frequency near 1e-5 with a wide spread.
   - Default `initial_state = "balance"`: every non-major state is centred
     instead on its mutation-selection balance with the majors next to it,
     what they would have fed it over a history `h`, with spread
     `balance_sd = 1` on the log scale. Its count at the first time point then
     informs its fitness. Majors stay free.
   - By default `h` is chosen by marginal likelihood among the interval and
     30, 100, 300, 1000 and 3000 divisions at `birth_rate` that are longer
     than the interval: the model (steps 6-8) is fitted at each and the fit
     with the lowest objective at the posterior mode it reports (step 8) is
     kept, among those whose scale optimisation gave a finite objective
     (`diagnostics$history_profile`: `objective`, nlminb's, and
     `posterior_objective`, the one compared; up to 1.1.1 the candidates were
     compared on nlminb's). A candidate at which the posterior mode cannot be
     found (step 8) is not kept; with none left the fit stops. A given
     `history` fixes it.

6. **Assemble the count model (compiled with TMB).**
   - Latent variables: per state, the log expected count at the first time point
     (zeta) and the log birth rate (g); the trend coefficients; the depth
     ratio (theta).
   - Fitness prior: g is centred on the trend with a residual scale; the
     intercept on the log of the supplied `birth_rate` (sd 0.05); dosage and
     curvature coefficients have their own scales; the shell coefficient has
     prior sd 0.5.
   - Net growth rate: d = (2 P_ii - 1) exp(g).
   - Dynamics in closed form: at the second time point each state has grown at
     its own rate and received influx from its one-step neighbours, each
     source growing at its own rate (first order in missegregation).
   - Likelihood: the counts at each time point are Poisson, with expectation
     exp(zeta) at the first and the dynamics, scaled by theta, at the second.

7. **Estimate the scale hyperparameters.** The dosage, curvature and residual
   scales maximise the Laplace-approximated marginal likelihood (`nlminb`
   within [1e-4, 0.5]). For each trial value an inner Newton solve finds the
   posterior mode of all latent variables. The point `nlminb` returns is
   checked (a healthy inner Hessian, the same objective from its own mode, no
   better inner mode from the initial values, a small finite-difference
   gradient): `diagnostics$convergence` is 0 when it is accepted, and 1 with a
   warning when not. If `nlminb` fails, the prior-centred scales (about 1% per
   step) are used, and `convergence` is 1 with the reason "nlminb failed" in
   `diagnostics$optimisation`. If the inner solve fails at the scales kept,
   also from the initial values (step 8), that candidate history is not
   kept; with one history, or none left, the fit stops with an error.

8. **Compute the posterior.** The posterior mode at those scales (TMB's
   Newton search from its best mode so far; if that fails, once more from the
   initial latent values; if that fails too, an error of class
   `karst_inner_failure`) and the objective there
   (`diagnostics$posterior_objective`), with the covariance taken as the inverse of the
   *expected* information: the Hessian
   with the counts replaced by their fitted expectations, which is positive
   definite where the observed Hessian may not be. It is factorised once
   (sparse Cholesky).

9. **Report relative fitness and its uncertainty.**
   - `fitness_mean`: net growth rate per unit of `dt`; its absolute level
     follows the assumed birth rate.
   - `relative_fitness`: minus the count-weighted mean over the majors, which
     is what two time points identify.
   - `relative_sd_model`: the model's own posterior standard deviation, too
     small on its own.
   - `relative_sd`: calibrated,
     `sqrt(model_sd^2 + (deficit * relative)^2 + (0.009 * birth_rate * (1 + shell))^2)`,
     with `deficit` 0.6 for a karyotype the likelihood has a cell of and 1.1
     for one it has none of (with more than two time points and
     `intermediate = "majors"`, a non-major counts as seen only at the first
     or the last time point). The constants were fitted on the development
     landscapes, and the unobserved one on simulated regimes where the
     predicted drop below the majors was too shallow; the 95% intervals then
     covered at least 90% of the karyotypes in every simulated regime tested.

10. **Store the fit.** `summary` (one row per domain state), `trend`
    (coefficients and covariance, centre, majors, scales, observed level,
    reference, calibration), `settings` and `posterior` (all settings and
    the posterior mode, from which the posterior can be rebuilt) and
    `diagnostics`.

## Part 2. Prediction (`predict(fit, karyotypes)`)

| Karyotype | Prediction |
| --- | --- |
| In the domain | Its posterior from the fit |
| Observed but outside the domain | The level of the observed non-major karyotypes |
| Anything else, e.g. two steps away | The trend, with the residual scale added to its variance |

Every row has `relative_fitness`, `relative_sd` and `relative_sd_model`.

## Part 3. alfak outputs (`alfak_compat` = fit + `as_alfak` + `write_alfak_outputs`)

11. **Read alfak's input.** `yi$x` has two or more columns, one per time
    point. The times are `passage_times` or `as.numeric(colnames(yi$x)) *
    yi$dt`, as in alfak, and the model of Part 1 is fitted to the intervals
    between them (more than two time points: see Variants) with
    `beta = pm`; further arguments, such as `birth_rate`, go to
    `fit_sparse_joint()`. The time of the first time point sets time 0 of the
    reported frequencies (step 15).

12. **Choose karyotypes as alfak does.** `fq`: at least `minobs` cells over
    all time points of the input. `nn`: their one-step neighbours without nullisomy. The
    landscape covers `fq`, `nn` and the neighbours of `nn`, in alfak's
    order. `minobs` changes only this choice, not the fit.

13. **Put fitness on alfak's scale.** `mean` is relative fitness plus
    `g0 = log(nb / n0) / dt`, the growth rate alfak assigns to the population
    (`dt` the first interval). With `correct_efflux = TRUE` it is the birth
    rate as alfak converts it, `(relative - k) / (2 P_ii - 1)`, where `k`
    makes the mean over `fq`, weighted by their frequencies at time 0
    (step 15, at the posterior mode), equal `g0`. The same `k` is used for
    the draws, the cross-validation and `landscape_data`.

14. **Draw from the posterior in place of the bootstrap.** The posterior is
    rebuilt from the stored settings and mode and all latent variables are
    drawn jointly. The posterior is taken at the stored mode, with no new
    search (karst 1.1.3); `as_alfak()` records the gradient there in
    `rebuild`. Draws are mapped to relative fitness (domain states
    through g, trend predictions through the coefficients plus the residual,
    with the reference drawn alongside). Each karyotype receives independent
    noise for the part of its calibrated standard deviation the model does
    not explain. Draws come in mirrored pairs and each row is scaled to the
    calibrated standard deviation, so every row's mean and median equal the
    prediction and its standard deviation equals `relative_sd` (divided by
    `2 P_ii - 1` with `correct_efflux`), while the model's correlations
    between karyotypes are kept.

15. **Assemble the files.**
    - `landscape_posterior_samples.Rds`: one row per landscape karyotype,
      one column per draw.
    - `landscape.Rds`: `k`, `mean`, `median`, `sd`, `fq`, `nn`, the row
      summaries of the draws.
    - `bootstrap_res.Rds`: `final_fitness` (draws for `fq`), `nn_fitness`
      (draws for `nn`), frequencies among `fq` from the drawn abundances at
      the first time point, moved back to time 0 of alfak's time axis with
      each draw's fitness, so that `x0 * exp(f * t)` gives them back at the
      first time point; `initial_*` equal `final_*`; alfak's retry
      attributes.

16. **Cross-validate.** alfak's folds: each `fq` karyotype with its one-step
    neighbours, each karyotype in the first fold that contains it. For each
    fold, the counts of its domain states leave the likelihood; those whose
    initial abundance has the free prior (majors, or every state under
    `initial_state = "free"`) get the prior of a never-seen karyotype, while
    the others keep the balance prior, which does not use their counts; the
    trend's centre and majors are recomputed without them, and the posterior mode is found again with the scales
    fixed. The predictions, aligned to the full fit on the majors outside
    the fold, are compared with the full fit:
    `R2R = 1 - sum((p - o)^2) / sum(o^2)` with both series centred. This
    measures how well a karyotype the model has not seen is predicted, not
    how accurate the full fit is, and is not comparable with alfak's R2R.
    `xval.Rds` holds `R2R`, `predictions`, `observations`, `xval_data` and
    `tmp`.

17. **Write and return.** The four alfak files, `karst_fit.Rds`, and
    `landscape_data.Rds` when `landscape_data_output = TRUE`; R2R is
    returned invisibly.

## Variants

- `fit_sparse_joint_rates()`: when the birth rate or the missegregation rate
  is unknown, fits at a few rates across a range (`birth_rate_range` or
  `beta_range`; default three, evenly spaced on the log scale). Point
  predictions are at the lowest rate, intervals cover the whole range, and
  the agreement of rankings across rates is reported.
- `fit_sparse_joint_cohort()`: several two-time-point datasets at once. Each keeps its domain
  and initial abundances; a karyotype shares its fitness across datasets with
  an estimated deviation; the trend is shared or per dataset, whichever has
  the higher marginal likelihood, compared as the candidate histories are
  (step 5). Its default is `initial_state = "free"`.
- More than two time points of one population: `fit_sparse_joint(counts, times = ...)`
  with one column per time point. The later time points are compared with the same
  dynamics started at the first time point, each with its own depth. By default
  the majors and the domain come from the first and the last time point, and the
  time points in between inform the majors only, with a weight
  `1 / max(1, dispersion)` from a first fit (`diagnostics$intermediate_weight`).
  `as_alfak()` and `alfak_compat()` take such fits.
- Every fit checks its scale optimisation (`diagnostics$optimisation`) and
  warns when the optimum is not accepted, naming the karyotype along which the
  Laplace approximation is flat. The check does not change the fit.

## Running times

| Case | Time |
| --- | --- |
| One simulated dataset | Seconds to about 4 minutes |
| One GBM dataset at three birth rates | Median 149 s |
| `alfak_compat()` with cross-validation, GBM P6 (226 karyotypes) | About 8 minutes |
| Simulated culture, 4-8 time points of 200-1000 cells (default) | Median 2.5-5 minutes, 2-4 times two time points |

## Computation

The TMB models `sparse_joint`, `sparse_joint_series`, `sparse_joint_pairwise`,
`sparse_joint_state_space` and `sparse_joint_cohort` are compiled into one
library (`src/karst_tmb.cpp`) and selected by name. Parsing of karyotype labels
and the missegregation probability of one chromosome are separate Rcpp code;
the simulators are written in R.
