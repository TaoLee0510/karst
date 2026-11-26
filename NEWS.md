# karst 1.1.4 (2026-11-23)

A second target for the cross-validation of `as_alfak()` and
`alfak_compat()`, so that their scores can be put on alfak's scale. The
default is unchanged: fits, landscapes, draws and scores are the same as in
1.1.3.

## New options

- `as_alfak(xval_target = "own")` compares each held-out prediction with the
  karyotype's fitness from its own counts, the counterpart of alfak's target
  (one bootstrap draw of its first-stage fitness):
  - every other latent variable stays at the fit's mode;
  - the karyotype's own terms of the objective are minimised over its birth
    rate and its own abundances, without the landscape prior on its birth
    rate. These are its counts, its initial abundance term and, in the
    state-space method, its own transitions.
- `xval_resample` (default `TRUE`) first resamples the counts once, cell by
  cell at each time point, as alfak's bootstrap resamples cells; `xval_seed`
  (default 1) seeds it, and the global random state is left as it was.
- Only karyotypes with at least one cell in the counts the fit uses are then
  tested. `xval` also holds `target = "own"`.
- `alfak_compat()` passes `xval_target` on.
- The default, `xval_target = "model"`, compares with the model's fitness from
  all the data, as before. That target comes from the same smooth model as
  the prediction, so the default scores are high and not on alfak's scale.

## Checks

- At the fit's mode, the R version of each karyotype's own terms (new file
  `R/sparse-joint-own.R`, all four methods) equals the TMB objective's to
  3e-13 (`tests/testthat/test-sparse-joint-own.R`).
- In round 6ae2a24 of the ALFA-K paper benchmark, 213 fits rescored with the
  full `as_alfak(xval_target = "own")`, refitting every fold, gave the same
  scores as the post-hoc computation (difference 0)
  (`benchmark/paper/docs/xv_target_6ae2a24_Ver2.md`).

# karst 1.1.3 (2026-11-20)

`as_alfak()` takes a fit's posterior at the mode the fit reports. Fits are
unchanged. `as_alfak()` and `alfak_compat()` change only in their draws, and
only for fits at whose stored mode TMB's search failed or moved on.

## Changed results

- `as_alfak()` rebuilds a fit's posterior to draw from it. Up to 1.1.2 the
  rebuild searched for the mode again, starting at the stored one, and in two
  cases did not end there:
  - Where the inner Hessian at the stored mode is numerically singular, TMB's
    search does not move but fails, because its Laplace value is not finite.
    1.1.1 and 1.1.2 then searched from the initial latent values and could
    end at another, worse mode; 1.1.0 stopped with an error.
  - Where the fit's own search had stopped short of the mode (the largest
    gradient over the latent variables up to 6.3 in round cd342db), the new
    search moved on.
- The draws of `bootstrap_res.Rds` and `landscape_posterior_samples.Rds` then
  came from the covariance at the other point: their correlations between
  karyotypes, and the spread of the drawn initial abundances. The landscape's
  `mean`, `median` and `sd` were not affected; they come from the fit.
- Now the fit's own posterior is taken at its stored mode, with no new search,
  as the fit computed it. Cross-validation folds, whose problem differs from
  the fit's, still search, from the stored mode and then from the initial
  values. So does a fit that did not store every latent variable.
- In round cd342db of the ALFA-K paper benchmark (5,985 fits rebuilt):
  - 7 rebuilds ended at another mode after TMB failed at the stored one. In
    the worst two the objective rose by 14.1 and 10.3.
  - 3 more moved by more than 1e-4 from a stored mode at which the fit's
    search had stopped early; the largest move was 0.40.

## New records

- `as_alfak()` returns `rebuild`, which `write_alfak_outputs()` does not
  write. `gradient` is the largest absolute gradient of the joint objective
  over the latent variables at the stored mode. `xval_folds` counts the
  cross-validation folds whose mode was found, of them those found only by the
  second search, and those left out because no mode was found.
- The solver of `sparse_joint_rebuild()` has `stored` and, for a fit's own
  posterior, `gradient`; its `value` is then the objective the fit reported.

# karst 1.1.2 (2026-11-19)

Candidate histories, and the two trends of a cohort, are compared on the
objective at the posterior mode each fit reports. A fit with a single history
is unchanged apart from the new record, unless it is a cohort fit that
chooses its trend (`trend = "select"`, the default).

## Changed results

- Up to 1.1.1 candidates were compared on nlminb's objective: the Laplace
  objective at the inner mode the scale optimisation ended with. The
  posterior a fit reports is at the mode found by the last search at the
  fitted scales (`sparse_joint_inner_mode()`). That mode can differ, and so can
  its objective: by more than 1e-3 for about a fifth of the candidates of the
  development benchmark and of the ABM simulations of the ALFA-K paper
  benchmark. A candidate whose mode was found only by the second search of
  1.1.1 is usually at a worse mode (median +8.95 on the ABM). It could still
  be chosen on nlminb's objective, and so was reported at a mode other than
  the one it was chosen for.
- Each candidate history is now compared on that objective. As before, only
  candidates whose scale optimisation gave a finite objective take part. With
  `trend = "select"`, `fit_sparse_joint_cohort()` compares its two trends the
  same way. The per-dataset histories of a cohort are chosen by
  `fit_sparse_joint()`, so they follow the new rule.
- Before the change, the effect was measured on 168 development datasets and
  377 ABM and missegregation fits (validation/sparse_joint/history_choice).
  The history changed at 5 and 9 of them, with no consistent change in
  accuracy: development RMSE ratio median 1.000; Spearman correlation on the
  random ABM set +0.005, 3 better and 0 worse.

## New records

- `diagnostics$posterior_objective` of every fit and cohort fit: the
  objective at the posterior mode it reports.
- `diagnostics$history_profile` has a column `posterior_objective`, the value
  compared. `objective` stays nlminb's.
- `diagnostics$trend_posterior_objective` of a cohort fit with
  `trend = "select"`: the value compared. `trend_objective` stays nlminb's.

# karst 1.1.1 (2026-11-19)

The posterior mode is searched once more before a fit fails. Fits in which
the search for the posterior mode at the fitted scales did not fail are
unchanged.

## Changed results

- When TMB's Newton search for the posterior mode at the fitted scales fails
  (the objective is not finite), it is made once more from the initial latent
  values, as the convergence check does, and TMB's own start is restored
  afterwards. Only when both fail is it an error (class
  `karst_inner_failure`). In a fit, TMB starts from its best mode so far,
  found during the scale optimisation, which can lie far from the mode at the
  fitted scales. When the posterior is rebuilt, including in each fold of
  the cross-validation of `as_alfak()`, TMB starts from the fit's stored mode
  and the second search from the problem's own initial latent values. 1.1.0
  left out a fold whose first search failed, so the cross-validation of such
  fits can change.
- A fit whose posterior mode was found only by the second search records the
  objective at that mode in `diagnostics$posterior_restart` (also a cohort
  fit). Its
  `objective`, by which candidate histories are compared, stays nlminb's, as
  in every earlier version. Its convergence check normally reports
  "objective not finite": the check's own evaluation starts where the first
  search did.
- A candidate history at which the posterior mode cannot be found is not
  chosen: in `diagnostics$history_profile` its objective and convergence are
  `NA`, and a column `inner_failure` marks it. A column `posterior_restart`
  gives the objective of each candidate whose mode needed the second search.
  Both columns appear only when needed. The fit stops only when no candidate
  is left, with class `karst_inner_failure` when every candidate failed so.
- In 1.1.0 one such candidate stopped the whole fit. This happened for 85 of
  the 6099 KARST fit tasks of the ALFA-K paper benchmark (round 2c60af1),
  84 of which had fitted in round e7e32f8. In the two of them checked, the
  failing candidate was not the one chosen, and the fit with the second
  search is the same as in round e7e32f8.
- When the posterior mode cannot be found, the first fit that sets the
  weight of the intermediate time points falls back to weight 1, as for a
  dispersion that is not finite. The check fit of `initial_state = "auto"`
  likewise keeps the balance, as for a level that is not finite.

# karst 1.1.0 (2026-11-18)

Fixes from the code review of 2026-11-18. Fits of valid input are unchanged
except where noted; `as_alfak()` and `alfak_compat()` with
`correct_efflux = TRUE` report different birth rates.

## Changed results

- `correct_efflux = TRUE` converts fitness to birth rates as alfak does
  (ALFA-K.R, dev2@f93b2d2): `(f - k) / (2 P_ii - 1)`, where `f` is relative
  fitness and the constant `k` makes the mean over the frequent karyotypes,
  weighted by their frequencies at time 0 (`final_frequencies` at the
  posterior mode), equal `log(nb / n0) / dt`. The landscape, the replicates,
  the cross-validation and `landscape_data` use the same `k`. Up to cfabe12
  the level was added first and the sum divided by `2 P_ii - 1`, which put
  that mean about `1 / (2 P_ii - 1)` times too high: 0.4% for 44 copies at a
  missegregation rate of 5e-5, 9% at 1e-3. With `n0 = nb = NULL` the model's
  own net growth rates are divided, as before.
- With more than two time points and `intermediate = "majors"`, a karyotype
  that is not major and was seen only at time points in between counts as
  unobserved both in the calibration of `relative_sd` (as before) and in the
  observed level that `predict()` gives observed karyotypes outside the
  domain (it counted as observed there).
- The pairwise and state-space cross-validation: a held-out state of the
  state-space method takes the drift weight of a state seen at no time
  point; its own counts set it before.
- The first fit that estimates `intermediate_weight` uses the shortest
  candidate history, as documented, also when `history` is not sorted.
- `simulate_sparse_counts()` draws each step's divisions with the mean that
  gives a karyotype its continuous-time growth over the step. The old mean,
  `count * b * step`, grew 5% slower than the birth rate `b` near 1 at the
  default step, and fitness differences came out 8-9% smaller than
  `truth_observed`. Simulations with the same seed differ.

## Errors where results were silently wrong

- `fit_karst()` passes `dt` and `beta` on only when given, so an
  `karst_data` object is fitted at its own missegregation rate; the default
  5e-5 replaced it before. Giving an argument both directly and in
  `sparse_joint_control` is an error.
- Karyotype labels must be unsigned integers separated by single dots.
  `prepare_karst_data()` and `predict()` rewrite them without leading zeros
  (`"2.02.3"` is `"2.2.3"`); two rows naming the same karyotype, or labels
  such as `"2.2.3."`, `" 2.2.3"` and `"+2.2.3"`, are an error. A trailing
  dot was read as a valid label and, matched as a string, made the karyotype
  a different one from its own counts.
- The missegregation rate of a fit must be positive. At zero the division
  kernel has no edges and the balance prior had no effect, while the fit
  still reported it.
- `balance_level = "free"` is an error for the pairwise and state-space
  methods, as for the anchored one; it was ignored there and stored as
  `"free"`.
- A failed inner optimisation when the posterior is computed or rebuilt is an
  error (class `karst_inner_failure`). The last point tried was used as the
  mode; in the cross-validation such a fold now fails instead of returning a
  point near the full-data mode.
- A fit with several candidate histories none of which has a finite
  objective stops with an error; it returned an object without a class.
- `fit_sparse_joint_cohort()`: `history` must give one time or one per
  dataset.

## Time points left out

- `summary` names its count columns after the time points of the input:
  `count_t0`, `count_t2` when `t1` was left out (they were numbered by
  position, so `count_t1` held `t2`).
- Karyotypes seen only at time points left out are in `dropped`, and
  `predict()` places them at the observed level; they were predicted from the
  trend as never seen.
- For the pairwise and state-space methods, `diagnostics$dropped_samples`
  has the same columns whether or not a time point was left out.

## Smaller fixes

- `predict()`, also for cohort fits, accepts a factor of labels and
  repeated karyotypes, also outside the domain, and reads labels as the fit
  writes them.
- `predict()` of `landscape_data` warns when it rounds non-integer copy
  numbers.
- `simulate_grf_landscape()`: a copy of a landscape whose parameters are
  changed no longer reads the cached fitness of the original;
  `centroid_min == centroid_max` and a karyotype with a single copy no
  longer hit `sample()` on a vector of length one.
- `simulate_sparse_counts()`: `detection_threshold = 0` no longer keeps
  karyotypes without a cell; an `overdispersion` so small that the
  Dirichlet proportions underflow is an error, and the argument is
  documented as the Dirichlet concentration it is (smaller values are more
  overdispersed).
- `print()` of a pairwise or state-space fit names the pair of each coverage
  row.
- Documentation: `alfak_compat()` follows alfak branch dev2 (f93b2d2), whose
  `xval.Rds` is a list; `major_count`, `max_states`, `dt`/`times` for
  prepared data, the state-space initial prior, `intermediate` and
  `intermediate_weight` with the other methods, and the pairwise expected
  counts are described as the code behaves.

## Development benchmark (`benchmark/sparse_two_timepoint`)

- The acceptance rule is an RMSE below the null predictor and a Spearman
  correlation above alfak (`minobs = 20`) at every distance. A criterion
  whose baseline was not compared on every landscape of the split is
  reported as not evaluated, not as passed: held-out had no alfak
  `minobs = 5` fits, and that criterion was skipped silently.
- A failed fit's rank correlation is 0, the null predictor's, as the method
  document states; it was left out of the mean (held-out, wavelength 1.6:
  alfak's Spearman 0.334 / 0.294 / 0.277 becomes 0.257 / 0.227 / 0.215 at
  distances 0 / 1 / 2).
- Rescoring merges with the other methods' rows by column instead of
  dropping them; the comparison cohorts keep the members of the original 12
  landscapes; the document generators read the wavelength classes from the
  data; a missing fit is NA in the direction accuracy, not a tie.
- `forward_prediction.R` and `rate_sensitivity.R` fit with the method's
  arguments (`initial_state = "free"` for `karst_sparse_joint`) instead of
  the package default.

## ALFA-K paper benchmark (`benchmark/paper`)

- Fig. 5b fits KARST on one minobs per trajectory, as Fig. 4 does: its three
  minobs come from one fit, and each karyotype counted up to three times. On
  the tables of round e7e32f8 the WGD interaction of KARST with all time
  points turns opposite to the paper's.
- A Salehi fit fails by its run status, not by a missing cross-validation
  score, in every table: 66 fits that ran and had forecasts but no CV score
  were left out of the paired forecasts (two time points: 371 pairs become
  394).
- The agreement with the paper's conclusions leaves out intercepts, compares
  only significance when neither side is significant, and describes Fig. 6
  instead of counting it (16 statements per arm).
- Fig. 4d records the share of frequent karyotypes the paper's `var(deltaf)`
  drops, and an `na.rm` version beside the paper's.
- The removed `karst_auto` no longer labels a fit; collecting checks every
  fit of the task tables and stops on a missing one; files are written under
  a temporary name and renamed; only complete forecasts are scored; a rerun
  of a task fits only the minobs not done.

## Validation scripts (`validation/sparse_joint`)

- The Salehi horizon test bins later passages by days, not by levels of the
  metadata's parent links, which in the p53 lines are 5 to 21 passages apart;
  the pooled results are unchanged. The comparison of the two designs on a
  common support needs the forecasts to be run again.
- The truth of the beta-choice and series-method studies is the net growth
  rate the model estimates, not the birth rate (factor 0.93 at 16 times and
  0.74 at 64 times the default missegregation rate); their results are marked
  as not recomputed.
- Several scripts stop instead of falling back silently: a missing build
  library, a missing score file, an expected error that no longer occurs.
- Numbers that need a rerun to change are marked "not recomputed with the
  current code" in the READMEs.

# karst 1.0.0 (2026-11-17)

alfak2 is renamed KARST (Karyotype Abundance-based Reconstruction of fitness
landscape from Sparse Time points). The methods are those of alfak2 0.2.0,
with the fixes below; names change, and old names have no aliases.

- The package is `karst`: `library(karst)`. Functions, classes, options and
  condition classes that carried `alfak2` carry `karst`: `fit_karst()`,
  `prepare_karst_data()`, `karst_data_times()`; classes such as
  `karst_sparse_joint_fit`; options such as
  `karst.sparse_joint_template`; conditions such as
  `karst_domain_coverage_error`.
- alfak, the R package `alfakR` of the ALFA-K paper, is called alfak in the
  documentation, labels and KARST's own function names: `as_alfak()` (was
  `as_alfakR()`), `write_alfak_outputs()` (was `write_alfakR_outputs()`),
  class `karst_alfak_output`; `alfak_compat()` keeps its name. Versions of alfak
  are named by branch and commit, such as alfak dev2@f93b2d2.
- Commit ids cited from before the rename map to the current history through
  `doc/commit_map.tsv`.
- The ALFA-K paper benchmark is part of the repository (`benchmark/paper`).

## Removed: the two-layer and joint-count pipelines

KARST fits the sparse joint count model only. The pipelines it replaced, kept
since 0.2.0 as options of `fit_karst()`, are removed with their code, tests
and validation studies; they remain in the git history.

- `fit_karst()` passes `counts`, `dt`, `beta` and `sparse_joint_control` to
  `fit_sparse_joint()`. Its `pipeline` argument and the arguments of the other
  two pipelines are gone; giving one is an error that says so.
- Removed functions: `fit_local_posterior()`, `build_karyotype_graph()`,
  `fit_extrapolation_layer()` with its second-layer methods,
  `fit_graph_posterior()`, `summarize_karst()` and `plot_karst()`, and the
  classes `karst_fit`, `karst_local_fit`, `karst_global_fit` and
  `karst_graph` with their methods. With them go the effective-depth
  preprocessing, the alfak-scale columns (`alfak_scale`) and the count
  attributes `observation_weights`, `soft_minobs` and `holdout_mode` that
  `prepare_karst_data()` kept for the local fit.
- Removed compiled code: the TMB local model, the graph Gaussian posterior,
  the graph builder, the transition operator, and a count simulator that no R
  function had called since May 2026. The TMB source is now `src/karst_tmb.cpp`.
- Removed validation studies of these pipelines: `validation/accuracy_round2`,
  `accuracy_upgrade`, `baseline_quadrature`, `observation_design` and
  `posterior_reference_4to8`.
- The development benchmark no longer fits the two-layer pipeline, and its
  acceptance rule no longer compares with it; the method document no longer
  describes or reports either pipeline.
- Sparse joint fits are unchanged: two- and multi-time-point, cohort,
  rate-range and alfak-format fits give identical results before and after
  the removal.

## Validation

Every KARST fit of the ALFA-K paper benchmark, the Salehi horizon test and the
development benchmark (held-out landscapes included) was redone with 1.0.0
(round e7e32f8, `benchmark/paper/docs/report_e7e32f8.md`):

- The fits are those of alfak2 0.2.0 (round 05a8855): the same landscapes
  (to 1e-16), histories, cross-validation predictions and forecasts, except
  the cross-validation of pairwise fits, which the fix above changes.
- The paper's tests, the Fig. 2d noise check and the forecasts against later
  passages give the numbers of the alfak2 era. Where the paper picks fits by the
  cross-validation score, the benchmark now resolves ties (scores within 1e-8,
  the lowest minobs kept): KARST derives every minobs from one fit, and the
  paper's rule counted many inputs two or three times, as often as rounding
  made their scores equal.
- On the development benchmark, coverage at one and two steps is higher than in
  the runs of 2026-11-10/11, whose fits predate the deficit of karyotypes seen
  at no time point; the method document's numbers are from this rerun.

## Fixes (code review of 2026-11-17)

- A missegregation rate given to `fit_sparse_joint()` or
  `fit_sparse_joint_cohort()` replaces the one stored in `karst_data`; before,
  `fit_sparse_joint_rates(beta_range =)` refitted prepared data at its stored
  rate only.
- `as_alfak()`: the cross-validation of a `series = "pairwise"` fit centres the
  unseen karyotypes of each pair on the depth at the start of the pair, as the
  fit does, not on the first time point's depth.

# karst 0.2.0 (2026-11-17, released as alfak2)

The first release built around the sparse joint count model: fitness
landscapes from sparse karyotype counts at two time points, the case the
package is meant for, with more time points as an extension.

## The default pipeline

- `fit_karst()` now fits the sparse joint count model (`fit_sparse_joint()`)
  by default. Both time points enter one count likelihood over the
  well-sampled ("major") karyotypes and all their one-step neighbours,
  observed or not, so that an unobserved neighbour of an abundant clone is
  evidence that it is less fit. The earlier two-layer pipeline remains
  available (`pipeline = "two_layer"`).
- Rare karyotypes start at mutation-selection balance with the clones next to
  them (`initial_state = "balance"`, the default), so that their level at the
  first time point informs their fitness. The time over which the balance
  built up (`history`) is chosen by marginal likelihood among the interval
  and 30 to 3000 divisions longer than it.
- Relative fitness comes with calibrated uncertainty (`relative_sd`), with a
  larger deficit for karyotypes seen at no time point.
- Fitness away from the domain is predicted from a fitted trend in chromosome
  dosage, curvature and distance to the major karyotypes (`predict()`).
- `diagnostics$domain_coverage` reports how much of each time point the domain
  holds; a time point with no cell in it stops the fit with an error of class
  `karst_domain_coverage_error`. `check_sparse_joint_counts()` summarises a
  count matrix before fitting.
- The scale optimisation is checked and a fit that is not accepted warns
  (`diagnostics$optimisation`).

## Extensions

- **More than two time points** (`times =`): `series = "anchored"` (default),
  `"pairwise"` or `"state_space"`. Two time points give the same fit with all
  three. Time points in between that carry no information under the model are
  left out with a warning (`diagnostics$dropped_samples`).
- **Cohorts**: `fit_sparse_joint_cohort()` fits several datasets together,
  with the trend shared or per dataset by marginal likelihood.
- **Uncertain rates**: `fit_sparse_joint_rates()` fits over a range of
  baseline birth rates or missegregation rates (`beta_range`), with point
  predictions at the lowest rate of the range.
- **alfak outputs**: `as_alfak()`, `write_alfak_outputs()` and
  `alfak_compat()` produce the objects and files of `alfakR::alfak()` (branch
  dev), with alfak's karyotype sets, scale, bootstrap layout and
  cross-validation folds.
- `initial_state = "free"` estimates every initial abundance independently;
  `initial_state = "auto"` (two time points) chooses between the two by the
  size of the cloud of rare karyotypes. `"auto"` misfires when the time points
  are one passage apart and is not the default.

## Validation

Scripts and reports are in `validation/sparse_joint/`; the method and its
results are in `doc/karst_method.tex`.

- On the development benchmark of simulated tumours (sparse two-time-point
  scenarios, 180–1080 days) and its held-out landscapes, the sparse joint
  model lowered the error of predicting no fitness differences to 0.47–0.74
  of that value, where alfak (minobs 20) exceeded it.
- On the ALFA-K paper's own tests (benchmark round of KARST 348addd, against
  alfak dev2@f93b2d2), KARST's landscapes were more accurate on the paper's
  simulations, and its forecasts beat the baseline more often.
- On the Salehi passages, forecasts scored against every later passage by a
  proper score (log-likelihood relative to no change): KARST beat alfak by
  1.01 per cell with two time points and by 1.12 with every passage
  (`validation/sparse_joint/salehi_horizon`). The paper's Fig. 2d ranking on
  these data is dominated by sampling noise and does not separate the two.

## Known limitations

- A missegregation rate other than the assumed one, or miscalled cells,
  cannot be told from the counts of two time points; fit over a range of
  rates when the rate is uncertain.
- Without a known division rate, rankings are robust but per-day magnitudes
  are conditional on the assumed rate.
- Intervals are too narrow after 1080 days and when the landscape changes at
  the first time point; unobserved neighbours are predicted slightly too fit.
- See "Open" in `validation/sparse_joint/README.md` for the full list.
