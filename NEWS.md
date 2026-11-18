# karst 1.0.0 (2026-11-17)

alfak2 is renamed KARST (Karyotype Abundance-based Reconstruction from Sparse
Time points). The methods are those of alfak2 0.2.0, with the fixes below;
names change, and old names have no aliases.

- The package is `karst`: `library(karst)`. Functions, classes, options and
  condition classes that carried `alfak2` carry `karst`: `fit_karst()`,
  `prepare_karst_data()`, `karst_data_times()`, `summarize_karst()`,
  `plot_karst()`; classes such as `karst_sparse_joint_fit`; options such as
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

## Validation

Every KARST fit of the ALFA-K paper benchmark, the Salehi horizon test and the
development benchmark (held-out landscapes included) was redone with 1.0.0
(round e7e32f8, `benchmark/paper/docs/report_e7e32f8.md`):

- The fits are those of alfak2 0.2.0 (round 05a8855): the same landscapes
  (to 1e-16), histories, cross-validation predictions and forecasts, except
  the cross-validation of pairwise fits, which the fix above changes.
- The paper's tests, the Fig. 2d noise check and the forecasts against later
  passages give the numbers of the alfak2 era. Tables that pick fits by the
  cross-validation score keep tied scores, as the paper does, and move slightly
  with rounding.
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
- The two-layer pipeline (`fit_local_posterior()`, `build_karyotype_graph()`)
  stops on more than two time points; before, it read the first two and
  rescaled fitness by their share of the whole interval.
- The default `anchor_covariance_inflation` covers every status a local fit
  reports; a status it lacked stopped the two-layer fit.
- `alfak_scale = TRUE` with a covariance reported for the near-field nodes only
  (`truncated_nearfield_gmrf`) no longer stops on mismatched dimensions.

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
