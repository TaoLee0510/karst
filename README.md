# KARST

**KARST: Karyotype Abundance-based Reconstruction of fitness landscape from
Sparse Time points.**
The R package `karst` infers fitness landscapes of aneuploid karyotypes from
sparse single-cell karyotype counts at two or more time points, and supports
analyses of karyotype evolution.

KARST was developed and released as `alfak2` up to version 0.2.0. The history
of this repository includes that work; commit ids that the documents and
result folders cite from before the history was cleaned of tool attribution
lines are mapped to the current ones in [doc/commit_map.tsv](doc/commit_map.tsv).
Here *alfak* is the R package `alfakR` of the ALFA-K paper; its versions are
named by branch and commit: alfak main@4ad82b0 (the version the paper's fits
were made with, inferred from their dates) and alfak dev2@f93b2d2 (the version
the benchmarks fit).

```r
remotes::install_github("TaoLee0510/karst")
```

The model, fitted by `fit_karst()` or `fit_sparse_joint()`, is the sparse
joint count model, described [below](#sparse-joint-model-for-sparse-counts-at-two-time-points).
`alfak_compat()` runs it with the arguments of `alfakR::alfak()` and writes
alfak's output files, so that analyses built on alfak read its results
unchanged.

## Method Components

- One count likelihood for the well-sampled karyotypes and all their one-step
  neighbours, observed or not, at two or more time points
  (`fit_sparse_joint()`).
- Initial abundances of rare karyotypes tied to mutation-selection balance
  with the clones they descend from (`initial_state = "balance"`, the
  default).
- Fitness of any karyotype outside the fitted domain from a fitted trend
  (`predict()`).
- Several datasets fitted together (`fit_sparse_joint_cohort()`), and fits over
  a range of birth or missegregation rates (`fit_sparse_joint_rates()`).
- A check of the counts before fitting (`check_sparse_joint_counts()`).
- Outputs in alfak's format (`alfak_compat()`, `as_alfak()`,
  `write_alfak_outputs()`).
- Benchmark truth landscapes can be generated as lightweight paper-style GRF
  oracles with `simulate_grf_landscape()`. Fitness is computed on demand for
  requested karyotypes, so synthetic simulations can use 22 chromosomes without
  enumerating the full copy-number lattice.

## Minimal Example

```r
library(karst)

landscape <- simulate_grf_landscape(
  n_chr = 4,
  n_centroids = 30,
  lambda = 0.8,
  min_cn = 1,
  max_cn = 4,
  seed = 1
)

sim <- simulate_sparse_counts(
  landscape,
  beta = 0.00005,
  dt = 1,
  n0 = 200,
  n1 = 200,
  seed = 2,
  initial_population = 1000
)

fit <- fit_karst(
  sim$counts,
  dt = 1,
  beta = 0.00005,
  sparse_joint_control = list(min_cn = 1, max_cn = 4)
)

head(fit$summary)
predict(fit, c("2.2.2.2", "2.3.2.2"))
```

## Sparse joint model for sparse counts at two time points

`fit_sparse_joint()`, the model `fit_karst()` fits, is built for the
case the package is meant for: a few hundred cells at each of two time points. It fits both time points with one
count likelihood on the well-sampled karyotypes and all their one-step
neighbours, observed or not. A rare karyotype is held at mutation-selection
balance by the clone it descends from, so its level, not its change in
frequency, carries its fitness; an unobserved neighbour of an abundant clone is
evidence that the neighbour is less fit.

```r
fit <- fit_sparse_joint(counts, dt = 180, beta = 0.00005, birth_rate = 0.5)
fit$summary                                  # karyotypes of the fitted domain
predict(fit, c("2.2.3.2.2", "2.2.3.3.2"))    # any karyotype, also two steps away

# the same model through the main entry point
fit <- fit_karst(counts, dt = 180, beta = 0.00005,
                  sparse_joint_control = list(birth_rate = 0.5))
```

Two time points identify fitness differences, not absolute fitness. Use
`relative_fitness`, which is measured from the count-weighted mean of the
major karyotypes, and `relative_sd`. `relative_sd` adds two terms of model
error, fitted to simulated landscapes, to the model's own posterior standard
deviation (`relative_sd_model`), which is too small on its own. The error grows
with the predicted difference, and more for karyotypes seen at neither
time point: how far below the clones they lie is not identified by zero counts.
Because `relative_sd` grows with the difference itself, its 95% intervals
always contain zero. They are for coverage, not for testing whether a
karyotype differs from the reference. On simulated neutral cultures, about
5% of karyotypes had |`relative_fitness`| / `relative_sd_model` above 4.5
([series methods](validation/sparse_joint/series_methods/README.md)).

The baseline birth rate has to be supplied; two time points do not identify it.
Rankings hardly depend on it. Fitness differences between well-sampled clones
do not depend on it either, but the deficit of a rare or unobserved karyotype
is fixed per division, so per unit of time it is proportional to the assumed
rate. When only a range is known, fit across it. Point predictions are taken at the
lowest rate of the range: too low a rate shrinks deficits towards zero, while
one four times too high predicts worse than no differences at all.

```r
fit <- fit_sparse_joint_rates(counts, dt = 180, beta = 0.00005,
                              birth_rate_range = c(0.05, 1))
predict(fit, karyotypes)   # relative fitness at the lowest rate of the range
                           # (point = "lowest", the default; "centre" for the
                           # central rate), an interval covering the range, and
                           # fitness per division
```

The missegregation rate `beta` acts the same way: it sets how many cells
missegregation puts into a neighbour, so a neighbour seen less often reads as
less fit. It is fixed, by default at 5e-5, the rate of the ALFA-K paper.
Choosing it per dataset by marginal likelihood was tried and not adopted
([beta choice](validation/sparse_joint/beta_choice/README.md)).
- On the paper's simulations the choice worked.
- With time points every 6 days (18–42 days from first to last) it chose
  rates up to 16 times too high, and errors rose by up to 3 times.
- Gates on the interval and cells, evidence thresholds and a prior were
  tested on simulations whose true rate is higher. None kept the default's
  accuracy where the default was right and still found a higher rate. With
  few cells, short intervals or strong drift, two time points do not
  separate a higher rate from drift and sampling noise.
- A drift term would not help. With two time points, selection, influx and
  drift each explain the one change seen. With more time points, an estimated
  drift (`series = "state_space"`, below) read drift as selection on neutral
  cultures. The rate has to come from outside the counts, for instance from a
  measurement.

When the rate is uncertain, fit across a range of rates, with the most likely
rate at the lowest end. Too high a rate costs more than too low a rate:

```r
fit <- fit_sparse_joint_rates(counts, dt = 180, beta_range = c(5e-5, 8e-4))
predict(fit, karyotypes)   # fitness at 5e-5, and an interval covering the range
```

On simulated cultures whose true rate was 4 or 16 times the default, the
range 1–16 times restored the coverage of the 95% intervals, which had been
as low as 0.85. The point predictions were those at the default. The
intervals were wider, most of all when the default was right.

Several datasets, such as several patients or cultures with two time points
each, can be fitted together. Each keeps its own domain and initial
abundances; a karyotype that is abundant in one dataset and an unobserved
neighbour in another is informed by both, and its fitness may deviate between
datasets by an estimated amount. The trend is either shared by all datasets or
fitted per dataset, whichever has the higher marginal likelihood, so that
datasets from different landscapes are not forced onto one trend:

```r
cohort <- fit_sparse_joint_cohort(list(P1 = counts1, P2 = counts2),
                                  dt = c(180, 360), birth_rate = 0.5)
cohort$landscape                       # shared fitness of every karyotype
predict(cohort, karyotypes, sample = "P1")
```

With the balance prior (`initial_state = "balance"`) rare karyotypes start at
mutation-selection balance with the clones next to them, so that their levels
at the first time point inform their fitness as well. This takes fitness before
the first time point to resemble fitness between the time points. On simulated
landscapes it was more accurate than `initial_state = "free"`, which uses the
first time point only as a starting point, even when the landscape changed at the
first time point. The time over which the balance built up (`history`) is chosen
by marginal likelihood among the interval and longer spans. With time points only
about 30 divisions apart, its errors were 4-10% larger than those of the free
model, but it ranked karyotypes better; at about 60 divisions the two were
equally accurate. Up to KARST 8327282 the balance built up over the interval
alone, and at about 30 divisions its errors were twice those of the free
model.

The balance prior reads the size of the cloud of rare karyotypes as fitness.
That is right when the missegregation rate is the assumed one. When the true
rate is higher, or cells are called a copy off, the cloud is larger for
reasons other than fitness, and on such simulations the balance prior ranked
neighbouring karyotypes worse than chance; `initial_state = "free"` did not.
`initial_state = "auto"` (two time points) tries to tell the two apart: a fit
at the interval with the cloud's level free measures how much larger the cloud
is than the balance over the interval allows, and above a factor of e it fits
`"free"`. It is an option, not the default:
- With time points one passage apart (the ALFA-K paper's simulations) the
  cloud is far larger than the balance over one interval even at the true
  rate. `"auto"` switched there in most fits, and their landscapes were
  clearly less accurate (correlation with the truth 0.58 against 0.75). It
  was the default in KARST 827c4b5 only.
- Where the rate is uncertain or calls are noisy, and the interval spans many
  divisions, `"free"` or `"auto"` can be the safer fit
  ([cloud level](validation/sparse_joint/cloud_level/README.md)).

A fit needs cells at both time points within one step of a well-sampled
karyotype. `fit$diagnostics$domain_coverage` gives, for each time point, how many
of its cells lie inside the fitted domain. When one time point has none, typically
because the two time points share no karyotype and all well-sampled ones come from
the other, the fit stops with an error of class
`karst_domain_coverage_error`, which carries the same table.
`check_sparse_joint_counts()` summarises a count matrix before fitting: the
share of cells on repeated karyotypes and how far each cell lies from the
nearest other cell at its time point.

```r
check_sparse_joint_counts(counts)
fit$diagnostics$domain_coverage
```

Time points from a patient, such as a primary tumour and its recurrence, are not
passages: the recurrence may grow from a few surviving cells, and only for part
of the time between the two time points. A few survivors cost little. Growth that started late
costs a great deal if the model is not told: the fitness differences are
spread over the whole interval and come out too large per unit of time. Give
the time the population actually grew, or `birth_rate` times the share of the
interval during which it grew; when that is not known, `fit_sparse_joint_rates()`
over a range of rates is the safer choice, at a cost when there was no delay
([patient regrowth study](validation/sparse_joint/patient_regrowth/README.md)).

The scale parameters are fitted by Laplace marginal likelihood, and where a
karyotype has very few cells the Laplace approximation can be nearly flat or
jump between explanations of its counts. Every fit checks the optimum it
reports (`diagnostics$optimisation`: whether it was accepted, and if not why)
and warns when it is not accepted, naming the latent variables along which
the approximation is flat (for instance the initial abundance of a karyotype
seen once). The fit itself is not changed by the check; on simulated data,
refitting such fits until they were accepted did not change their accuracy
([convergence study](validation/sparse_joint/README.md)).

### More than two time points

`fit_sparse_joint()` also fits one population sampled more than twice: a count
matrix with one column per time point and `times` (or `dt`, one interval per pair
of consecutive time points). Every later time point is compared with the same
dynamics started at the first time point, each with its own depth, and fitness is
taken to be constant over the series. By default (`intermediate = "majors"`)
the majors and the domain come from the first and the last time point, as for two
time points, and the time points in between inform the majors only, with a weight
that falls with the dispersion of their counts (`intermediate_weight`), so a
population that drifts between time points, for instance through small passaging
bottlenecks, leans on the first and the last.

```r
fit <- fit_sparse_joint(counts, times = c(0, 6, 12, 18), beta = 0.00005, birth_rate = 0.9)
fit$diagnostics$intermediate_weight
```

On simulated passaged cultures (12 landscapes; 4 or 8 time points 6 or 60 days
apart; bottlenecks of 200, 2000 and 20000 cells; 200 or 1000 cells per time point)
the default was as accurate as fitting the first and the last time point alone:
median RMSE ratios of 0.97-1.01, within 3% in 21 of 24 settings, with the
same coverage. Under this model the time points in between add little. Using them for
every karyotype (`intermediate = "all"`) was 4-6% more accurate with little
drift (bottleneck 20000) but up to 2.4 times less accurate with strong drift
(bottleneck 200), where clones that rise and vanish between time points read as
fit; it suits large populations passaged through large bottlenecks only. The
default takes 2-4 times as long as two time points, `"all"` up to 20 times
([stage 2](validation/sparse_joint/multi_timepoint/README.md)). Cohort fits
and the two older pipelines take two time points.

`series` chooses how the time points are tied together. All three methods share
one landscape, and with two time points all three give the same fit.

| `series` | dynamics | majors and domain | time points in between |
| --- | --- | --- | --- |
| `"anchored"` (default) | one trajectory from the first time point | from the first and the last time point | as above |
| `"pairwise"` | each pair of consecutive time points from its earlier time point as observed | each pair's own | end of one pair, start of the next |
| `"state_space"` | latent abundances at every time point, linked by the dynamics plus a drift whose scale is estimated (`fit$diagnostics$drift`) | each pair's own | each counted once |

```r
fit <- fit_sparse_joint(counts, times = c(0, 120, 240, 360), beta = 0.00005, birth_rate = 0.9,
                        series = "state_space")
```

The three were compared on passaged cultures (as above), on clonal sweeps
(the first time point while clones still replace one another, time points 120 days
apart) and on neutral cultures. Results are the median RMSE ratio to the
first and last time point alone
([series methods](validation/sparse_joint/series_methods/README.md)).
- **Clonal sweeps**: `anchored` 0.89–0.94, `pairwise` 0.68–0.85,
  `state_space` 0.79–0.82.
- **`pairwise` under drift**: up to 1.8 with bottlenecks of 200 cells, and
  1.1–2.0 on neutral cultures. Each pair restarts from the counts observed at its
  earlier time point, so passaging drift within a pair reads as selection.
- **`state_space`**: 0.93–1.02 on the passaged cultures, but up to 1.2 at
  one and two steps on neutral cultures. It takes 3–12 times as long as
  `anchored`; with 8 time points some fits take up to 18 hours and 16 GB.

`anchored` stays the default. The other two suit clonal sweeps through large
bottlenecks.

An automatic choice among the three was tried and not adopted. It took
`pairwise` or `state_space` when a time point in between had less than 30% of
its cells on the majors of the first and the last time point (a clonal sweep),
and `anchored` otherwise. Its thresholds were tuned on the simulated sets
above. On the ALFA-K paper benchmark it was then less accurate than
`anchored`:
- on missegregation landscapes with 8 passages;
- on Salehi forecasts.

So the method is named in `series`
([series methods](validation/sparse_joint/series_methods/README.md); the
diagram of the whole fit is in [doc/flowcharts](doc/flowcharts/flow_overview.pdf)).

The first and the last time point must each have a cell in the domain. Some
time points in between have no cell on the states they count for: no cell on a
major by default, or none in the domain with `"all"`. A clonal sweep can do
this, when a passage in the middle lies two or more steps from both ends'
majors.
- Such a time point is left out with a warning, and the fit is the fit to the
  other time points. Its depth is free, so under the model it carries no
  information.
- `fit$diagnostics$dropped_samples` lists the time points left out.
- `as_alfak()` still picks its karyotypes from every time point, as alfak
  does.
- With `"pairwise"` and `"state_space"`, both time points of each consecutive
  pair need a cell in that pair's own domain. A time point in between without
  one is left out, and its neighbours become a pair. If the first or the last
  time point has none, its neighbour is left out instead, so the series keeps its
  span. At worst this gives the two-time-point fit of the ends.
- Results on the Salehi trajectories of the ALFA-K paper benchmark (187 with
  more than two passages that the default can fit):
  - the default left out 103 time points in 53 trajectories;
  - the pair rule left out one time point in each of 42;
  - no trajectory stops.

### Outputs in the format of alfak

`alfak_compat()` takes the arguments of `alfakR::alfak()` (branch `dev`), fits
the sparse joint model, writes the same files to `outdir` and returns the
cross-validation R-squared invisibly. A pipeline switches by calling it where
it called `alfak()`:

```r
r2r <- alfak_compat(yi, outdir, minobs = 20, nboot = 45, n0 = 1e5, nb = 1e7, pm = 0.00005,
             birth_rate = 0.05)   # further arguments go to fit_sparse_joint()
readRDS(file.path(outdir, "landscape.Rds"))   # k, mean, median, sd, fq, nn
```

| File | Content |
| --- | --- |
| `landscape.Rds` | Frequent karyotypes (`fq`, at least `minobs` cells), their one-step neighbours (`nn`) and the neighbours' neighbours, in alfak's order, with `mean`, `median` and `sd` of fitness |
| `landscape_posterior_samples.Rds` | One row per landscape karyotype, one column per draw: joint draws from the model's posterior, in place of alfak's Kriging bootstrap |
| `bootstrap_res.Rds` | `initial_fitness`, `final_fitness`, `initial_frequencies`, `final_frequencies` and `nn_fitness`, one row per draw, with alfak's retry attributes; the frequencies are at time 0 of alfak's time axis, so `x0 * exp(f * t)` gives them at time `t` as in alfak |
| `xval.Rds` | `R2R`, `predictions`, `observations` and `xval_data`, plus `tmp`, the paired matrix earlier alfak versions saved |
| `landscape_data.Rds` | With `landscape_data_output = TRUE`: `mean` and `median` models whose `predict()` gives the fitness of any karyotype |
| `karst_fit.Rds` | The KARST fit itself |

`yi$x` may have more than two columns, as in alfak; they are fitted as
described [above](#more-than-two-time-points).

The scale is alfak's: net growth rate per unit of time, with the
count-weighted mean of the major karyotypes set to `log(nb / n0) / dt`, `dt`
the first interval between time points as in alfak.
`correct_efflux = TRUE` reports birth rates. `minobs` only chooses which
karyotypes are listed as frequent; the fit uses every cell. The
cross-validation holds out each frequent karyotype with its neighbours, as in
alfak's folds, but by removing their counts and refitting the model, so its
R-squared measures how well KARST predicts a karyotype it has not seen.
`as_alfak()` and `write_alfak_outputs()` do the same for an existing fit.

[`doc/karst_workflow.md`](doc/karst_workflow.md) follows a fit
step by step, and [`doc/karst_method.tex`](doc/karst_method.tex) gives the
mathematics. Design, benchmark results and checks on real data are recorded in
[`validation/sparse_joint/README.md`](validation/sparse_joint/README.md); the
benchmark is [`benchmark/sparse_two_timepoint`](benchmark/sparse_two_timepoint/README.md).
