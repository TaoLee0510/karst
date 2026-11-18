# Sparse two-timepoint benchmark

Scores fitness-landscape estimates from sparse karyotype counts at two time
points against a known truth, separately for observed karyotypes and their one-
and two-step neighbours. It is the development scoreboard for KARST: every
candidate change is compared here against the original ALFA-K (`alfak`) and
the current KARST pipeline before it is considered for the 22-chromosome HPC
confirmation.

## Design

The simulated data are built to resemble real two-timepoint tumour datasets. The
reference is a set of 20 GBM inputs (22 chromosomes, copy numbers 1–4). Among
the 15 with at least 100 cells at both timepoints, and for the simulated
`typical_500_700_360d` datasets of the development landscapes:

| Statistic (median, quartiles) | Real GBM | Simulated |
| --- | --- | --- |
| Karyotypes observed | 65 (55–182) | 162 (45–328) |
| Fraction seen at both timepoints | 0.31 (0.12–0.53) | 0.13 (0.08–0.27) |
| Total-variation distance between time points | 0.82 (0.47–0.91) | 0.71 (0.10–0.85) |
| Largest clone | 0.13 (0.09–0.21) | 0.25 (0.13–0.81) |
| Karyotypes seen in one cell only | 0.27 (0.04–0.34) | 0.57 (0.53–0.61) |
| SD of log frequency change, shared karyotypes with ≥20 cells | 1.73 (0.98–2.25) | 2.12 (1.71–2.50) |

`run.R --stage=profile --real-dir=DIR` recomputes both columns. The simulated
data have more single-cell karyotypes and correspondingly fewer shared ones
than the real inputs; the other statistics overlap.

Choices that follow from this comparison:

- **22 chromosomes.** With 8 chromosomes the population is trapped on a local
  peak within about 700 days and the cells drawn at a time point are 98% one
  karyotype.
- **First time point at day 1440.** Earlier, clones replace each other
  completely within a year (no karyotype is shared between the time points);
  later, the smoother landscapes freeze as well.
- **Wavelengths 0.4, 0.8, 1.2, 1.6**, three landscapes each. The 0.4 class is
  the frozen, near-monoclonal regime that about one real dataset in ten
  resembles. The others give ongoing adaptation. A fifth class, **0.2**, the
  most rugged landscapes of the ALFA-K paper's simulations, was added on
  2026-11-18 after the others (landscapes L13-L15 of each split), so that the
  first twelve keep their identifiers and seeds. So far it is generated and
  fitted for the held-out split only.
- **Nullisomy is lethal.** Daughters with a copy number outside 1–8 die. The
  original ALFA-K simulator lets nullisomic cells live; its populations then
  fix a zero-copy chromosome, which never occurs in the real inputs and which
  `alfak` cannot fit because it discards every nullisomic neighbour.

Truth is the Gaussian-random-field landscape of the ALFA-K simulations, and the
simulator (`abm.cpp`) keeps their division and missegregation rules.
`run.R --stage=validate` checks it against closed-form results: the daughter
distribution against the missegregation kernel, viability losses, growth,
selection between two clones, passaging, and reproducibility.

A landscape has 8 "patients": independent trajectories from the diploid founder
on the same landscape. Patient 1 makes the single-dataset scoreboard (96
development datasets); all eight form the cohort for cross-dataset learning.

| Scenario | Cells at t0 / t1 | Interval (days) |
| --- | --- | --- |
| `very_sparse_50_50_180d` | 50 / 50 | 180 |
| `sparse_200_200_180d` | 200 / 200 | 180 |
| `sparse_200_200_360d` | 200 / 200 | 360 |
| `typical_500_700_360d` | 500 / 700 | 360 |
| `unbalanced_100_3000_360d` | 100 / 3000 | 360 |
| `unbalanced_3000_100_360d` | 3000 / 100 | 360 |
| `long_500_700_1080d` | 500 / 700 | 1080 |
| `bottleneck_500_700_360d` | 500 / 700, regrown from 200 surviving cells | 360 |

These eight core scenarios make the acceptance verdict. Six stress scenarios
break an assumption on purpose and are reported beside it:

| Stress scenario | Cells at t0 / t1 | Interval (days) | What it breaks |
| --- | --- | --- | --- |
| `short_500_700_30d`, `short_500_700_60d` | 500 / 700 | 30, 60 | Few generations between the time points (27 and 54), as for a slowly dividing tumour |
| `shift_500_700_360d`, `shift_3000_100_360d` | 500 / 700, 3000 / 100 | 360 | Half of the landscape's waves are redrawn at the first time point, as when treatment starts then; fitness before and after stays correlated |
| `turned_500_700_360d`, `turned_3000_100_360d` | 500 / 700, 3000 / 100 | 360 | The landscape is turned about the dominant karyotype at the first time point: that karyotype keeps its fitness, but which neighbour is fit is drawn anew |

In the last four, truth is the landscape between the time points.

## Scoring

Scored nodes are defined from the counts alone: every observed karyotype
(distance 0) and uniform samples of up to 2000 one-step and 2000 two-step
neighbours with copy numbers 1–8. All methods are scored on the same nodes; a
node a method cannot predict is scored as "no fitness difference" and counted
in `prediction_fraction`.

Truth is the per-day growth rate the simulated counts follow,
`log(1 + fitness * delta_t) / delta_t`. Truth and prediction are both centred
on the count-weighted mean of the karyotypes with at least 5 cells, so a
method's absolute scale does not enter.

Per distance: RMSE, `rmse_zero` (the RMSE of predicting no differences),
regression slope of prediction on truth, Spearman correlation, AUC for
"fitter than the reference", sign accuracy of one-step fitness differences
into the shell, and 95% interval coverage and interval score where a method
reports uncertainty.

Datasets are averaged within landscape and scenario, then landscapes are
weighted equally. Comparisons are paired by landscape with a landscape
bootstrap interval.

**Acceptance rule for a candidate.** Pooled over scenarios, at each of
distances 0, 1 and 2:

1. RMSE below the no-information answer and below the current KARST pipeline
   (95% interval of the paired difference below zero);
2. Spearman correlation above the current KARST pipeline and above the
   original ALFA-K (interval above zero);
3. no scenario with RMSE more than 10% above the current KARST pipeline.

Coverage is reported separately, with 90–98% as the target. The bar for error
is the no-information answer rather than ALFA-K because ALFA-K's RMSE is
dominated by a few unconstrained fits and is far above it, while the current
KARST pipeline sits on it. `run.R --stage=report` evaluates the rule for every
method that is not one of these baselines.

## Running

Both packages must be installed, or `STT_R_LIB` must point to a library that
holds the builds under test.

```sh
R=benchmark/sparse_two_timepoint/run.R
Rscript $R --stage=validate
Rscript $R --stage=truth    --split=dev --cohort
Rscript $R --stage=datasets --split=dev --cohort
Rscript $R --stage=fit      --split=dev --methods=zero,alfak_minobs5,karst_two_layer_fixed
Rscript $R --stage=score    --split=dev
Rscript $R --stage=report   --split=dev
```

Results go to `benchmark/results/sparse_two_timepoint/` (not tracked). Fits are
cached per dataset and method, so a stage can be rerun after adding a method.

`--workers` defaults to 2. A KARST two-layer fit at 22 chromosomes needs
about 1 GB. `alfak` needs up to 4.5 GB and 15 minutes on the largest datasets
(`minobs = 5` at 3000 cells), so raise the worker count only with memory to
spare. `karst_two_layer_cv`, the package's default leave-one-anchor-out
selection of graph hyperparameters, is defined but not part of a normal round:
one dataset with 97 observed karyotypes did not finish in 25 minutes, against
10 seconds with the hyperparameters fixed.
The held-out landscapes are for the final confirmation only and need
`--split=heldout --unseal-heldout`.

### Adding landscapes

`add_landscapes.R` makes the truth population, the datasets and the scored
nodes of named landscapes only, for a split whose other datasets exist (run.R
would need the truth populations of every landscape of the split), and then
adds them to the split's dataset index:

```sh
S=benchmark/sparse_two_timepoint/add_landscapes.R
L=--landscapes=heldout_L13,heldout_L14,heldout_L15
Rscript $S --split=heldout --unseal-heldout $L --task=1     # one landscape per task (or SLURM_ARRAY_TASK_ID)
Rscript $S --split=heldout --unseal-heldout $L --index-only # once all are made
```

### Rerunning the fits of an earlier run on the HPC

The fits of the KARST 1.0.0 rerun (2026-11-17, the alfak2 era's fits
refitted, held-out landscapes included) were run one fit per Slurm task in
the image of the paper benchmark (`../paper/hpc`):

```sh
# from the karst clone; FITS_FROM is the earlier run's fits/ folder
Rscript benchmark/sparse_two_timepoint/make_tasks.R FITS_FROM TASKS.tsv
# one fit per task: hpc/submit.sh dev TASKS.tsv 8G 06:00:00 (from benchmark/paper)
Rscript benchmark/sparse_two_timepoint/run_task.R TASKS.tsv ID
# afterwards: --stage=score and --stage=report of each split, and
Rscript benchmark/sparse_two_timepoint/compare_eras.R FITS OLD_FITS OUT.csv
```

`make_tasks.R` maps the alfak2 era's method names (`alfak2_*`, `alfakR_*`)
to the current ones, copies `alfak`'s and the null predictor's fits instead of
refitting them, and leaves out methods the configuration no longer has.
`compare_eras.R` compares each new fit with the earlier one (status,
predicted fitness and its sd). The distance 3-4 sweep
(`validation/sparse_joint/far_distance.R`) runs the same way, one fit per
task (`--ids=<dataset> --methods=<method> --fit-only=yes`; `hpc/submit.sh
far`), then once without these options to score every fit.

## Limits

The simulation has no spatial sampling bias and no karyotype-calling error,
and treatment is represented only by the bottleneck scenario. Passing this
benchmark is necessary for a candidate, not evidence of accuracy on patient
data.
