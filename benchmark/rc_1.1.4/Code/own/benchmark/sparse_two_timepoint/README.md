# Sparse two-timepoint benchmark

Scores fitness-landscape estimates from sparse karyotype counts at two time
points against a known truth, separately for observed karyotypes and their one-
and two-step neighbours. It is the development scoreboard for KARST: every
candidate change is compared here against the original ALFA-K (`alfak`) and
the no-information answer before it is considered for the 22-chromosome HPC
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
- **100 held-out landscapes per wavelength** (2026-11-20), as the ALFA-K paper
  simulates 100 landscapes per wavelength. The first three of each class are
  L01-L15 as before (identifiers, seeds and data unchanged); the other 97 of
  each class are L16-L500, class by class (0.4, 0.8, 1.2, 1.6, 0.2), with
  seeds 82016-82500. The development split keeps its fifteen landscapes.
- **Nullisomy is lethal.** Daughters with a copy number outside 1–8 die. The
  ALFA-K paper's simulations (alfakR 4ad82b0, "fix 0 CN bug") also kill a
  daughter with a zero-copy chromosome; none of their 525 simulations has one.
  alfakR versions before 4ad82b0 let nullisomic cells live: their populations
  fix a zero-copy chromosome, which never occurs in the real inputs and which
  `alfak` cannot fit because it discards every nullisomic neighbour. The upper
  limit of 8 is this benchmark's own; the paper's simulator has none (up to 8
  copies in its main set, 10 at the highest missegregation rates).

Truth is the Gaussian-random-field landscape of the ALFA-K simulations, with
the same formula and generating settings (10 centroids with coordinates −10
to 20, a diploid founder with fitness above 0.4), and the simulator
(`abm.cpp`) keeps their division and missegregation rules. It differs from the
paper's simulator in three rules: the upper copy-number limit above; a clone
with fitness ≤ 0 does not divide but survives until passaging thins it out,
where the paper's simulator kills it at once (this changes no dataset at
wavelengths ≥ 0.4 and little at 0.2: `doc/notes/nonpositive_fitness.md`); and
the number of dividing cells is binomial rather than a capped Poisson. The
KARST report (`doc/report`) compares the two simulation set-ups item by item.
`run.R --stage=validate` checks it against closed-form results: the daughter
distribution against the missegregation kernel, viability losses, growth,
selection between two clones, passaging, and reproducibility.

A landscape has 8 "patients": independent trajectories from the diploid founder
on the same landscape. Patient 1 makes the single-dataset scoreboard (96
development datasets); all eight form the cohort for cross-dataset learning.
The control cohort of a landscape (`other_landscapes`) is patient 1 of the
next seven, cyclically, of the twelve landscapes of the four original
wavelength classes of its split, so that landscapes added later do not change it.

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
neighbours with copy numbers 1–8. The one-step shell of a near-monoclonal
dataset can be smaller (33 of the 96 held-out datasets, all 24 at wavelength
0.4, and 36 of the 96 development datasets have fewer than 2000 one-step
neighbours; at least 118); the two-step shell always has more. All methods
are scored on the same nodes; a node a method cannot predict is scored as "no
fitness difference" and counted in `prediction_fraction`. A failed fit, or one
without a prediction for every reference karyotype, is scored as the
no-information answer on the whole dataset: its error is that answer's, and
its Spearman correlation is zero, the rank correlation of a prediction that
orders nothing (`stt_rank_correlation()`; `stt_null_rank()` in `report.R`
applies this to metrics scored before 2026-11-18).

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

1. RMSE below the no-information answer (95% interval of the paired
   difference below zero);
2. Spearman correlation above the original ALFA-K with `minobs = 20`
   (`acceptance$alfak`, the setting fitted on both splits; interval above
   zero).

Coverage is reported separately, with 90–98% as the target. The bar for error
is the no-information answer rather than ALFA-K because ALFA-K's RMSE is
dominated by a few unconstrained fits and is far above it. `run.R
--stage=report` evaluates the rule for the candidates named in the
configuration (`acceptance$candidates`), and reports only the methods the
configuration lists. A criterion is evaluated only when the candidate and its
baseline are compared on every landscape of the split; one whose baseline was
not fitted, or was fitted on fewer landscapes, is reported as "not evaluated",
never counts as passed, and keeps the candidate from passing the rule. Other
ALFA-K settings (`acceptance$original`) appear in the paired tables only.
Cohort methods are compared, on the cohort scenario, with the single-dataset
model of the same initial state (`acceptance$cohort`).

## Running

Both packages must be installed, or `STT_R_LIB` must point to a library that
holds the builds under test.

```sh
R=benchmark/sparse_two_timepoint/run.R
Rscript $R --stage=validate
Rscript $R --stage=truth    --split=dev --cohort
Rscript $R --stage=datasets --split=dev --cohort
Rscript $R --stage=fit      --split=dev --methods=zero,alfak_minobs20,karst_sparse_joint_balance
Rscript $R --stage=score    --split=dev
Rscript $R --stage=report   --split=dev
```

Results go to `benchmark/results/sparse_two_timepoint/` (not tracked). Fits are
cached per dataset and method, so a stage can be rerun after adding a method.
`--stage=score --methods=a,b` replaces those methods' rows of
`metrics_<split>.csv` and keeps the others'; if the scoring has gained or
lost a column since, the kept rows get NA there and a warning says so.

`--workers` defaults to 2. `alfak` needs up to 4.5 GB and 15 minutes on the
largest datasets (`minobs = 5` at 3000 cells), so raise the worker count only
with memory to spare.
The held-out landscapes are for the final confirmation only and need
`--split=heldout --unseal-heldout`.

### A round in a folder of its own

`STT_OUTPUT_DIR` names the results folder that `stt_config()` uses when no
folder is given (every script here and `validation/sparse_joint/far_distance.R`;
the HPC image receives it from `../paper/hpc/env.sh`), so that a round can be
run beside the earlier ones. `round_tasks.R` writes the task table of the fits a
round still needs, listed from the configuration (the datasets need not exist
yet) and leaving out fits already in the folder:

```sh
export STT_OUTPUT_DIR=benchmark/results/sparse_two_timepoint_<round>
T=benchmark/sparse_two_timepoint/round_tasks.R
Rscript $T --split=heldout --methods=karst_sparse_joint_balance,karst_sparse_joint --roles=core,stress karst.tsv
Rscript $T --split=heldout --methods=alfak_minobs10,alfak_minobs20 --landscapes=16-500 alfak.tsv
Rscript $T --split=heldout --methods=karst_sparse_joint_balance --sweep=far far.tsv
```

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

## Reproducibility

The simulator draws its random numbers while iterating over a
`std::unordered_map` of clones and through `std::binomial_distribution`, both
of which the C++ standard leaves to the implementation. The same seed
therefore gives different populations with different compilers or standard
libraries (the Mac's libc++ against the HPC's libstdc++). The generated truth
populations and datasets, not the seeds, are the record: generate them once
and copy them. Every fit starts from its dataset's seed (`run.R --stage=fit`,
`validation/sparse_joint/far_distance.R`), so `alfak`'s bootstrap does not
depend on the worker; sparse joint fits draw no random numbers.

## Limits

The simulation has no spatial sampling bias and no karyotype-calling error,
and treatment is represented only by the bottleneck scenario. Passing this
benchmark is necessary for a candidate, not evidence of accuracy on patient
data.
