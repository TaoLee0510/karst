# Clones with fitness ≤ 0: the paper's simulator and the own benchmark's

2026-11-20; rerun on 2026-11-21 on round 6ae2a24. Scripts are in
`doc/notes/nonpositive_fitness/`, and the results of each round are in
`round_<round>/` there.

## Question

Both simulation benchmarks draw the true landscape from the same Gaussian
random field, f(k) = Σᵢ sin(‖k − cᵢ‖ / λ) / (π √10). This fitness can be
negative, and the two agent-based simulators treat a clone with f ≤ 0
differently:

- **ALFA-K paper** (alfakR 4ad82b0, `src/abm_core.cpp`, which made the paper's
  simulations): the whole clone dies in that step.
- **KARST's own benchmark** (`benchmark/sparse_two_timepoint/abm.cpp`): the
  clone does not divide, and its cells survive until passaging thins them out.

The note answers three questions:

1. How large is the effect?
2. Can it be measured?
3. Would the own benchmark's results change if its simulator followed the
   paper's rule?

## Rounds

The results below are those of round 6ae2a24: 500 held-out landscapes (100
per wavelength), fitted with KARST 1.1.3 and alfak `dev2` f93b2d2. They are
the data of the report's Part 3A.

The first run, on round e7e32f8, used 15 held-out landscapes (KARST 1.0.0).
Its paired simulation had 8 trajectories per landscape. It reached the same
conclusions; its results are in `round_e7e32f8/`.

## Mechanism

Under both rules a cell with f ≤ 0 never divides, so it leaves no descendants,
including descendants by missegregation. Its lineage ends either way. The rules
differ only in two respects:

- the cells occupy places in the population, which can move a passage forward
  by a step;
- the cells can be drawn into the data.

Both effects scale with the share of such cells in the population.

## Results

### 1. Prevalence in the own benchmark's data (`prevalence.R`)

The 7,000 scored held-out datasets: 500 landscapes × 14 designs (8 core, 6
stress). Medians per dataset, with maxima in brackets:

| λ | Cells with f ≤ 0 in the population at t1 | Observed cells with f ≤ 0 | Scored nodes one step away | Two steps away |
|---|---|---|---|---|
| 0.2 | 4×10⁻⁶ (4×10⁻⁴) | 0 (0.20, a "turned" stress design) | 12% (31%) | 30% (40%) |
| 0.4 | 0 | 0 (1×10⁻³) | 0 (4.5%) | 0.05% (8.1%) |
| 0.8 | 0 | 0 | 0 | 0 (0.15%) |
| 1.2, 1.6 | 0 | 0 | 0 | 0 |

- **More than 10⁻³ of a population** at f ≤ 0 occurs only at t0 of the
  "turned" stress designs, 12 datasets in all. Turning the landscape at t0
  puts about 20% of the t0 cells of landscape L14 at f ≤ 0 (less than 0.6%
  elsewhere). Under the benchmark's rule these cells are thinned to almost
  nothing by t1, 360 days later.
- **At t1** no population has more than 4×10⁻⁴ of its cells at f ≤ 0.
- **Core designs:** in four of them a single sampled cell has f ≤ 0.

### 2. The truth of a lethal karyotype (`rescore.R`, `rescore_summary.R`)

**The current convention.** The benchmark's truth is the per-day rate the
counts follow, log(1 + max(f, 0) Δt) / Δt. A clone with f ≤ 0 therefore gets
rate 0. That is below every positive rate, so for ranks (Spearman) these nodes
already rank last, tied.

**Two alternatives.** Under the paper's rule such a clone dies at once, which
has no finite rate. The existing fits of the 4,000 held-out core datasets (KARST
and alfak at minobs 5, 10 and 20) were rescored under two other readings:

- **death**: rate log(1 + f Δt) / Δt, that is, dying at rate |f|. This is the
  signed scale the paper scores on.
- **exclude**: these nodes are left out.

With the current convention (clip), the rescoring reproduces the benchmark's
stored scores. The largest differences over 48,000 rows are 5×10⁻¹³ in RMSE
and 6×10⁻¹⁶ in Spearman.

The table gives the error relative to the null predictor ("no fitness
differences") at λ = 0.2, pooled over designs and landscapes as in the
benchmark. Pooled over these datasets, 12% of the one-step and 29% of the
two-step nodes have f ≤ 0.

| | clip (current) | death | exclude |
|---|---|---|---|
| KARST, one step | 0.39 | 0.41 | 0.39 |
| KARST, two steps | 0.32 | 0.37 | 0.35 |
| alfak, minobs 5–20, one step | 3.2–4.0 | 3.1–3.9 | 3.4–4.2 |
| alfak, minobs 5–20, two steps | 2.9–3.6 | 2.7–3.3 | 3.1–4.0 |

- **Observed karyotypes (distance 0)** hardly change: KARST 0.537, 0.538 and
  0.536.
- **KARST's error** moves by at most 16% (two steps, death).
- **At λ = 0.2**, no Spearman correlation changes by more than 0.022.
- **At λ ≥ 0.4**, no ratio or Spearman correlation changes by more than 0.003.

### 3. The own benchmark's simulator with the paper's rule

**A plain switch cannot isolate the rule (`unpaired_check.R`).** Switching the
rule in the benchmark's simulator and keeping the seed gives an unrelated
trajectory.

- On held-out L13 the total variation between the rules is 1.0 at days 1440,
  1620, 1800 and 2520, and the largest clone differs.
- All clones share one random stream, and the first clone removed shifts every
  later draw at passaging.
- The comparison would then be between two random trajectories.

**The paired simulator (`abm_paired.cpp`).** It isolates the rule:

- every clone draws from its own random stream for each step and purpose
  (division, passaging), keyed by seed, step and karyotype;
- the population is an ordered map.

Removing a clone then leaves every other clone's draws unchanged. Two runs
differ only through the rule and what follows from it.

**The experiment (`paired_experiment.R`, `paired_summary.R`).**

- Every held-out landscape × 2 rules: 500 paired trajectories, with the
  seeds of the trajectory the benchmark scores.
- For each trajectory, the 8 core designs and the 2 "turned" stress designs
  were drawn with the benchmark's own dataset seeds (5,000 datasets per rule).
- 566 s on 16 cores (HPC).

| λ | Cells removed by the paper's rule over 2,520 days | Populations, paper's rule vs benchmark's | Datasets with identical counts |
|---|---|---|---|
| 1.2, 1.6 | none | identical | all (1,000 of 1,000 per λ) |
| 0.8 | none in 99 of 100 runs; 2 cells in one (L131) | identical | all (1,000 of 1,000) |
| 0.4 | median 3; at least 1 in 86 of 100 runs; at most 48,421 (L20, all in its "turned" designs) | identical | core 800 of 800; "turned" 199 of 200 (L20, 16 cells) |
| 0.2 | median about 14,800, at most 302,132 | total variation median 6×10⁻⁶; 4 of 100 runs drifted onto another course | core 187 of 800 |

**Why counts differ at λ = 0.2.** The population under the benchmark's rule
holds a few extra clones at frequencies near 10⁻⁶. The multinomial draw takes
the categories one after another, so one extra category shifts the rest of the
draw: the same population, sampled again.

Of the 800 core datasets:

- **187 are identical.**
- **586 are such redraws.** Their true t1 population moved by a total
  variation of at most 0.04. Their counts differ by a median of 8 cells (at
  most 113), that is, 0.9% of the sampled cells (at most 8%).
- **27 come from a population that took another course by chance:**
  - all 8 designs of three trajectories that drift before day 1440 (L412,
    L425, L436);
  - the long design of L451, which drifts by day 2520;
  - the bottleneck designs of L476 and L492, whose regrowth from 200 cells
    takes another course.

## Conclusions

1. **The effect can be measured, and it is small.** Under both rules, cells
   with f ≤ 0 leave no descendants. They hardly occur at λ ≥ 0.4, and make up
   about 10⁻⁶ to 4×10⁻⁴ of the population at λ = 0.2.
2. **At λ ≥ 0.4 the paper's rule gives the same core datasets, count for
   count** (3,200 of 3,200). The fits and scores of Part 3A would therefore be
   identical. One "turned" stress dataset of 800 differs, by 16 cells.
3. **At λ = 0.2 the data change only as much as another random seed would
   change them:** mostly a redraw of the same population, and in 27 of 800
   core datasets a trajectory that drifts away by chance. Neither is
   systematic.
4. **What does change systematically is the truth given to lethal karyotypes
   when scoring, at λ = 0.2.** On the paper's signed scale, KARST's error ratio
   rises from 0.39 to 0.41 (one step) and from 0.32 to 0.37 (two steps).
   KARST stays well below the null predictor and alfak well above it.
5. **Switching the benchmark's simulator in practice would still change its
   held-out scores.** Its shared random stream would turn every λ = 0.2
   landscape, and most λ = 0.4 ones, into new trajectories. The scores would
   move by as much as reseeding moves them, not because of the rule. The
   benchmark therefore keeps its rule, and the report states the difference
   (`doc/report`, page "The two simulated truths").

## Caveats

- **Platform.** Linux (libstdc++) and macOS (libc++) give different
  trajectories, because the standard library's random distributions differ.
  Every comparison here is within one platform. The paired runs of round
  6ae2a24 ran on the HPC, and those of round e7e32f8 ran locally.
- **Inputs.** The rescoring uses the fits of round 6ae2a24 (KARST 1.1.3). The
  rule question concerns the data and the truth, not the fitting code.

## Reproducing

Run from the repository root. `KARST_OWN_ROUND` names the round (default
6ae2a24), and `NF_CORES` sets the number of workers (default
`SLURM_CPUS_PER_TASK`, else 3).

```sh
Rscript doc/notes/nonpositive_fitness/prevalence.R          # round_<round>/prevalence.rds
Rscript doc/notes/nonpositive_fitness/rescore.R             # round_<round>/rescore.rds (about 1 h, 16 cores, 4,000 datasets)
Rscript doc/notes/nonpositive_fitness/rescore_summary.R
Rscript doc/notes/nonpositive_fitness/paired_experiment.R   # round_<round>/paired_experiment.rds (about 10 min, 16 cores)
Rscript doc/notes/nonpositive_fitness/paired_summary.R
Rscript doc/notes/nonpositive_fitness/unpaired_check.R      # the shared-stream comparison, seconds
```

**Inputs.**

- `prevalence.R` and `rescore.R` read the round's held-out results: datasets,
  scored nodes, fits of KARST and alfak, and `dataset_index_heldout.csv`.
- Their location is `OWN_RESULTS`. The default is
  `doc/report/data/own_<round>`, copied by `doc/report/build/fetch_own.sh`.
  For round 6ae2a24 that copy is the light one, without nodes and alfak fits,
  so the two scripts ran on the HPC: `OWN_RESULTS` was the round's folder, and
  `NF_HERE` was a copy of this folder.
- The simulation scripts need only this repository. `paired_experiment.R` runs
  on the held-out landscapes of `config.R`, one trajectory each, for the
  patient the benchmark scores.
