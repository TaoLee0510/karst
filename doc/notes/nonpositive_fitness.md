# Clones with fitness ≤ 0: the paper's simulator and the own benchmark's

2026-11-20. Scripts and results are in `doc/notes/nonpositive_fitness/`.

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

The 210 scored held-out datasets (round e7e32f8). Medians per dataset, with
maxima in brackets:

| λ | Cells with f ≤ 0 in the population at t1 | Observed cells with f ≤ 0 | Scored nodes one step away | Two steps away |
|---|---|---|---|---|
| 0.2 | 3×10⁻⁶ (1×10⁻⁴) | 0 (0.20, a "turned" stress design) | 11% (31%) | 29% (40%) |
| 0.4 | 0 | 0 | 0 (1.9%) | 0.05% (2.8%) |
| 0.8 | 0 | 0 | 0 | 0 (0.15%) |
| 1.2, 1.6 | 0 | 0 | 0 | 0 |

The only datasets where such cells make up more than 10⁻³ of a population are
the two "turned" stress designs of landscape L14. Turning the landscape at t0 puts about 20% of
the t0 cells at f ≤ 0. Even there, under the benchmark's rule, these cells have
been thinned to nothing by t1 (360 days later).

### 2. The truth of a lethal karyotype (`rescore.R`, `rescore_summary.R`)

**The current convention.** The benchmark's truth is the per-day rate the
counts follow, log(1 + max(f, 0) Δt) / Δt. A clone with f ≤ 0 therefore gets
rate 0. That is below every positive rate, so for ranks (Spearman) these nodes
already rank last, tied.

**Two alternatives.** Under the paper's rule such a clone dies at once, which
has no finite rate. The existing fits of the 120 held-out core datasets were
rescored under two other readings:

- **death**: rate log(1 + f Δt) / Δt, that is, dying at rate |f|. This is the
  signed scale the paper scores on.
- **exclude**: these nodes are left out.

With the current convention (clip), the rescoring reproduces the benchmark's
stored scores: the largest difference is 2×10⁻¹³ in RMSE and 4×10⁻⁶ in
Spearman.

Error relative to the null predictor ("no fitness differences"), pooled over
designs and landscapes as in the benchmark, at λ = 0.2. Pooled over these
datasets, 11% of the one-step and 28% of the two-step nodes have f ≤ 0.

| | clip (current) | death | exclude |
|---|---|---|---|
| KARST, one step | 0.41 | 0.43 | 0.40 |
| KARST, two steps | 0.32 | 0.38 | 0.33 |
| alfak, minobs 5–20, one step | 35–37 | 34–36 | 38–40 |
| alfak, minobs 5–20, two steps | 30–32 | 28–29 | 34–35 |

- Observed karyotypes (distance 0) are not affected.
- At λ = 0.2, no Spearman correlation changes by more than 0.015.
- At λ ≥ 0.4, no ratio or Spearman correlation changes by more than 0.001.

### 3. The own benchmark's simulator with the paper's rule

**A plain switch cannot isolate the rule (`unpaired_check.R`).** Switching the
rule in the benchmark's simulator and keeping the seed gives an unrelated
trajectory. On held-out L13 the total variation between the rules is 1.0 at
days 1440, 1620, 1800 and 2520, and the largest clone differs. All clones share
one random stream, and the first clone removed shifts every later draw at
passaging. The comparison would then be between two random trajectories.

**The paired simulator (`abm_paired.cpp`).** It isolates the rule:

- every clone draws from its own random stream for each step and purpose
  (division, passaging), keyed by seed, step and karyotype;
- the population is an ordered map.

Removing a clone then leaves every other clone's draws unchanged. Two runs
differ only through the rule and what follows from it.

**The experiment (`paired_experiment.R`, `paired_summary.R`).**

- 15 held-out landscapes × 8 patients × 2 rules (120 paired trajectories).
- For each, the 8 core designs and the 2 "turned" stress designs were drawn
  with the benchmark's own dataset seeds (2,400 datasets).
- About 7 minutes on 3 cores.

| λ | Cells removed by the paper's rule over 2,520 days | Populations, paper's rule vs benchmark's | Datasets with identical counts |
|---|---|---|---|
| 0.8, 1.2, 1.6 | none | identical | all (240 of 240 per λ) |
| 0.4 | median 1, at most 83; ≥ 1 in 17 of 24 runs | identical | all (240 of 240) |
| 0.2 | median 9,073, at most 29,935 | total variation median 1×10⁻⁶; 1 of 24 runs (L14, patient 6) drifted onto another course | core designs 54 of 192 |

**Why counts differ at λ = 0.2.** The population under the benchmark's rule
holds a few extra clones at frequencies near 10⁻⁶. The multinomial draw takes
the categories one after another, so one extra category shifts the rest of the
draw: the same population, sampled again. In the 184 core datasets that did not
drift, the counts differ by:

- a median of 4 cells (at most 48);
- that is, 0.6% of the sampled cells (at most 5%).

## Conclusions

1. **The effect can be measured, and it is small.** Under both rules, cells
   with f ≤ 0 leave no descendants. They do not occur at λ ≥ 0.8, hardly occur
   at λ = 0.4, and make up about 10⁻⁶ to 10⁻⁴ of the population at λ = 0.2.
2. **At λ ≥ 0.4 the paper's rule gives the same datasets, count for count.**
   The fits and scores would therefore be identical.
3. **At λ = 0.2 the data change only as much as another random seed would
   change them:** a redraw of the same population, and rarely a trajectory that
   drifts away by chance. Neither is systematic.
4. **What does change systematically is the truth given to lethal karyotypes
   when scoring, at λ = 0.2.** On the paper's signed scale, KARST's error ratio
   rises from 0.41 to 0.43 (one step) and from 0.32 to 0.38 (two steps).
   KARST stays well below the null predictor and alfak far above it.
5. **Switching the benchmark's simulator in practice would still change its
   held-out scores.** Its shared random stream would turn every λ = 0.2
   landscape, and most λ = 0.4 ones, into new trajectories. The scores would
   move by as much as reseeding moves them, not because of the rule. The
   benchmark therefore keeps its rule, and the report states the difference
   (`doc/report`, page "The two simulated truths").

## Caveats

- **Platform.** Local runs (macOS, libc++) do not reproduce the HPC
  trajectories (Linux, libstdc++), because the standard library's random
  distributions differ. Every comparison here is local against local.
- **Inputs.** The rescoring uses the fits of round e7e32f8 (KARST 1.0.0). The
  rule question concerns the data and the truth, not the fitting code.

## Reproducing

Run from the repository root:

```sh
Rscript doc/notes/nonpositive_fitness/prevalence.R          # prevalence.rds
Rscript doc/notes/nonpositive_fitness/rescore.R             # rescore.rds (about 4 min, 3 cores)
Rscript doc/notes/nonpositive_fitness/rescore_summary.R
Rscript doc/notes/nonpositive_fitness/paired_experiment.R   # paired_experiment.rds (about 7 min, 3 cores)
Rscript doc/notes/nonpositive_fitness/paired_summary.R
Rscript doc/notes/nonpositive_fitness/unpaired_check.R      # the shared-stream comparison, seconds
```

**Inputs.** `prevalence.R` and `rescore.R` read the own benchmark's held-out
results: datasets, scored nodes, fits and `dataset_index_heldout.csv`. Their
location is `OWN_RESULTS`, by default `doc/report/data/own_e7e32f8`, copied by
`doc/report/build/fetch_own.sh`. The simulation scripts need only this
repository.
