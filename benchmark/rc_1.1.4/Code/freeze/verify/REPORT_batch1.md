# Verification of batch 1 (2026-11-26)

**What was rerun:** sampled fits and every light analysis step, from `Code/`
and `data/` only. Outputs went to `verify/batch1`. They were compared with
`results_rc` by `compare.R`:
- R objects must be `identical` after the run-time fields (seconds, peak
  memory, host) are dropped;
- text tables must be equal line by line;
- anything else must have the same bytes.

**Tables:** `verify/batch1/compare/` (`summary.tsv`, one table per pair,
`own_subset.tsv`, `rc_check.tsv`).

**Static check:** no `.R` or `.sh` file in `Code/` names a path outside
`RC_HOME`. The only exceptions are the freeze tooling's source paths in
`Code/freeze`.

## Identical to the frozen results

| Check | Rerun | Result |
|---|---|---|
| E reproduction check | all 60 refits (paper Salehi, rate grid, own) with 50e5889 against the frozen fits | 60 / 60 identical (largest difference 0) |
| B1 own round fits | KARST default and independent, ALFA-K minobs 10, a distance 3–4 sweep fit | identical |
| B2 own multiples | KARST ×2 (2 fits), ×64, ALFA-K ×8 | identical |
| B2 scoring | `own_score.R` and `own_direction.R` on 20 datasets (every multiple, method, stratum) | 1,485 + 840 rows, every column equal (largest difference 0) |
| B2 summaries | `own_grid_summary.R`, `own_grid_diag.R` from the frozen tables | 10 / 10 files identical |
| A1 ABM analysis | `analyse_abm.R` from the frozen tables | 16 / 16 identical |
| A1 Salehi analysis | `analyse_salehi.R` | 27 / 27 identical, plus `sel_eval.Rds`, which the current code also writes (added for the rate grid) |
| A1 ABM forecasts | the forecast replicates of the rerun KARST fits | identical (the ABM is seeded) |
| A3 | `xv_alfak_target.R` (KARST's CV against alfak's targets, every fit) | identical |
| A4 | one version B landscape (fit on list E, forecasts, scores) | 37 / 37 identical |
| A2 paper grid | `paper_grid_summary.R` | 6 / 6 identical |
| D Kester inputs | `setup.R`: inputs, folds, task table | identical |
| D Kester log scores | `kester_score.R` on the frozen fits | every log score and NLL identical |

## Equal up to the last digit of floating point

**KARST's landscapes in the ALFA-K format:**
- Files: `landscape.Rds` and `xval.Rds` of the paper and Kester fits.
- Size: relative differences of 1e-16.
- Cause: the rerun uses the 1.1.4 library, the round the 1.1.3 library.
  Their fits are identical (reproduction check); the format conversion differs
  in the last bit.
- Consequence: in `score.Rds`, Pearson correlations agree to 1e-16, while
  Spearman correlations move by 2e-5 to 2e-3 (relative), because karyotypes
  with tied fitness are ranked differently.

**Rate-grid comparison tables:**
- Files: `fig2d_rule.tsv` and `per_line_angles.tsv`.
- Size: the 15th to 16th significant digit.
- The PDF figures differ only in their embedded metadata.

**Kester W1 / Wasserstein distances:** at most 1.3e-8, from the transport
solver's numerical tolerance.

## Not reproducible bit for bit (random draws without a fixed seed)

The original code sets no seed for:
- the posterior draws and bootstrap of the ALFA-K-format output
  (`landscape_posterior_samples.Rds`, `bootstrap_res.Rds`), for KARST and for
  ALFA-K;
- ALFA-K's own fit, which rests on its bootstrap.

A rerun therefore reproduces these only in distribution, for example on the
sampled Salehi ALFA-K fit:
- landscape mean: 0.3% relative difference;
- landscape sd: 7%.

Downstream results that use these draws reproduce in distribution as well. The
frozen files in `results_rc` are the record. Seeding them would change the
code and the results, which a freeze should not do.

## Problems found and fixed in `Code/` during the check

These are path and folder issues only, with no effect on any result:
- **`own/.../run_task.R` must create its `fits/` folder.** The round's submit
  script had made it beforehand.
- **`hpc/rc_job.sh` lets `STT_DATA_DIR` be overridden.** This was needed for
  the subset check.
- **A plot device wrote `Rplots.pdf` into the read-only working folder.**
  `rc_container_rscript` now runs scripts in a symbolic-link mirror of `Code/`
  under `RC_RESULTS/work`.
