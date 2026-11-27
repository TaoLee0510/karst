# Verification of batch 2 (2026-11-26)

**What was rerun:** a sample of every step of the batch 2 studies
(`results_rc/part4_rate_choice`, C1–C7 and C8a), from `Code/` and `data/` only.
- Each step ran on the frozen upstream results (`RC_UPSTREAM` = `results_rc`).
- A step that reads frozen outputs of its own study got them through a mirror
  in `verify/batch2`. In a mirror, folders are real and files are symbolic
  links to `results_rc`. Mirrors exist only in `verify/`, never in
  `results_rc`.
- Every output of the rerun step itself was computed fresh.
- Scripts: `batch2_run.sh` (submits the dev part, `batch2_run_dev.sh`, at its
  end), `batch2_rerun.sh` (see "First run" below) and `batch2_rerun2.sh` (C2's
  figure with the Salehi reference table; C5's evaluation as of its run,
  `batch2_beta_asof.sh`).

**How it was compared:** `batch2_compare.sh` runs `compare.R --pairs`. Every
file a rerun wrote (mirror links and logs left out) is compared with the same
path in `results_rc`, under the rules of batch 1:
- R objects must be `identical` after the run-time fields are dropped;
- text tables must be equal line by line;
- anything else must have the same bytes.

**Tables:** `verify/batch2/compare/` (`summary.tsv`, one table per pair).

**Read only held:** `results_rc` is mounted read only in the image. The three
analyses of the first run that were given their arguments in the wrong order
tried to write there and failed with "Read-only file system". Nothing in
`results_rc` was changed.

**Detail of every difference:** `diff_detail.R` walks each differing file
down to its leaves (components of R objects, cells of tables). Its output is
`verify/batch2/compare/detail.tsv`, one row per differing leaf, labelled:
- **path**: a stored file path, as-run root against freeze root (same file name);
- **float**: numbers with a relative difference below 1e-6;
- **numeric**: larger numeric differences;
- **structure**: components or columns present on one side only.

PDFs are compared without their creation and modification dates.

## Results

**Every one of the 749 files the reruns wrote was compared; 472 are identical
or equal, and all 277 others are accounted for below.** (Files of the linked
input folders, `A/` and `CE/`, are not outputs and were left out.)

| Pair | Files | Identical or equal | Differ |
|---|---|---|---|
| C1 setup, fit and rows, analyses | 6 + 135 + 9 | 5 + 130 + 7 | 1 + 5 + 2 |
| C2 setup, fit, rows, analyses | 10 + 43 + 110 + 16 | 10 + 40 + 109 + 10 | 0 + 3 + 1 + 6 |
| C3 gate, task tables | 38 + 7 | 27 + 7 | 11 + 0 |
| C4 forecast_rate | 13 | 7 | 6 |
| C6 outside_domain | 6 | 2 | 4 |
| C7 marginal_pl | 7 | 2 | 5 |
| C5 and C8a dev: fits from inputs, downstream, screen, salehi_fit, evaluations | 250 + 55 + 24 + 7 + 11 | 53 + 53 + 0 + 6 + 4 | 197 + 2 + 24 + 1 + 7 |
| C5 evaluation as of its run (`batch2_beta_asof.sh`) | 2 | 0 | 2 |

The differing files fall into these classes (`detail.tsv`):
- **stored paths:** 233 leaves. File paths kept inside records, as-run root
  against freeze root: `truth_file` 180, `source` 24, `package` 19, and
  `karst_path`, `fit_path`, `root`, `support_roots` 10. In addition, C7's
  `round_dir`, which the freeze's task table names as `data/own_benchmark_dev`.
- **run times:** the `seconds` column of the C1 engine check; the time and file
  names stored in the two v7 truth records.
- **floating point:** at most a relative 4e-10, except C4's W1 (2.4e-8, the
  transport solver's tolerance, as Kester's W1 in batch 1).
- **PDFs:** 8 differ in their dates only. Two differ in the drawing: C1's
  `fig2d_truth_cells.pdf` in one full-page clipping instruction, and C3's
  `step0.pdf` (see the C3 row).
- **newer code with extra columns:** C1 and C2 rows (below).
- **the full-folder rerun of C5's evaluation:** below.
- **not results:** the Rcpp compile cache (below).

### Identical, or identical apart from stored paths and run times

| Study | Rerun | Result |
|---|---|---|
| C1 fig2d_truth | setup (cells, tasks, engine check), one fit (KARST, 3 minobs), step 1 and 2 rows, both analyses | identical apart from: KARST's `run.Rds` naming its library folder (`package`); the engine check's `seconds`. Rows: see "Newer code" below. Analysis tables identical; figure PDFs differ in their dates only, and `fig2d_truth_cells.pdf` also in one full-page clipping instruction (no visible effect) |
| C2 fig2d_misseg | setup (cells, tasks, dispersion), one fit, rows, both analyses | as C1. `fig2d_misseg.pdf` (drawn by the beta_true analysis from both arms) needs `F2T_SALEHI_TSV`, the round's Salehi reference table `analysis/fig2d_noise/fig2d_noise.tsv`, which adds the dashed Salehi lines and a caption line; the first rerun left it unset and drew the figure without them. Rerun with it (`batch2_rerun2.sh`), the figure differs in its dates only. `RUN.md` now names the variable |
| C3 dispersion gate | step 0, the three task generators, one fit each of v7, pick and vm, the step 1 analysis | task tables identical; fits identical apart from the stored library path (`karst_path`, `package`); step 1 tables equal to a relative 1.3e-11 (`vm_fits.tsv`, `rfq`); step 0's statistics (35,350 rows) equal to 1e-14. These statistics are ratios on a discrete grid with many ties, so values at a boundary fall either way: the log prints medians of exactly 0.12875 as 0.1287 or 0.1288, and in `step0.pdf` a few points at a boxplot's 1.5 IQR fence are drawn as outliers in one and inside the whisker in the other. `cells.tsv` and `s0.tsv` identical; verdict "Step 0 PASSES", as frozen |
| C4 forecast_rate | one simulation (main and varyMisseg), ALFA-K arm on both, one Salehi input, the analysis | equal up to the transport solver's tolerance: W1 of the Salehi input 2.4e-8 relative; the analysis tables follow (same verdicts: Q1 "a correct rate does not improve W1") |
| C6 outside_domain | one simulation, step 0 engine check, step 1 analysis | equal to a relative 3.7e-10 (engine check `det_vs_abm`), the stored `fit_path`; verdict "NOT REPRODUCED: stop", as frozen |
| C7 marginal_pl | one own-benchmark dataset, one paper input each of the main ABM and varyMisseg, Salehi PL, the step 0/1 analysis | equal to a relative 2.7e-13; stored paths (`fit_path`, and `round_dir`, which the freeze's task table names as `data/own_benchmark_dev`); verdict "STOP at step 0", as frozen |
| C8a noise (N1–N3, v2–v7) | one task of every fit, simulation, screen and downstream script; every evaluation | fits and simulations identical apart from stored paths (`truth_file` 180, `source` 24, `package`, `karst_path`, `root`, `support_roots`) and the time and file names stored in the v7 truth records; evaluations equal to a relative 4.4e-14 |
| C5 real-data β (fits, horizon, decomp, disp) | one task of each; `decomp_eval.R`, `disp_eval.R` | identical apart from stored paths; evaluations equal to a relative 4.4e-14 |

### Newer code: the same values, more columns

**C1 and C2 rows (`step*/rows`, `beta_fixed/rows`):**
- The rerun's rows carry columns the frozen ones lack: `pmis`, `beta`,
  `in_list_daughter_fit`, and `wdist_fallbacks`.
- Every column the two share is identical.
- The cause: C1's rows were made on 2026-11-23 before `scripts/fig2d_truth.R`
  was extended for fig2d_misseg. The extension added the rate columns and the
  counter of the transport solver's fallback (the networkflow hang fix).
  `Code/paper` holds the later version, which is the one C2 ran.

### C5's evaluation: the frozen file predates the extension fits

**What the frozen file is:**
- `beta/evaluation.rds` was written on 2026-11-20 at 05:06, and `summary.rds`
  15 s later.
- At 07:27 the decomp study's extension (`beta_salehi_ext`, tasks
  `decomp/tasks/fit_ext.tsv`) added 310 Salehi fits at ×128 and ×256 to the
  same `beta/fits/salehi/<line>/<id>/` folders.
- `beta_eval.R` reads every record in those folders. A rerun on the final
  folders therefore reads 2,803 Salehi records, where the frozen evaluation
  read 2,493.

**Rerun on the final folders (`dev_eval`):**
- The Salehi picks and their counts change (for example, E1 by single
  trajectory, one development line: "0.5:1 1:2" against "0.5:1 1:4 2:1"), and
  so does the S2 value (largest relative change 0.11).
- The verdicts (`pass` of every criterion) are identical: E1–E3 fail and E4
  passes, as reported.

**Rerun on the records as they were at 05:06 (`batch2_beta_asof.sh`):**
- The mirror holds the 18,045 `beta/` files not newer than the frozen
  `evaluation.rds`.
- Both files then reproduce to floating point (largest relative difference
  3e-15).
- So the frozen evaluation is reproducible from the freeze; it is the
  evaluation before the extension fits.

### Not results

`noise/v3/sims/*/cpp_cache/`: Rcpp's compile cache of the ABM, which the v3
simulations write next to their outputs. Its folder name is a hash per build.

## First run (18:04) and what was corrected

The first run of `batch2_run.sh` had five problems. One was in `Code/`; the
others were in the verification scripts. The affected steps were rerun by
`batch2_rerun.sh` (19:15).

1. **`Code/paper/scripts/fig2d_truth_analyse.R` failed on the truth
   simulations (C1 steps 1 and 2).**
   - The cause: "undefined columns selected", because the landscape-accuracy
     block asks for a `pmis` column, which only the varyMisseg simulations
     have.
   - The history: the script was extended for fig2d_misseg on 2026-11-23 at
     13:19. C1's analyses had been run at 11:59 with the earlier version,
     which was not kept. The 13:19 version is the as-run copy in
     `results_rc/part4_rate_choice/fig2d_misseg/code`.
   - The fix: the block now selects only the columns present. The misseg runs
     are not affected.
2. **The C4, C6 and C7 analyses were given their arguments in the wrong
   order** (results folder first, input folder second). They found no inputs
   and could not write; see "Read only held" above.
3. **`salehi_fit.R` (C8a v4) ran in the root for fits from inputs only.**
   It reads the v4 screen output. It now runs in its own root, `dev_sfit`,
   with only `noise/v4/screen/salehi` mirrored.
4. **The dev studies' downstream and evaluation steps were submitted twice.**
   - `batch2_run.sh` on the login node was taken as stopped when its ssh
     session ended, but it kept running. `batch2_run_dev.sh` (Slurm job
     22568223) submitted the same steps again.
   - Every task in `dev_down`, `dev_screen` and `dev_eval` therefore had two
     writers.
   - These folders were deleted and rerun once. `batch2_run_dev.sh` now starts
     them empty.
5. **The evaluations ran side by side, but some read each other's outputs.**
   - `disp_eval.R` reads `decomp_eval.R`'s `evaluation.rds`, and
     `beta_summary.R` reads `beta_eval.R`'s.
   - One `disp_eval` run failed for this reason, and one `beta_summary` run
     failed because the image did not start ("unknown userid").
   - The evaluations now run one after another, in that order.
