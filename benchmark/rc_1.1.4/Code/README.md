# Code of the KARST 1.1.4 rc freeze

**What this is:** the code that produced `../results_rc`, changed in its paths
and in one column selection (see "Changes against the as-run code"). Every
path comes from `rc_env.sh`, which sets `RC_HOME` to the
`rc_1.1.4` folder. Nothing reads or writes outside `RC_HOME`. The run order is
in `RUN.md`.

**Where the code came from:**
- The paper benchmark comes from `karst/benchmark/paper`: the working tree of
  2026-11-26, the version the HPC staging copies ran.
- The own benchmark comes from `benchmark/sparse_two_timepoint`,
  `validation/sparse_joint/far_distance.R` and
  `doc/report/build/own_*.R`, at commit 50e5889.
- The as-run copies are kept with the results:
  - `results_rc/part1_paper/salehi_rate_grid/code`;
  - `results_rc/part1_paper/round_6ae2a24/xv_target_code`;
  - `results_rc/logs/rc_work/code`;
  - the round scripts in `results_rc/part2_own/round_6ae2a24`.

## Layout

| Folder | Contents |
|---|---|
| `rc_env.sh` | paths, Slurm settings, `rc_container_rscript` (the image with `Code`, `data`, `results_rc` read only, `RC_RESULTS` writable), `rc_monitor` |
| `env/` | the image `cancer_kfl_r44_hpc_amd64_alfakr_dev2_f93b2d2.sif` (HPC copy only; R 4.4.2, alfakR dev2 f93b2d2); `src/`: karst source archives 50e5889 (1.1.4 rc) and 6ae2a24 (1.1.3), and `karst_src` (the staged 1.1.4 work tree of A3/A4); `lib/`: the libraries built from them in the image (`karst_50e5889`, `karst_6ae2a24`, `karst_wip`) and `extra_r` (CRAN packages for Fig. 4–6); `locks/` |
| `paper/` | paper benchmark (A1–A4): `R/`, `scripts/`, `hpc/` (`env.sh`, `submit.sh`, `array_task.sh`, `collect_analyse.sh`, `bg_job.sh` + `submit_rate_grid.sh` for A2, `xv_job.sh` + `submit_step2.sh` + `submit_align.sh` for A3/A4), `tasks/` (task tables of the round, the rate grid and version B) |
| `own/` | own benchmark (B1), in the repository's layout so that its relative `source()` calls work: `benchmark/sparse_two_timepoint/`, `validation/sparse_joint/`, `doc/report/build/`, `tasks/` (the round's task tables) |
| `rc/` | the rc scripts: reproduction check (part 0), the paper-grid summary (A2), the own-benchmark multiples (B2: fit, score, direction, summary, diagnostics); `tasks/` |
| `kester/` | Kester (D): setup, fits, scores, summary, ALFA-K rate diagnostic; `R/` (helpers; `R/paper/` = the paper benchmark's fit functions) |
| `hpc/` | `rc_job.sh` with `submit_rc.sh` (arrays, `TASK_FILE ID OUT`) and `submit_rc_args.sh` (one job, free arguments): the runner of `rc/`, `kester/` and `patch/` scripts |
| `dev_repo/` | batch 2's dev studies (C5 real-data β, C8a noise) and the dispersion gate (C3): the karst_dev repository copy they ran from (dev/, the validation setup chain, the older benchmark/paper/R they sourced, benchmark/sparse_two_timepoint, the package R/ and src/), pruned to these studies; `dev/hpc/env.sh` maps every root into the freeze (`DEV_RESULTS`, `DEV_BACKUP`, `PANCAN_DATA`, `V7_UPSTREAM`); the gate's code in `validation/sparse_joint/dispersion_gate` with its three task generators |
| `patch/` | batch 2's patch studies (C4 forecast_rate, C6 outside_domain, C7 marginal_pl): `R/`, `scripts/`, `tasks/`; `rc_paths.R` holds their upstream roots (`RC_UPSTREAM`); run through `hpc/rc_job.sh` |
| `monitor/` | the add-hours monitor (copies of the two scripts from the user's home folder) |
| `freeze/` | how the freeze was assembled: copy task tables, `copy_task.sh`, `submit_copy.sh`, `make_manifest.sh`, the manifest parts and their logs; `verify/`: the checks of both batches (`batch1_*.sh`, `batch2_*.sh`, `compare.R`, `diff_detail.R`) and their reports `REPORT_batch1.md`, `REPORT_batch2.md` |

## Changes against the as-run code

**Paths only:**
- `R/config.R` and `hpc/env.sh` take their roots from `RC_HOME`.
- The paper's published outputs, read as a reference by some analyses, moved
  from `PAPER_DATA/processed` to `PAPER_REF/processed`
  (`data/reference/paper`; 8 lines).
- The own benchmark reads its inputs from `STT_DATA_DIR` (`config$data_dir`).
- The rc and Kester scripts read their inputs from `data/` and the ×1 fits
  from `OWN_ROUND_FITS`.
- `rc_check`'s task table names the frozen fits relative to `RC_HOME`.
- `own/benchmark/sparse_two_timepoint/run_task.R` creates its `fits/` folder
  (the round's submit script made it beforehand).

**Removed:**
- The `horizon` mode (Salehi horizon test, outside the freeze).
- `hpc/link_alfak.sh` (linked the alfak2-era ALFA-K fits).

**Argument forms:** `kester_summary.R` and `paper_grid_summary.R` take
`TASK_FILE ID OUT`, like the other rc scripts.

**One fix that is not a path (found by the batch 2 check):**
- `paper/scripts/fig2d_truth_analyse.R`: the landscape-accuracy block now
  selects only the columns present.
- The as-run version (2026-11-23 13:19, written for fig2d_misseg) asked for
  `pmis` also on the truth simulations and stopped there.
- C1's analyses had run at 11:59 with the earlier version, which was not kept.
- No value changes (`freeze/verify/REPORT_batch2.md`).

**rc additions:** `dev_repo/dev/hpc/args_task.sh` and `submit_args.sh` run
the dev evaluations as Slurm jobs; the as-run evaluations were run with a local
Rscript. `args_task.sh` retries an image that does not start, as
`array_task.sh` does.

**Batch 2, paths only, as in batch 1:**
- Upstream results come from `RC_UPSTREAM`.
- Paths stored inside saved records or task tables are mapped from the as-run
  roots into the freeze where they are read:
  - `rc_as_run_path()` in `dev/R/common.R`;
  - `v7_check_read()`;
  - `forecast_score.R`.
- Task tables with absolute dataset folders now name them relative to
  `RC_HOME` (marginal_pl's `own_*.tsv`).
- Inputs that the as-run code read from copies (Salehi maps, daughter angles,
  window fits) are read from the copied originals in `results_rc`; each was
  checked byte-identical.

## Not runnable from here (kept as records in `results_rc`)

- **Slurm accounting steps:** `collect_analyse.sh` steps `sacct` and
  `collect_resources`.
- **`scripts/compare_rounds.R`:** needs the alfak2-era rounds.
- **Section 2 of `scripts/align_check.R`:** compares with a dev study's
  saved landscapes; it is skipped unless `DEV_EXTEND_ALFAK` is set.
- **The Fig. 2c sub-analyses:** under
  `round_6ae2a24/analysis/fig2c_two_timepoints/*` and
  `analysis/fig2c_native_20261124`. They have their own scripts, stored with
  their outputs.
- **The sources of `env/lib/extra_r`:** no longer kept. The built library is
  here; `paper/hpc/build_extra_r.sh` needs the tarballs in `env/src/extra`.
- **`patch/forecast_rate/diag/info_check.R`:** a local diagnostic that read
  the out-of-scope outside_rank study.
- **`dev_repo/dev/realdata/noise/v7/isolation_check.R`:** it searched projects
  outside the freeze for earlier copies of the v7 simulations.
- **C8b cloud_level and C9 beta_choice:** left out of the freeze (user
  decision, 2026-11-26). They are alfak2-era studies, their as-run code reads
  about 1 TB of alfak2-era fits, and the repo copies were renamed since.
