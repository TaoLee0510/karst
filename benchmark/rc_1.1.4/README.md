# KARST 1.1.4 rc: frozen benchmark results, code and inputs

**The release candidate:**
- karst 50e5889, tag `v1.1.4-rc1`.
- Default fits identical to 1.1.3 (6ae2a24); see `results_rc/part0_freeze`.
- The freeze was assembled on 2026-11-26 in two batches:
  - batch 1: the paper benchmark, KARST's own benchmark and Kester;
  - batch 2 (`results_rc/part4_rate_choice`): the studies of the
    missegregation rate: on simulations with a known rate (fig2d_truth,
    fig2d_misseg, the dispersion gate, the forecast rate, out-of-domain
    fitness, marginal PL), on real data (real-data β), and the attempts to
    estimate β from the data (noise N1–N3, v2–v7).
- Left out by decision (2026-11-26):
  - the alfak2-era cloud-level and `beta_choice` studies (their as-run code
    reads about 1 TB of alfak2-era fits);
  - dev/improve t2c.

| Folder | What it holds |
|---|---|
| `Code/` | standalone code that produces `results_rc` from `data/`; every path comes from `Code/rc_env.sh` (`RC_HOME` = this folder). See `Code/README.md` and `Code/RUN.md`. |
| `data/` | the inputs, and only those: `paper/` (ABM and varyMisseg simulations, Salehi count matrices, lineages, metadata), `reference/` (the ALFA-K paper's published outputs, read as a reference by some analyses; HPC only), `own_benchmark/` (datasets, scored nodes, index, far nodes), `kester/` (karyotype count tables); for batch 2: `own_benchmark_dev/` (the own benchmark's 168 development datasets), `dispersion_gate_step0/` (the gate's step 0 inputs), `pancan/` (the ALFA-K inputs and karyotype counts of the GBM, HNSC, ovarian and melanoma cohorts, read by the real-data studies), `dev/` (`backup/`: the alfak2-era simulated inputs e1–e4; `realdata/`: `A/salehi.rds`, the table of Salehi trajectories, and `CE/c2/grid`, the c2 simulation grid). No sequencing data. |
| `results_rc/` | every result, as real copies (no symbolic links); read only |
| `MANIFEST.tsv` | path, bytes and sha256 of every file under `Code/`, `data/` and `results_rc/`, and whether it is in the local copy |
| `PLAN.md`, `REPORT.md`, `*.R`, `hpc/`, `results/` | the working files of 2026-11-26 (plan, report and the scripts as run). They are superseded by `Code/` and `results_rc/`, which hold copies of them. |

## Where the copies are

**HPC (complete: 2.46 TB, 1.53 million files; batch 2's `part4_rate_choice` is
0.51 TB and 424,000 files of it):**
`/share/lab_crd/taoli/Project/karst/benchmark/rc_1.1.4`

**Local copy (option C; 10.6 GB, 197,191 files, every sha256 checked against `MANIFEST.tsv`):**
- Holds `Code/` (without the image and the built libraries), `data/`
  (without `reference/`), and `results_rc/` without the fit folders.
- That is, the scores, tables, analyses, figures, summaries and reports are
  local; the fits stay on the HPC.
- `MANIFEST.tsv` lists every file of the HPC copy, with a column saying
  whether it is local.

## `results_rc/`

| Folder | Item | Contents | Source (read only) |
|---|---|---|---|
| `part0_freeze/` | E | reproduction check `rc_check2` (60 refits, all identical), its superseded first version, the own-benchmark pilot | `karst_rc_1.1.4/results` |
| `part1_paper/round_6ae2a24/` | A1 | the paper benchmark at the default rate: fits (ABM, varyMisseg, Salehi; the ALFA-K fits of the alfak2-era round, same alfak build), scores, tables, analyses (including the Fig. 2c sub-analyses), figures, job registry | `benchmark/paper/results_6ae2a24` |
| ″ `analysis_xv_*`, `figures_xv_*`, `xv_target_code/` | A3 | the tests with KARST's cross-validation against alfak's targets and against its own counts | ″ |
| ″ `align_karst/` | A4 | version B: both methods on the same extended list E | ″ |
| `part1_paper/salehi_rate_grid/` | A2 | Salehi refits at 5e-5 × 2–64 (`x02`–`x64`; `x01` = the round's scores and analyses, its fits are the round's `fits/salehi`), the comparison (`analysis/`), `dispersion.tsv`, tasks, logs, the as-run code | `results_6ae2a24/beta_grid` |
| `part1_paper/paper_grid/` | A2 | the Fig. 2c / 2d tables by multiple (rc report, part 1) | local `rc_1.1.4/results` |
| `part1_paper/logs/`, `tasks/` | A1 | the paper benchmark's job logs and task tables | `benchmark/paper` |
| `part2_own/round_6ae2a24/` | B1 | the own benchmark's held-out round (500 landscapes): every method's fits, main and far table, scores, summaries, truths, the f ≤ 0 check, round scripts | `benchmark/results/sparse_two_timepoint_6ae2a24` |
| `part2_own/rate_grid/` | B2 | KARST and ALFA-K refits at 5e-5 × 2–64 (`x02`–`x64`), score tables, summaries | `karst_rc_1.1.4/results/own_grid` |
| `part3_kester/` | D | Kester: inputs, folds, fits, ABM forecasts, scores, compositions, summary | `Kester/results/forecast_50e5889` |
| `reports/` | | the reports and study documents (paper, own, rc, Kester) | the repository's docs |
| `logs/rc_work/` | | the rc work folder's code, task tables, job registry and logs | `karst_rc_1.1.4` |
| `part4_rate_choice/fig2d_truth/` | C1 | Fig. 2d on the paper's main ABM with the true landscape: cells, tasks, engine check, fits, rows, analyses, as-run code | `results_6ae2a24/fig2d_truth` |
| `part4_rate_choice/fig2d_misseg/` | C2 | Fig. 2d on varyMisseg, fitted at 5e-5 and at the true rate: as C1, plus dispersion and the one-step audit | `results_6ae2a24/fig2d_misseg` |
| `part4_rate_choice/dispersion_gate/` | C3 | the dispersion gate: step 0, the step 1 fits (v7, pick, vm) and analysis, the staged code (`results/stage1`) | `validation/sparse_joint/dispersion_gate/results` |
| `part4_rate_choice/forecast_rate/`, `outside_domain/`, `marginal_pl/` | C4, C6, C7 | the patch studies: plans, code, tasks, outputs, results, reports; `marginal_pl/upstream/` holds the 35 Salehi window fits it read from the out-of-scope outside_rank study | `karst_patch/*` |
| `part4_rate_choice/dev_results/realdata/` | C5, C8a | `beta/` (real-data β: fits, horizon, decomp, disp, evaluations) and `noise/` (N1–N3, v2–v7), in the dev code's `DEV_RESULTS` layout; `dev_results/jobs.tsv`, `logs/` | `karst_dev/results/realdata` |
| `reports/rate_choice/` | | the documents of the batch 2 studies | the repository's docs |

## How it was assembled and checked

**Copying (`Code/freeze/`):**
- `rsync -aL` in Slurm array tasks, one list of entries per task.
- Symbolic links were followed, so the ALFA-K fits that the round had linked
  from the alfak2 era are real files here.
- Each task recorded the file count and bytes of source and copy (they match
  for every task) and a sha256 per copied file.

**Results that existed only on the local machine:** added with checksums
compared. Among them are analyses run locally (CV-target report, alignment
comparison, rate-grid comparison, Kester summary tables, `report_extra`).

**One file existed in two versions:**
`round_6ae2a24/analysis/fig2c_two_timepoints/oracle_both/manifest.R`. The
local, later one is kept; the HPC one is beside it as
`manifest.R.hpc_version`.

**Batch 2 (`Code/freeze/make_copy_tasks_batch2.sh`, `_batch2b.sh`):**
- 63 copy tasks: 361,696 result files (503 GB), plus the inputs, the dev code
  and the two dev libraries. Counts and bytes match for every task.
- One task (`data_dev_backup`) ended with rsync code 23. The reason: the
  alfak2-era backup holds two dangling links, `inputs/e2/truth` and
  `inputs/e3/truth`. The truths are stored inside the datasets, which were
  copied.
- After copying, the dev results were moved into the dev code's layout:
  `realdata_beta/` became `dev_results/realdata/beta/` and
  `beta_from_data/noise/` became `dev_results/realdata/noise/`. Their sha256
  parts were rewritten to the new paths.
- The copies of the cloud-level and `beta_choice` studies were deleted once
  they were left out. Their sha256 parts are in `Code/freeze/manifest/dropped`.
- Results that existed only on the local machine were merged in, with
  checksums compared. Among them: the dev evaluations, which ran locally, and
  the patch and gate studies' locally merged outputs.
- Two sessions' scratch folders, which held files the studies read, are kept
  as `dispersion_gate/session_scratch/gate` and
  `fig2d_misseg/session_scratch/vm`.

**Verification:** `Code/freeze/verify/REPORT_batch1.md` (the tables are in
`verify/batch1/compare/`). Sampled fits and every light analysis step were
rerun from `Code/` and `data/` and compared with `results_rc`.
- **Identical:** the reproduction check (60 / 60), the own benchmark (fits,
  scores, summaries), the A1 analyses, A3, A4, the paper grid, the Kester
  inputs and log scores.
- **Equal to the last floating-point digit:** KARST's landscapes in ALFA-K
  format (1e-16) and the rate-grid comparison tables.
- **Equal only in distribution:** the unseeded random draws (posterior samples
  and bootstrap of the ALFA-K format, ALFA-K's own fits). The frozen files are
  the record.

**Verification of batch 2:** `Code/freeze/verify/REPORT_batch2.md` (tables in
`verify/batch2/compare/`, `detail.tsv` per differing value). A sample of every
step of C1–C7 and C8a was rerun on the frozen upstream results. Of the 749
files the reruns wrote:
- **472 are identical or equal.**
- **The other 277 differ** only by:
  - stored file paths;
  - run times;
  - floating point (at most 4e-10, and C4's W1 2.4e-8);
  - PDF dates;
  - extra columns of the newer fig2d code, with every shared column identical;
  - boxplot points at a tie in C3's step 0 figure.
- **C5's `beta/evaluation.rds`** predates 310 later ×128/×256 fits in the same
  folders. On the records as of its run it reproduces to 3e-15. On the final
  folders its Salehi picks change but no verdict does.
- **One code fix was needed:** a column selection in
  `fig2d_truth_analyse.R`, with no value change.
