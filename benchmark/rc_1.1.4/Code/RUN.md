# Running the rc freeze from `Code/` and `data/`

Every command runs on the HPC login node from the rc folder
(`RC_HOME = .../karst/benchmark/rc_1.1.4`).

**Output:**
- All output goes to `RC_RESULTS`, by default `RC_HOME/rerun`, laid out like
  `results_rc/`.
- `results_rc/` itself is never written: it is mounted read only, and
  `rc_env.sh` refuses an `RC_RESULTS` inside it.

**Checking a rerun:** compare a rerun folder with the same path in
`results_rc/`, for example with `rsync -rnci rerun/X/ results_rc/X/` or the
sha256 column of `MANIFEST.tsv`.

**Resources and inputs:**
- The memory and time of every step are those of the original runs (their
  `jobs.tsv` registries are in `results_rc`).
- Each Slurm job gets an add-hours monitor (`Code/monitor`).

```bash
cd /share/lab_crd/taoli/Project/karst/benchmark/rc_1.1.4
export RC_RESULTS=$PWD/rerun          # any folder outside results_rc
source Code/rc_env.sh
```

## Part 0: reproduction check of the release candidate

Refits 60 stored KARST fits (paper Salehi, rate grid, own benchmark) with
50e5889 and compares them with the frozen fits in `results_rc`.

```bash
bash Code/hpc/submit_rc.sh rc_check rc/rc_check.R Code/rc/tasks/rc_check.tsv $RC_RESULTS/part0_freeze/rc_check2
```

## Part 1: paper benchmark (A1–A4)

All steps run from `Code/paper`. `hpc/env.sh` sets:
- `PAPER_DATA=data/paper` and `PAPER_REF=data/reference/paper`;
- `BENCH_RESULTS=$RC_RESULTS/part1_paper/round_6ae2a24`;
- `BENCH_TASKS=Code/paper/tasks/6ae2a24`;
- the karst 1.1.4 rc library.

To use the round's own 1.1.3 library, set
`KARST_LOCK=Code/env/locks/karst_6ae2a24.lock`. Its default fits are
identical (part 0).

**A1, round 6ae2a24 (every paper test, default rate):**

```bash
cd Code/paper
# KARST fits, by memory tier (as results_rc/part1_paper/round_6ae2a24/jobs.tsv)
for t in abm_mem4G:4G:06:00:00 abm_mem8G:8G:12:00:00 abm_mem16G:16G:12:00:00 abm_mem32G:32G:12:00:00 \
         misseg_mem4G:4G:06:00:00 misseg_mem8G:8G:12:00:00 misseg_mem16G:16G:12:00:00 misseg_mem32G:32G:12:00:00 \
         salehi_mem8G:8G:12:00:00 salehi_mem16G:16G:12:00:00 salehi_mem32G:32G:12:00:00; do
  IFS=: read -r f mem h m s <<< "$t"; bash hpc/submit.sh run tasks/6ae2a24/tasks_${f}.tsv $mem $h:$m:$s
done
# ALFA-K fits: the round took them from the alfak2-era round (same alfak build,
# settings and inputs); its tables come from make_tasks.R with method alfak
SCRIPT_ARGS="scripts/make_tasks.R $RC_RESULTS/part1_paper/tasks_alfak alfak" bash hpc/submit.sh script 1 4G 00:30:00
#   then: bash hpc/submit.sh run $RC_RESULTS/part1_paper/tasks_alfak/tasks_<source>.tsv MEM TIME
bash hpc/submit.sh score_salehi 277 8G 02:00:00          # after the fits (DEPEND=...)
bash hpc/submit.sh fig4 70 8G 04:00:00
bash hpc/submit.sh fig6 7 8G 06:00:00
SCRIPT_ARGS="scripts/fig2d_noise.R --stage=chunk" bash hpc/submit.sh script 28 8G 03:00:00
sbatch --qos=xxlarge --cpus-per-task=4 --mem=64G --time=06:00:00 --exclude=$RC_EXCLUDE \
  --output=$RC_RESULTS/logs/paper/%x_%j.out --export=ALL --job-name=collect_analyse hpc/collect_analyse.sh \
  check_ksets collect_abm analyse_abm example_abm collect_salehi analyse_salehi cv_score_analysis summary_comparison
SCRIPT_ARGS="scripts/fig456_summary.R" bash hpc/submit.sh script 1 16G 02:00:00
```

**Steps of A1 kept as records only (not rerun from `Code/`):**
- `sacct` and `collect_resources` read the Slurm accounting.
- `compare_rounds.R` compares with the alfak2-era rounds, which are outside
  the freeze.
- The Fig. 2c sub-analyses under `analysis/fig2c_two_timepoints/*` and
  `analysis/fig2c_native_*` were run by their own scripts, which are stored
  with their outputs.

**A2, Salehi rate grid ×2–×64:**
- ×1 is the round.
- Its fits are A1's `fits/salehi`, copied or linked into
  `salehi_rate_grid/x01/fits` before `score` and `collect` for ×1.

```bash
for M in 2 4 8 16 32 64; do
  bash Code/paper/hpc/submit_rate_grid.sh fit $M        # then score, noise, collect with DEPEND=
done
cd Code/paper && Rscript scripts/beta_grid_compare.R --root=$RC_RESULTS/part1_paper/salehi_rate_grid \
  --dispersion=$RC_RESULTS/part1_paper/salehi_rate_grid/dispersion.tsv     # (inside the image: hpc/submit.sh script)
bash Code/hpc/submit_rc.sh paper_grid rc/paper_grid_summary.R Code/rc/tasks/one.tsv $RC_RESULTS/part1_paper/paper_grid
```

`dispersion.tsv` comes from `scripts/fig2d_truth.R --stage=dispersion` (the
Salehi rows; batch 2, study fig2d_misseg).

**A3, cross-validation targets:**
- These steps run on A1's round folder, so run A1 first.
- They use library `karst_wip`, the staged 1.1.4 work tree that A3 and A4 ran
  with.

```bash
SCRIPT_ARGS="scripts/check_ksets.R" bash Code/paper/hpc/submit.sh script 1 48G 03:00:00
SCRIPT_ARGS="scripts/xv_alfak_target.R" bash Code/paper/hpc/submit.sh script 1 48G 03:00:00
XVS=alfak_target bash Code/paper/hpc/submit_step2.sh
# own target: xv_job.sh own_chunk (OWN_CHUNKS=200, array 1-200) then own_combine, then
XVS=own_target bash Code/paper/hpc/submit_step2.sh
```

**A4, version B (`ALIGN=karst`):**

```bash
bash Code/paper/hpc/submit_align.sh tasks
DEPEND=<align array> bash Code/paper/hpc/submit_align.sh downstream
# ABM: xv_job.sh check_ksets_abm collect_abm analyse_abm (16 CPUs, 96G); then align_summary, align_report
```

## Part 2: own benchmark (B1, B2)

**B1, round 6ae2a24:**
- Fits come from the task tables in `Code/own/tasks`, through
  `Code/paper/hpc/submit.sh` modes `dev` (main table) and `far` (distance
  3–4 sweep).
- Run with `STT_OUTPUT_DIR=$RC_RESULTS/part2_own/round_6ae2a24` and
  `STT_DATA_DIR=$RC_HOME/data/own_benchmark`.
- Then:
  - scoring with `run.R --stage=score` and `far_distance.R --score`;
  - the summaries `doc/report/build/own_benchmark_summary.R` and
    `own_coverage.R`.
- All of these run from `Code/own` in the image.

```bash
cd Code/paper
export STT_OUTPUT_DIR=$RC_RESULTS/part2_own/round_6ae2a24 STT_DATA_DIR=$RC_HOME/data/own_benchmark
bash hpc/submit.sh dev ../own/tasks/karst_main.tsv 16G 12:00:00
bash hpc/submit.sh dev ../own/tasks/alfak_main.tsv 16G 12:00:00        # alfak5_main(_big): 16G / 64G
bash hpc/submit.sh far ../own/tasks/karst_far.tsv 16G 12:00:00         # alfak*_far likewise
# zero predictor: run.R --stage=fit --methods=zero; scores and summaries:
CONTAINER_PWD=$RC_CODE/own SCRIPT_ARGS="benchmark/sparse_two_timepoint/run.R --stage=score --split=heldout --workers=8" \
  bash hpc/submit.sh script 1 64G 12:00:00
```

The datasets and scored nodes are inputs (`data/own_benchmark`). The ABM
truths and the generation code (`add_landscapes.R`, `generate.R`, `abm.cpp`)
are included to show how the datasets were made; a rerun does not need them.

**B2, missegregation multiples ×2–×64 (main table):**

```bash
for M in 2 4 8 16 32 64; do
  bash Code/hpc/submit_rc.sh own_karst_x$M rc/own_task.R Code/rc/tasks/own_karst_sparse_joint_balance_x$M.tsv $RC_RESULTS/part2_own/rate_grid
  bash Code/hpc/submit_rc.sh own_alfak_x$M rc/own_task.R Code/rc/tasks/own_alfak_minobs10_x$M.tsv $RC_RESULTS/part2_own/rate_grid
done
CPUS=16 MEM=64G bash Code/hpc/submit_rc.sh own_score rc/own_score.R Code/rc/tasks/one.tsv $RC_RESULTS/part2_own/rate_grid
CPUS=16 MEM=64G bash Code/hpc/submit_rc.sh own_direction rc/own_direction.R Code/rc/tasks/one.tsv $RC_RESULTS/part2_own/rate_grid
MEM=32G bash Code/hpc/submit_rc.sh own_summary rc/own_grid_summary.R Code/rc/tasks/one.tsv $RC_RESULTS/part2_own/rate_grid
MEM=32G bash Code/hpc/submit_rc.sh own_diag rc/own_grid_diag.R Code/rc/tasks/one.tsv $RC_RESULTS/part2_own/rate_grid
```

`own_score`, `own_direction` and `own_grid_diag` read the ×1 fits from
`OWN_ROUND_FITS`:
- default: `$RC_RESULTS/part2_own/round_6ae2a24/fits`, that is B1's rerun;
- to score against the frozen round instead, set
  `OWN_ROUND_FITS=$RC_FROZEN/part2_own/round_6ae2a24/fits`.

## Part 3: Kester (D)

```bash
K=$RC_RESULTS/part3_kester
bash Code/hpc/submit_rc.sh kester_setup kester/setup.R Code/rc/tasks/one.tsv $K
bash Code/hpc/submit_rc.sh kester_fit kester/kester_task.R $K/tasks.tsv $K          # 476 tasks (DEPEND=setup)
CPUS=16 OR_CORES=16 MEM=64G bash Code/hpc/submit_rc.sh kester_score kester/kester_score.R Code/rc/tasks/one.tsv $K
bash Code/hpc/submit_rc.sh kester_summary kester/kester_summary.R Code/rc/tasks/one.tsv $K
bash Code/hpc/submit_rc.sh kester_diag kester/diag_alfak.R Code/rc/tasks/one.tsv $K
```

## Part 4: missegregation-rate studies (batch 2: C1–C7, C8a)

**Upstream results:** these studies read other studies' results through
`RC_UPSTREAM`.
- Default: the run's own output root, i.e. a full rerun in the order below.
- `export RC_UPSTREAM=$RC_FROZEN` reruns a single study on the frozen results.
- Every study used the round's 1.1.3 library (`karst_6ae2a24`), except the
  dev studies, whose libraries are named per step.
- The sources of all libraries are in `Code/env/src`.

**C1 fig2d_truth and C2 fig2d_misseg** (paper code, 1.1.3 library):

```bash
export KARST_LOCK=$RC_CODE/env/locks/karst_6ae2a24.lock
bash Code/paper/hpc/submit_fig2d_truth.sh setup            # cells, tasks, engine check
bash Code/paper/hpc/submit_fig2d_truth.sh fit pilot        # then: fit rest
bash Code/paper/hpc/submit_fig2d_truth.sh rows 1 pilot 20  # step 1 rows (the round's fits); rows 1 rest 380
bash Code/paper/hpc/submit_fig2d_truth.sh rows 2 pilot 20  # step 2 rows; rows 2 rest 380
cd Code/paper && F2T_ROOT=$RC_RESULTS/part4_rate_choice/fig2d_truth F2T_BOOT=1000 \
  SCRIPT_ARGS="scripts/fig2d_truth_analyse.R 1" bash hpc/submit.sh script 1 32G 06:00:00   # and 2
bash Code/paper/hpc/submit_fig2d_misseg.sh setup           # cells, tasks, dispersion
bash Code/paper/hpc/submit_fig2d_misseg.sh fit beta_fixed pilot   # beta_fixed / beta_true x pilot / rest
bash Code/paper/hpc/submit_fig2d_misseg.sh rows beta_fixed pilot 10
cd Code/paper && F2T_ROOT=$RC_RESULTS/part4_rate_choice/fig2d_misseg F2T_SOURCE=misseg F2T_BOOT=1000 \
  F2T_SALEHI_TSV=$RC_UPSTREAM/part1_paper/round_6ae2a24/analysis/fig2d_noise/fig2d_noise.tsv \
  SCRIPT_ARGS="scripts/fig2d_truth_analyse.R beta_fixed" bash hpc/submit.sh script 1 32G 06:00:00   # then beta_true
```

The second of these runs also draws `analysis/figures/fig2d_misseg*.pdf` from
both arms; `F2T_SALEHI_TSV` adds the Salehi reference lines (dashed), as in the
as-run figure.

**C3 dispersion gate** (`Code/dev_repo/validation/sparse_joint/dispersion_gate`; 1.1.3):

```bash
G=dev_repo/validation/sparse_joint/dispersion_gate
bash Code/hpc/submit_rc_args.sh gate_step0 $G/step0_dispersion.R "--paper-main=$RC_DATA/paper/raw/ABM \
  --paper-vm=$RC_DATA/paper/raw/ABM_varyMisseg --v7-e=$RC_DATA/dispersion_gate_step0/v7_e \
  --v7-c2=$RC_DATA/dispersion_gate_step0/v7_c2 --own=$RC_DATA/dispersion_gate_step0/own \
  --out=$RC_RESULTS/part4_rate_choice/dispersion_gate/results/step0 --cores=4"
export GATE_RESULTS=$RC_RESULTS/part4_rate_choice/dispersion_gate/results
bash Code/hpc/submit_rc_args.sh mk_v7 $G/task_generators/mk_v7_tasks.R "-"   # then mk_pick_tasks.R; mk_vm_tasks.R <fig2d_misseg tasks>
GATE_TASKS=$GATE_RESULTS/step1/tasks bash Code/$G/submit_gate.sh v7 1-12000   # pick 1-3052; vm 2|4|8|16 1-375
bash Code/hpc/submit_rc_args.sh gate_analyse $G/step1_analyse.R "--v7=... --pick=... --vm=... --vm-x1=... --vm-raw=... --cells=... --out=..."
```

Notes on step 0 and the generators:
- Step 0 read subsets of the v7 and own datasets, which are kept as they were
  in `data/dispersion_gate_step0`.
- The three task generators are the ones from the session that built the
  tables.

**C4 forecast_rate, C6 outside_domain, C7 marginal_pl** (`Code/patch`; 1.1.3;
array scripts take `TASK_FILE ID OUT`):

```bash
export KARST_LIB=$RC_CODE/env/lib/karst_6ae2a24
O=$RC_RESULTS/part4_rate_choice
bash Code/hpc/submit_rc.sh fr_sim patch/forecast_rate/scripts/sim_task.R Code/patch/forecast_rate/tasks/sims.tsv $O/forecast_rate/out
bash Code/hpc/submit_rc.sh fr_alfak patch/forecast_rate/scripts/alfak_task.R Code/patch/forecast_rate/tasks/sims.tsv $O/forecast_rate/out
bash Code/hpc/submit_rc.sh fr_salehi patch/forecast_rate/scripts/salehi_task.R Code/patch/forecast_rate/tasks/salehi.tsv $O/forecast_rate/out
bash Code/hpc/submit_rc_args.sh fr_analyse patch/forecast_rate/scripts/analyse.R "$O/forecast_rate/results/final $O/forecast_rate/out"
bash Code/hpc/submit_rc.sh od_sim patch/outside_domain/scripts/sim_task.R Code/patch/outside_domain/tasks/sim_tasks.tsv $O/outside_domain/out
bash Code/hpc/submit_rc_args.sh od_step0 patch/outside_domain/scripts/step0_engine_check.R "3"
bash Code/hpc/submit_rc_args.sh od_step1 patch/outside_domain/scripts/step1_analyse.R "$O/outside_domain/results/step1 Code/patch/outside_domain/tasks/sim_tasks.tsv $O/outside_domain/out"
bash Code/hpc/submit_rc.sh mp_own patch/marginal_pl/scripts/own_task.R Code/patch/marginal_pl/tasks/own_dev.tsv $O/marginal_pl/out
bash Code/hpc/submit_rc_args.sh mp_paper patch/marginal_pl/scripts/paper_task.R "Code/patch/marginal_pl/tasks/sim_tasks.tsv $O/marginal_pl/out"
bash Code/hpc/submit_rc_args.sh mp_salehi patch/marginal_pl/scripts/salehi_pl.R "-"
bash Code/hpc/submit_rc_args.sh mp_step01 patch/marginal_pl/scripts/step01_analyse.R "$O/marginal_pl/results/step01 Code/patch/marginal_pl/tasks/sim_tasks.tsv $O/marginal_pl/out"
```

**C5 real-data β and C8a noise studies** (`Code/dev_repo`, the karst_dev
repository as run):
- `dev/hpc/env.sh` sets `DEV_RESULTS=$RC_RESULTS/part4_rate_choice/dev_results`.
- Before the first job, link the CE grid and the A table into it:
  `source Code/dev_repo/dev/hpc/env.sh; dev_link_inputs`. Add study names
  such as `beta` to link a study that is not being rerun, read from
  `RC_UPSTREAM`.
- Libraries:
  - `karst_cd342db` (1.1.2): C5, and noise n1, v2, v3 and the v4 screen;
  - `karst_v4`: the v4 fits, Salehi fit and horizon, v5 and v6;
  - `karst_6ae2a24`: the v4 base check and v7.
- Each step is `DEV_LIB=<lib> bash Code/dev_repo/dev/hpc/submit.sh SCRIPT
  TASK_FILE MEM TIME`, in the order and with the resources of
  `results_rc/part4_rate_choice/dev_results/jobs.tsv`, for example:

```bash
DEV_LIB=karst_cd342db bash Code/dev_repo/dev/hpc/submit.sh dev/realdata/beta/beta_sim_fit.R Code/dev_repo/dev/realdata/beta/tasks/beta_sim.tsv 16G 12:00:00
```

The evaluations (`beta_eval.R`, `beta_summary.R`, `decomp_eval.R`,
`disp_eval.R`, `n23_eval.R`, `v2_eval.R`, `v3_eval.R`, `v4_eval.R`,
`v5/step12_eval.R`, `v6/calib_eval.R`, `v7/calib_eval.R`) take no arguments:

```bash
DEV_LIB=karst_cd342db bash Code/dev_repo/dev/hpc/submit_args.sh dev/realdata/beta/beta_eval.R 32G 08:00:00
```

They were run on a partial local copy of the dev results; the list of the
files that copy held is `Code/freeze/verify/dev_eval_files.txt`.

**Kept as records only:**
- `patch/forecast_rate/diag/info_check.R`, which read the out-of-scope
  outside_rank study;
- `dev/realdata/noise/v7/isolation_check.R`, which searched projects outside
  the freeze for earlier copies of the v7 simulations.
