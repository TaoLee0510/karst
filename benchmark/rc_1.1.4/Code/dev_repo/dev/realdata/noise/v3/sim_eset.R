#!/usr/bin/env Rscript
# New e-set simulations of dev/realdata/noise/PLAN_v3.md: the generators of
# validation/sparse_joint/multi_timepoint (e1) and series_methods (e2, e3, e4),
# unchanged, on new landscapes nz_L01-nz_L12 (seed base 111000; wavelength
# classes as mt_config()), with 200, 300, 500 and 1000 cells per time point.
# One landscape per task, each in its own output folder.
#   Rscript dev/realdata/noise/v3/sim_eset.R TASK_FILE ID
# Task columns: set (e1 | e2 | e3 | e4), landscape (1-12), stage (truth | datasets).
# truth: e1 (shared by e2, e3) and e4 (own rate). datasets: after the truth.
# Output: DEV_RESULTS/realdata/noise/v3/sims/<set>/<landscape_id>/{truth,datasets,nodes}
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
dev_source_validation()
Sys.setenv(OMP_NUM_THREADS = "1")
v3_root <- file.path(dev_results, "realdata", "noise", "v3", "sims")
base <- mt_config(tempfile(), "series")
new_landscapes <- base$landscapes
new_landscapes$landscape_id <- sprintf("nz_L%02d", seq_len(nrow(new_landscapes)))
new_landscapes$landscape_seed <- 111000L + seq_len(nrow(new_landscapes))
land <- new_landscapes[task$landscape, , drop = FALSE]
set_dir <- function(set) file.path(v3_root, set, land$landscape_id)

make_config <- function(set) {
  cfg <- if (set == "e1") mt_config(set_dir("e1"), "series") else sm_config(set)
  cfg$output_dir <- set_dir(set)
  cfg$input_dir <- set_dir(set)
  cfg$dataset_dir <- file.path(set_dir(set), "datasets")
  cfg$landscapes <- land
  cfg$depths <- c(200, 300, 500, 1000)
  if (set == "e4") cfg$beta_multiples <- if (task$landscape %% 2L == 1L) 4 else 16
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  cfg
}

started <- proc.time()[["elapsed"]]
cfg <- make_config(task$set)
if (task$stage == "truth") {
  if (task$set == "e1") si_stage_truth(cfg, 1L) else if (task$set == "e4") sm_stage_truth(cfg, 1L) else
    stop("truth only for e1 and e4")
} else if (task$stage == "datasets") {
  if (task$set == "e1") {
    mt_stage_datasets(cfg, 1L)
  } else {
    if (task$set %in% c("e2", "e3")) {
      link <- file.path(cfg$output_dir, "truth")
      if (!file.exists(link)) file.symlink(file.path(set_dir("e1"), "truth"), link)
    }
    sm_stage_datasets(cfg, 1L)
  }
} else stop("unknown stage ", task$stage)
cat(task$set, land$landscape_id, task$stage, "seconds", round(proc.time()[["elapsed"]] - started), "\n")
