#!/usr/bin/env Rscript
# Step 0 of dev/realdata/noise/PLAN_v7.md, e sets: for one landscape and one
# missegregation multiple, the truth population (stt_generate_truth(), grown
# from day 0 at that rate, as e4) and the windows of the three designs with
# the validation code unchanged (e1 selection: mt_make_window(); e2 clonal
# sweeps and e3 neutral: sm_make_window()); first and last time point of 4 or
# 8 sampled with stt_sample_counts() at the five depth pairs.
#   Rscript dev/realdata/noise/v7/sim_eset.R TASK_FILE ID   (columns: landscape, multiple)
# Output: DEV_RESULTS/realdata/noise/v7/sims/e/x<m>/truth/<landscape>_P01.rds and
#         DEV_RESULTS/realdata/noise/v7/sims/e/datasets/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library(karst)
Sys.setenv(OMP_NUM_THREADS = "1")
dev_source_validation()
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
land <- v7_landscapes[task$landscape, , drop = FALSE]
m <- as.numeric(task$multiple); mi <- match(m, v7_multiples) - 1L
started <- proc.time()[["elapsed"]]

# compile both simulators in a private cache first (parallel tasks)
cache_cfg <- mt_config(file.path(tempdir(), "v7_cpp"), "series")
sm_load_abm(cache_cfg)
sources <- c(file.path(dev_repo, "benchmark", "sparse_two_timepoint", c("abm.cpp", "generate.R", "landscape.R")),
             file.path(dev_repo, "validation", "sparse_joint", "series_methods", c("abm_constant.cpp", "setup.R")),
             file.path(dev_repo, "validation", "sparse_joint", "multi_timepoint", "setup.R"),
             file.path(dev_repo, "dev", "realdata", "noise", "v7", "sim_eset.R"))

configure <- function(set) {
  cfg <- if (set == "e1") mt_config(tempfile(), "series") else sm_config(set)
  cfg$output_dir <- v7_path("sims", "e", sprintf("x%02d", m))
  cfg$input_dir <- cfg$output_dir
  cfg$landscapes <- land[, c("landscape_id", "split", "wavelength_class", "wavelength", "landscape_seed")]
  cfg$beta <- v7_beta_default * m
  cfg
}

## truth
cfg1 <- configure("e1")
truth_file <- stt_truth_path(cfg1, land$landscape_id, 1L)
if (!file.exists(truth_file)) {
  stt_generate_truth(cfg1, cfg1$landscapes[1, ], 1L)
  tr <- readRDS(truth_file)
  tr$v7 <- c(list(multiple = m, half = land$half, landscape_seed = land$landscape_seed), v7_stamp(sources))
  dev_save_rds(tr, truth_file)
}
truth <- readRDS(truth_file)
if (!isTRUE(all.equal(truth$beta, v7_beta_default * m))) stop("truth rate mismatch")

## windows and datasets
made <- 0L
for (design in c("e1", "e2", "e3")) {
  cfg <- configure(design)
  plan <- if (design == "e1") mt_window_plan(cfg) else sm_window_plan(cfg)
  for (item in plan) {
    item$seed <- as.integer(item$seed + mi * 1000000L)
    start <- if (design == "e1") cfg$t0 else item$start
    ids <- vapply(cfg$samples, function(k) vapply(seq_len(nrow(v7_depth_pairs)), function(p)
      sprintf("%s_%s_x%02d_t%04d_b%05d_s%03d_k%d_d%04d-%04d", land$landscape_id, design, m, as.integer(start),
              as.integer(item$bottleneck), as.integer(item$spacing), k, as.integer(v7_depth_pairs$d0[p]),
              as.integer(v7_depth_pairs$d1[p])), ""), character(nrow(v7_depth_pairs)))
    paths <- v7_path("sims", "e", "datasets", paste0(ids, ".rds"))
    if (all(file.exists(paths))) next
    window <- if (design == "e1") list(populations = mt_make_window(cfg, truth, item), constant_fitness = 0) else
      sm_make_window(cfg, truth, item)
    pops <- window$populations
    for (ki in seq_along(cfg$samples)) for (p in seq_len(nrow(v7_depth_pairs))) {
      k <- cfg$samples[[ki]]; id <- ids[p, ki]; path <- v7_path("sims", "e", "datasets", paste0(id, ".rds"))
      if (file.exists(path)) next
      interval <- item$spacing * (k - 1L)
      counts <- stt_sample_counts(pops[[1L]], pops[[k]], v7_depth_pairs$d0[p], v7_depth_pairs$d1[p], interval, v7_seed(id))
      storage.mode(counts) <- "integer"
      dev_save_rds(list(
        dataset_id = id, family = "e", design = design, landscape_id = land$landscape_id, wavelength = land$wavelength,
        half = land$half, multiple = m, beta = truth$beta, t0 = start, interval = interval, samples = k,
        spacing = item$spacing, bottleneck = item$bottleneck, depth_t0 = v7_depth_pairs$d0[p],
        depth_t1 = v7_depth_pairs$d1[p], seed = v7_seed(id), window_seed = item$seed,
        population_sizes = c(sum(pops[[1L]]), sum(pops[[k]])), counts = counts,
        truth = list(centroids = truth$centroids, wavelength = truth$wavelength, abm_delta_t = truth$abm_delta_t,
                     constant_fitness = if (window$constant_fitness > 0) window$constant_fitness),
        truth_file = truth_file), path)
      made <- made + 1L
    }
  }
}
cat(land$landscape_id, "x", m, ": datasets made", made, " seconds", round(proc.time()[["elapsed"]] - started), "\n")
