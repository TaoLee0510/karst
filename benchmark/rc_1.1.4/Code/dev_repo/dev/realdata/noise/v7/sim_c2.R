#!/usr/bin/env Rscript
# Step 0 of dev/realdata/noise/PLAN_v7.md, C2: for one landscape, one
# missegregation multiple and one bottleneck, the Salehi-like trajectory of
# dev/realdata/CE/c2_sim.R unchanged (c2_simulate(): passages every 15 days to
# the bottleneck, cap 2e5 between them, 720 days before t0, copy numbers 1-8;
# dev/abm/abm_c2.cpp), and its datasets at t0 + 15, 30, 45, 60 days sampled
# with stt_sample_counts() at the five depth pairs.
#   Rscript dev/realdata/noise/v7/sim_c2.R TASK_FILE ID   (columns: landscape, multiple, bottleneck)
# Output: DEV_RESULTS/realdata/noise/v7/sims/c2/grid/<setting>/truth/<landscape>.rds and
#         DEV_RESULTS/realdata/noise/v7/sims/c2/datasets/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library(karst)
Sys.setenv(OMP_NUM_THREADS = "1")
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
# the globals of c2_sim.R that c2_simulate() reads, pointed at v7
c2_root <- v7_path("sims", "c2")
c2_landscapes <- v7_landscapes[, c("landscape_id", "wavelength", "landscape_seed")]
land <- v7_landscapes[task$landscape, , drop = FALSE]
m <- as.numeric(task$multiple); B <- as.numeric(task$bottleneck)
setting <- list(setting = sprintf("x%02d_B%d", m, B), round = 7L, target = "v7", beta_sim = v7_beta_default * m,
                bottleneck = B, burn_days = 720, omega = 0, gain = 0, beta_factor = 2.5, ploidy_aware = 0L,
                max_cn = 8L, seed_offset = 1L)
c2_bench_source()
started <- proc.time()[["elapsed"]]
truth_file <- c2_truth_path(setting$setting, land$landscape_id)
if (!file.exists(truth_file)) {
  truth <- c2_simulate(setting, land$landscape_id)
  sources <- c(file.path(dev_repo, "dev", "abm", "abm_c2.cpp"), file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"),
               file.path(dev_repo, "benchmark", "sparse_two_timepoint", c("generate.R", "landscape.R")),
               file.path(dev_repo, "dev", "realdata", "noise", "v7", "sim_c2.R"))
  truth$v7 <- c(list(multiple = m, half = land$half, landscape_seed = land$landscape_seed), v7_stamp(sources))
  dev_save_rds(truth, truth_file)
}
truth <- readRDS(truth_file)
made <- 0L
for (interval in c2_intervals) for (p in seq_len(nrow(v7_depth_pairs))) {
  id <- sprintf("%s_c2_x%02d_B%04d_i%03d_d%04d-%04d", land$landscape_id, m, as.integer(B), as.integer(interval),
                as.integer(v7_depth_pairs$d0[p]), as.integer(v7_depth_pairs$d1[p]))
  path <- v7_path("sims", "c2", "datasets", paste0(id, ".rds"))
  if (file.exists(path)) next
  p0 <- truth$snapshots[["0"]]; p1 <- truth$snapshots[[as.character(interval)]]
  if (!length(p0) || !length(p1)) stop("Empty population in ", id, call. = FALSE)
  counts <- stt_sample_counts(p0, p1, v7_depth_pairs$d0[p], v7_depth_pairs$d1[p], interval, v7_seed(id))
  storage.mode(counts) <- "integer"
  dev_save_rds(list(dataset_id = id, family = "c2", design = "c2", landscape_id = land$landscape_id,
                    wavelength = land$wavelength, half = land$half, multiple = m, beta = setting$beta_sim,
                    t0 = truth$t0, interval = interval, bottleneck = B, depth_t0 = v7_depth_pairs$d0[p],
                    depth_t1 = v7_depth_pairs$d1[p], seed = v7_seed(id), max_cn = setting$max_cn,
                    population_sizes = c(sum(p0), sum(p1)), counts = counts, setting = setting$setting,
                    truth_file = truth_file), path)
  made <- made + 1L
}
cat(land$landscape_id, setting$setting, ": datasets made", made, " seconds", round(proc.time()[["elapsed"]] - started), "\n")
