#!/usr/bin/env Rscript
# New C2 simulations of dev/realdata/noise/PLAN_v3.md: dev/realdata/CE/c2_sim.R
# unchanged (its ABM and c2_simulate()), on new landscapes nz_c2_L01-L12
# (seed base 112000; wavelengths 0.8 and 1.2, six each), the six round-1
# settings, and 300, 500 and 1000 cells at both time points.
#   Rscript dev/realdata/noise/v3/sim_c2.R TASK_FILE ID   (columns: setting, landscape)
# Output: DEV_RESULTS/realdata/noise/v3/sims/c2/grid/<setting>/{truth,datasets}/
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
Sys.setenv(OMP_NUM_THREADS = "1")
# the globals of c2_sim.R that c2_simulate() reads
c2_root <- file.path(dev_results, "realdata", "noise", "v3", "sims", "c2")
c2_landscapes <- data.frame(landscape_id = sprintf("nz_c2_L%02d", 1:12), wavelength = rep(c(0.8, 1.2), each = 6L),
                            landscape_seed = 112000L + 1:12, stringsAsFactors = FALSE)
depths <- c(300, 500, 1000)
setting <- c2_setting(task$setting)
landscape_id <- c2_landscapes$landscape_id[task$landscape]
c2_bench_source()
started <- proc.time()[["elapsed"]]
truth <- c2_simulate(setting, landscape_id)
for (interval in c2_intervals) for (d in seq_along(depths)) {
  id <- sprintf("%s_i%03d_n%04d", landscape_id, as.integer(interval), as.integer(depths[d]))
  path <- file.path(c2_root, "grid", setting$setting, "datasets", paste0(id, ".rds"))
  if (file.exists(path)) next
  p0 <- truth$snapshots[["0"]]; p1 <- truth$snapshots[[as.character(interval)]]
  if (!length(p0) || !length(p1)) stop("Empty population in ", setting$setting, " ", id, call. = FALSE)
  seed <- as.integer(((c2_landscapes$landscape_seed[task$landscape] %% 1000L) * 1000L + interval +
                        setting$seed_offset * 7L) * 10L + d)
  counts <- stt_sample_counts(p0, p1, depths[d], depths[d], interval, seed)
  storage.mode(counts) <- "integer"
  dev_save_rds(list(dataset_id = id, landscape_id = landscape_id, setting = setting$setting, wavelength = truth$wavelength,
                    t0 = truth$t0, interval = interval, depth = depths[d], seed = seed, max_cn = as.integer(setting$max_cn),
                    counts = counts, truth = list(beta_sim = truth$beta_sim, bottleneck = truth$bottleneck)), path)
}
cat(task$setting, landscape_id, "seconds", round(proc.time()[["elapsed"]] - started), "\n")
