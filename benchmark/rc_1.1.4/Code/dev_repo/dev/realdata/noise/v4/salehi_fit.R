#!/usr/bin/env Rscript
# Step 4 of dev/realdata/noise/PLAN_v4.md: one fit (T0, T1, T2 or T2r) of one Salehi
# trajectory of the noise class (first and last passage, the paper
# benchmark's settings, library karst_v4), written in alfak's format at
# minobs 10 for forecast_score.R, as beta/beta_salehi.R does.
#   Rscript dev/realdata/noise/v4/salehi_fit.R TASK_FILE ID   (columns: id, method)
# Output: DEV_RESULTS/realdata/noise/v4/salehi/fits/<id>/ends/<method>/minobs_10/
#         and DEV_RESULTS/realdata/noise/v4/salehi/records/<method>/<id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
v4 <- file.path(dev_results, "realdata", "noise", "v4")
alfak_dir <- file.path(v4, "salehi", "fits", task$id, "ends", task$method, "minobs_10")
rec_file <- file.path(v4, "salehi", "records", task$method, paste0(task$id, ".rds"))
if (file.exists(rec_file) && file.exists(file.path(alfak_dir, "run.Rds"))) quit(save = "no")
scr <- readRDS(file.path(v4, "screen", "salehi", paste0(task$id, ".rds")))
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(v4, "salehi", "bench_unused"))
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
settings <- bench_settings$salehi
Sys.setenv(OMP_NUM_THREADS = "1")
yi <- salehi_input(task$id, ends = TRUE, paths)
x <- as.matrix(yi$x)
storage.mode(x) <- "double"
times <- input_times(yi)
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
beta <- switch(task$method, T0 = 5e-5, T1 = scr$beta_hat, T2 = 5e-5, T2r = 5e-5)
arguments <- c(list(x, dt = diff(times), beta = beta), settings$karst_args)
if (task$method == "T2") arguments$noise_phi <- scr$phi_hat
phi_r <- NA_real_
if (task$method == "T2r") {
  # robust variance inflation at beta_hat (PLAN_v4.md correction 16:24), as the screen: rho 1 per day
  source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
  source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
  colv <- function(j) stats::setNames(as.numeric(x[, j]), rownames(x))
  phi_r <- env_phi_robust(colv(1L), colv(ncol(x)), diff(times[c(1L, length(times))]), 1, scr$beta_hat)$phi_r
  arguments$noise_phi <- phi_r
}
set.seed(seed_of(task$id))
started <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
fit_seconds <- proc.time()[["elapsed"]] - started
key <- list(id = task$id, method = task$method, line = scr$line, group = scr$group, class = scr$class, beta = beta,
            phi_hat = scr$phi_hat, phi_r = phi_r, karst = version, fit_seconds = fit_seconds)
if (inherits(fit, "error")) {
  write_run(alfak_dir, "failed", fit_seconds, conditionMessage(fit), list(condition_class = class(fit)[[1L]]))
  dev_save_rds(c(key, list(status = "fit failed", message = conditionMessage(fit))), rec_file)
  quit(save = "no")
}
started <- proc.time()[["elapsed"]]
result <- tryCatch({
  alfak_out <- karst::as_alfak(fit, minobs = 10, nboot = bench_settings$nboot, n0 = settings$n0, nb = settings$nb,
                               correct_efflux = bench_settings$correct_efflux, first_time = times[1])
  karst::write_alfak_outputs(alfak_out, alfak_dir)
  NULL
}, error = function(e) conditionMessage(e))
write_run(alfak_dir, if (is.null(result)) "ok" else "failed", proc.time()[["elapsed"]] - started,
          if (is.null(result)) NA_character_ else result,
          list(fit_seconds = fit_seconds, freq_origin = "time0", package = find.package("karst"), beta = beta,
               noise_phi = switch(task$method, T2 = scr$phi_hat, T2r = phi_r, NA_real_)))
dev_save_rds(c(key, list(status = if (is.null(result)) "ok" else "as_alfak failed", message = result,
                         domain_coverage = fit$diagnostics$domain_coverage)), rec_file)
