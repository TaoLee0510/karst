#!/usr/bin/env Rscript
# D3 (dev/realdata/PLAN.md, section C–E), W3: the B default Salehi fit
# (benchmark/paper settings, exactly as dev/realdata/B/fit_task.R) with
# wgd = omega, written in alfak's format at minobs 10 for forecast_score.R.
# The base of the comparison is B's own default fit (same arguments, no wgd).
#   Rscript dev/realdata/CE/d3_salehi.R TASK_FILE TASK_ID   (columns: id, method, omega)
# Output: DEV_RESULTS/realdata/CE/d3/fits/<id>/ends/<method>/minobs_10/
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
root <- file.path(rd_results, "CE", "d3", "fits")
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(rd_results, "CE", "d3", "bench_unused"))
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
settings <- bench_settings$salehi
dir <- file.path(root, task$id, "ends", task$method, "minobs_10")
if (file.exists(file.path(dir, "run.Rds"))) quit(save = "no")
cat("karst", version, find.package("karst"), "| task", paste(unlist(task), collapse = " "), "\n")
Sys.setenv(OMP_NUM_THREADS = "1")
yi <- salehi_input(task$id, ends = TRUE, paths)
x <- as.matrix(yi$x)
storage.mode(x) <- "double"
times <- input_times(yi)
arguments <- c(list(x, dt = diff(times), beta = settings$pm), settings$karst_args, list(wgd = task$omega))
started <- proc.time()[["elapsed"]]
fit <- tryCatch(do.call(karst::fit_sparse_joint, arguments), error = function(e) e)
fit_seconds <- proc.time()[["elapsed"]] - started
if (inherits(fit, "error")) {
  write_run(dir, "failed", fit_seconds, conditionMessage(fit),
            list(condition_class = class(fit)[[1L]], package = find.package("karst")))
  quit(save = "no")
}
started <- proc.time()[["elapsed"]]
result <- tryCatch({
  alfak_out <- karst::as_alfak(fit, minobs = 10, nboot = bench_settings$nboot, n0 = settings$n0, nb = settings$nb,
                               correct_efflux = bench_settings$correct_efflux, first_time = times[1])
  karst::write_alfak_outputs(alfak_out, dir)
  NULL
}, error = function(e) conditionMessage(e))
write_run(dir, if (is.null(result)) "ok" else "failed", proc.time()[["elapsed"]] - started,
          if (is.null(result)) NA_character_ else result,
          list(fit_seconds = fit_seconds, freq_origin = "time0", package = find.package("karst"), omega = task$omega))
