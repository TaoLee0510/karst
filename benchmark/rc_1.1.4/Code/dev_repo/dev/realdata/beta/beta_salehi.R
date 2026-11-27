#!/usr/bin/env Rscript
# beta0/beta1 of dev/realdata/beta/PLAN.md on Salehi: one trajectory (first
# and last passage, the paper benchmark's settings), fitted at every rate of
# the grid, with the three criteria of beta1; when `horizon` is 1, also
# written in alfak's format at minobs 10 for forecast_score.R (E4, S1).
#   Rscript dev/realdata/beta/beta_salehi.R TASK_FILE ID [MULTIPLES, e.g. 128:256]   (columns: id, line, group, horizon)
# Output: DEV_RESULTS/realdata/beta/fits/salehi/<line>/<id>/<b...>.rds and
#         DEV_RESULTS/realdata/beta/salehi/fits/<id>/ends/<b...>/minobs_10/
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
source(file.path(dev_repo, "dev", "realdata", "beta", "R", "beta.R"))
# decomp/PLAN.md: more rates on the same trajectories as a third argument, e.g. 128:256
# (SCRIPT_EXTRA; environment variables do not reach R inside the container, and
# sbatch --export splits values at commas, so the rates are separated by colons)
if (length(args) >= 3L) beta_grid <- as.numeric(unlist(strsplit(args[-(1:2)], "[,:]")))
out_dir <- file.path(beta_results, "fits", "salehi", task$line, task$id)
alfak_root <- file.path(beta_results, "salehi", "fits")
setwd(file.path(dev_repo, "benchmark", "paper"))
Sys.setenv(BENCH_RESULTS = file.path(beta_results, "salehi", "bench_unused"))
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
settings <- bench_settings$salehi
Sys.setenv(OMP_NUM_THREADS = "1")
yi <- salehi_input(task$id, ends = TRUE, paths)
x <- as.matrix(yi$x)
storage.mode(x) <- "double"
times <- input_times(yi)
seed <- beta_seed(task$id)
for (m in beta_grid) {
  label <- beta_label(m)
  path <- file.path(out_dir, paste0(label, ".rds"))
  alfak_dir <- file.path(alfak_root, task$id, "ends", label, "minobs_10")
  if (file.exists(path) && (task$horizon != 1 || file.exists(file.path(alfak_dir, "run.Rds")))) next
  set.seed(seed)
  run <- beta_fit_timed(c(list(x, dt = diff(times), beta = beta_default * m), settings$karst_args))
  key <- list(kind = "salehi", line = task$line, group = task$group, id = task$id, beta = beta_default * m,
              label = label, multiple = m, seconds = run$seconds, karst = version,
              cells = colSums(x), dt = diff(times))
  if (inherits(run$fit, "error")) {
    record <- c(key, beta_failure(run$fit))
    if (task$horizon == 1) write_run(alfak_dir, "failed", run$seconds, conditionMessage(run$fit),
                                     list(condition_class = class(run$fit)[[1L]], package = find.package("karst")))
  } else {
    fit <- run$fit
    record <- c(key, list(status = "ok", criteria = beta_criteria(fit, seed),
                          domain_coverage = fit$diagnostics$domain_coverage))
    if (task$horizon == 1) {
      started <- proc.time()[["elapsed"]]
      result <- tryCatch({
        alfak_out <- karst::as_alfak(fit, minobs = 10, nboot = bench_settings$nboot, n0 = settings$n0, nb = settings$nb,
                                     correct_efflux = bench_settings$correct_efflux, first_time = times[1])
        karst::write_alfak_outputs(alfak_out, alfak_dir)
        NULL
      }, error = function(e) conditionMessage(e))
      write_run(alfak_dir, if (is.null(result)) "ok" else "failed", proc.time()[["elapsed"]] - started,
                if (is.null(result)) NA_character_ else result,
                list(fit_seconds = run$seconds, freq_origin = "time0", package = find.package("karst"), beta = beta_default * m))
    }
  }
  dev_save_rds(record, path)
}
