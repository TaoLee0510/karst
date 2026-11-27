#!/usr/bin/env Rscript
# Step 1 of dev/realdata/noise/PLAN_v4.md on the real data, first and last
# time point: the 155 Salehi trajectories with later passages
# (beta/decomp/tasks/fit_ext.tsv; rho 1 per day) or the NUMBAT datasets of N1
# (rho 0.141 per day, as N1; report only).
#   Rscript dev/realdata/noise/v4/screen_real.R TASK_FILE ID   (column: kind = salehi | numbat)
# Output: DEV_RESULTS/realdata/noise/v4/screen/<kind>/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
rd_load_karst()
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
out_dir <- file.path(dev_results, "realdata", "noise", "v4", "screen", task$kind)
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
screen <- function(id, x0, x1, dt, rho, meta) {
  out_file <- file.path(out_dir, paste0(id, ".rds"))
  if (file.exists(out_file)) return(invisible())
  rec <- tryCatch({
    env <- env_analyse(x0, x1, dt, rho, seed = seed_of(id))
    bb <- env_by_beta(x0, x1, dt, rho)
    cl <- env_v4_class(env, bb)
    list(status = "ok", summary = env$summary, by_beta = bb, class = cl$class, beta_hat = cl$beta_hat,
         phi_hat = cl$phi_hat)
  }, error = function(e) list(status = "failed", message = conditionMessage(e)))
  dev_save_rds(c(list(kind = task$kind, dataset_id = id, dt = dt, rho = rho, n0 = sum(x0), n1 = sum(x1)), meta, rec),
               out_file)
}
if (task$kind == "salehi") {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f)
  paths <- bench_paths()
  sal <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "decomp", "tasks", "fit_ext.tsv"),
                           stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sal))) {
    yi <- salehi_input(sal$id[i], ends = TRUE, paths)
    x <- as.matrix(yi$x); tt <- input_times(yi)
    screen(sal$id[i], col(x, 1), col(x, ncol(x)), diff(tt[c(1L, length(tt))]), 1,
           list(line = sal$line[i], group = sal$group[i]))
  }
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(nb))) {
    r <- rd_restore(nb$cancer[i], nb$sample[i])
    if (r$status != "ok") next
    x <- rd_ends(r$restored)
    screen(nb$sample[i], col(x, 1), col(x, 2), diff(range(r$times)), 0.141, list(cancer = nb$cancer[i]))
  }
} else stop("unknown kind ", task$kind)
cat(task$kind, "done\n")
