#!/usr/bin/env Rscript
# Steps 1-2 of dev/realdata/noise/PLAN_v5.md for one screen folder of v4 (the
# simulations' first and last time point, rho 1) or for the real data
# (Salehi: 155 trajectories, rho 1; NUMBAT: rho 0.141): the high-frequency
# karyotypes with their envelope test, two-proportion z and (simulations)
# true net rate, and the bounds of versions a and b.
#   Rscript dev/realdata/noise/v5/step12_task.R TASK_FILE ID   (columns: kind, tag)
# Output: DEV_RESULTS/realdata/noise/v5/anchor/<tag>.rds, a list of records.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
rd_load_karst()
Sys.setenv(OMP_NUM_THREADS = "1")
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "anchor.R"))
v4 <- file.path(dev_results, "realdata", "noise", "v4")
out_file <- file.path(dev_results, "realdata", "noise", "v5", "anchor", paste0(task$tag, ".rds"))
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
one <- function(id, x0, x1, dt, rho, meta, truth_fun = NULL) {
  rec <- tryCatch({
    r <- anc_dataset(x0, x1, dt, rho, id)
    if (!is.null(truth_fun)) {
      f <- truth_fun(r$k_all$karyotype)
      w <- r$k_all$x0 + r$k_all$x1
      r$f_mean <- sum(w * f) / sum(w)
      r$high$f_true <- f[r$k_all$high]
    }
    r$k_all <- NULL
    c(list(status = "ok"), r)
  }, error = function(e) list(status = "failed", message = conditionMessage(e)))
  c(list(dataset_id = id), meta, rec)
}

started <- proc.time()[["elapsed"]]
if (task$kind %in% c("eset", "c2")) {
  dev_source_validation()
  source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
  c2_bench_source()
  files <- list.files(file.path(v4, "screen", task$tag), pattern = "[.]rds$", full.names = TRUE)
  recs <- lapply(files, function(p) {
    scr <- readRDS(p)
    meta <- scr[c("kind", "set", "variant", "landscape_id", "half", "depth", "bottleneck", "truth", "dt")]
    if (!identical(scr$status, "ok")) return(c(list(dataset_id = scr$dataset_id, status = "no screen"), meta))
    d <- readRDS(scr$source)
    truth_fun <- if (task$kind == "eset") function(labels) sm_true_net_rate(labels, d$truth, d$beta) else {
      d_truth <- d
      d_truth$truth <- readRDS(file.path(dirname(dirname(scr$source)), "truth", paste0(d$landscape_id, ".rds")))
      function(labels) c2_true_rate(labels, d_truth)
    }
    one(scr$dataset_id, scr$x0, scr$x1, scr$dt, 1, meta, truth_fun)
  })
} else if (task$kind == "salehi") {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f)
  paths <- bench_paths()
  sal <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "decomp", "tasks", "fit_ext.tsv"),
                           stringsAsFactors = FALSE)
  recs <- lapply(seq_len(nrow(sal)), function(i) {
    yi <- salehi_input(sal$id[i], ends = TRUE, paths)
    x <- as.matrix(yi$x); tt <- input_times(yi)
    one(sal$id[i], col(x, 1), col(x, ncol(x)), diff(tt[c(1L, length(tt))]), 1,
        list(kind = "salehi", line = sal$line[i], group = sal$group[i]))
  })
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  recs <- lapply(seq_len(nrow(nb)), function(i) {
    r <- rd_restore(nb$cancer[i], nb$sample[i])
    if (r$status != "ok") return(NULL)
    x <- rd_ends(r$restored)
    one(nb$sample[i], col(x, 1), col(x, 2), diff(range(r$times)), 0.141, list(kind = "numbat", cancer = nb$cancer[i]))
  })
  recs <- Filter(Negate(is.null), recs)
} else stop("unknown kind ", task$kind)
dev_save_rds(recs, out_file)
cat(task$tag, length(recs), "records,", sum(vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "ok,",
    round(proc.time()[["elapsed"]] - started), "s\n")
