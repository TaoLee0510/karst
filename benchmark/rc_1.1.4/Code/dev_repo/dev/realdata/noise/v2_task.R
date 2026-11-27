#!/usr/bin/env Rscript
# P1-P4 of dev/realdata/noise/PLAN_v2.md on one chunk of inputs: the envelope
# analysis of R/envelope.R on the two-time-point input, after an optional
# noise injection (variant none | dm5 | dm20 | reg05 | reg10).
#   Rscript dev/realdata/noise/v2_task.R TASK_FILE ID
# Task columns: kind (eset | c2 | salehi | numbat), set, variant, rho, chunk, chunks.
# Output: DEV_RESULTS/realdata/noise/v2/<kind>_<set>_<variant>_rho<rho>/<dataset_id>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope.R"))
out_dir <- file.path(dev_results, "realdata", "noise", "v2", sprintf("%s_%s_%s_rho%s", task$kind, task$set, task$variant, task$rho))
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
seed_of <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L
sims <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_sim.tsv"), stringsAsFactors = FALSE)

if (task$kind == "eset") {
  source(file.path(dev_repo, "dev", "R", "eset.R")); dev_source_validation()
  ids <- sims$dataset_id[sims$kind == "eset" & sims$set == task$set]
  load_one <- function(id) {
    d <- dev_eset_load(task$set, id)$dataset; x <- mt_counts(d, "ends")
    list(x0 = col(x, 1), x1 = col(x, 2), dt = diff(d$times[c(1L, length(d$times))]), truth = d$beta)
  }
} else if (task$kind == "c2") {
  source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
  ids <- sims$dataset_id[sims$kind == "c2" & sims$set == task$set]
  load_one <- function(id) {
    d <- c2_dataset(task$set, id); x <- d$counts
    list(x0 = col(x, 1), x1 = col(x, ncol(x)), dt = d$interval, truth = d$truth$beta_sim)
  }
} else if (task$kind == "salehi") {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f)
  paths <- bench_paths()
  sal <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_salehi.tsv"), stringsAsFactors = FALSE)
  ids <- sal$id
  load_one <- function(id) {
    yi <- salehi_input(id, ends = TRUE, paths); x <- as.matrix(yi$x)
    list(x0 = col(x, 1), x1 = col(x, 2), dt = diff(input_times(yi)), truth = NA_real_)
  }
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  ids <- paste(nb$cancer, nb$sample, sep = ":")
  load_one <- function(id) {
    cs <- strsplit(id, ":", fixed = TRUE)[[1]]
    r <- rd_restore(cs[1], cs[2])
    if (r$status != "ok") return(NULL)
    x <- rd_ends(r$restored)
    list(x0 = col(x, 1), x1 = col(x, 2), dt = diff(range(r$times)), truth = NA_real_)
  }
} else stop("unknown kind ", task$kind)

ids <- ids[(seq_along(ids) - 1L) %% task$chunks + 1L == task$chunk]
for (id in ids) {
  out_file <- file.path(out_dir, paste0(gsub("[:/]", "_", id), ".rds"))
  if (file.exists(out_file)) next
  inp <- load_one(id)
  if (is.null(inp)) next
  x0 <- inp$x0; x1 <- inp$x1
  if (task$variant != "none") {
    set.seed(seed_of(paste(id, task$variant)))
    kind <- if (startsWith(task$variant, "dm")) "dm" else "region"
    level <- switch(task$variant, dm5 = 5, dm20 = 20, reg05 = 0.5, reg10 = 1.0)
    fam <- if (kind == "region") env_families(x0[x0 > 0], x1[x1 > 0]) else NULL
    x0 <- env_inject(x0[x0 > 0], kind, level, fam); x1 <- env_inject(x1[x1 > 0], kind, level, fam)
  }
  started <- proc.time()[["elapsed"]]
  res <- tryCatch(env_analyse(x0, x1, inp$dt, task$rho, seed = seed_of(id)), error = function(e) e)
  rec <- if (inherits(res, "error")) list(status = "failed", message = conditionMessage(res)) else
    c(list(status = "ok"), res)
  dev_save_rds(c(list(kind = task$kind, set = task$set, dataset_id = id, variant = task$variant, rho = task$rho,
                      truth = inp$truth, seconds = proc.time()[["elapsed"]] - started), rec), out_file)
}
cat("done", nrow(task), "\n")
