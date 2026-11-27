#!/usr/bin/env Rscript
# N1 of dev/realdata/noise/PLAN.md (corrections of 2026-11-20 11:40): the
# measures of R/noise.R on every input of one kind, for the first and last
# time point ("ends", the judged design) and, where there are more, for each
# pair of adjacent time points.
#   Rscript dev/realdata/noise/n1_measure.R TASK_FILE ID   (column: kind = eset | c2 | salehi | numbat)
# Output: DEV_RESULTS/realdata/noise/n1_<kind>.rds
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "R", "realdata.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
version <- rd_load_karst()
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "noise.R"))
out_file <- file.path(dev_results, "realdata", "noise", paste0("n1_", task$kind, ".rds"))
col <- function(x, j) stats::setNames(as.numeric(x[, j]), rownames(x))
harmonic <- function(v) if (length(v)) length(v) / sum(1 / v) else NA_real_
sims <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_sim.tsv"), stringsAsFactors = FALSE)
rows <- list()
add <- function(info, m) rows[[length(rows) + 1L]] <<- cbind(info, m)

if (task$kind == "eset") {
  source(file.path(dev_repo, "dev", "R", "eset.R"))
  dev_source_validation()
  for (i in which(sims$kind == "eset")) {
    loaded <- dev_eset_load(sims$set[i], sims$dataset_id[i])
    d <- loaded$dataset
    x <- mt_counts(d, "all"); tt <- d$times
    pas <- tryCatch(dev_eset_passages(sims$set[i], d), error = function(e) NULL)
    between <- function(a, b) if (is.null(pas)) data.frame(after = numeric(0)) else pas[pas$day > a + 1e-9 & pas$day <= b + 1e-9, ]
    pairs <- rbind(c(1L, ncol(x)), if (ncol(x) > 2L) cbind(seq_len(ncol(x) - 1L), seq_len(ncol(x) - 1L) + 1L))
    for (p in seq_len(nrow(pairs))) {
      a <- pairs[p, 1]; b <- pairs[p, 2]
      pb <- between(tt[a], tt[b])
      add(data.frame(kind = "eset", set = sims$set[i], dataset_id = sims$dataset_id[i], design = if (p == 1) "ends" else "adjacent",
                     from = a, to = b, true_beta = d$beta, true_bottleneck = harmonic(pb$after), replay = !is.null(pas),
                     stringsAsFactors = FALSE),
          noise_measure(col(x, a), col(x, b), tt[b] - tt[a], 1, if (is.null(pas)) NA_real_ else nrow(pb)))
    }
  }
} else if (task$kind == "c2") {
  source(file.path(dev_repo, "dev", "realdata", "CE", "c2_sim.R"))
  for (i in which(sims$kind == "c2")) {
    d <- c2_dataset(sims$set[i], sims$dataset_id[i])
    x <- d$counts
    b <- as.numeric(sub("^.*_B", "", sims$set[i]))
    add(data.frame(kind = "c2", set = sims$set[i], dataset_id = sims$dataset_id[i], design = "ends", from = 1L, to = ncol(x),
                   true_beta = d$truth$beta_sim, true_bottleneck = b, replay = TRUE, stringsAsFactors = FALSE),
        noise_measure(col(x, 1), col(x, ncol(x)), d$interval, 1, floor(d$interval / 15 + 1e-9)))
  }
} else if (task$kind == "salehi") {
  setwd(file.path(dev_repo, "benchmark", "paper"))
  Sys.setenv(BENCH_RESULTS = file.path(dev_results, "realdata", "noise", "bench_unused"))
  for (f in list.files("R", full.names = TRUE)) source(f)
  paths <- bench_paths()
  sal <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_salehi.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(sal))) {
    yi <- salehi_input(sal$id[i], ends = TRUE, paths)
    x <- as.matrix(yi$x); tt <- input_times(yi)
    passages <- diff(as.numeric(colnames(x)))
    for (rho in c(1, 0.141)) {
      add(data.frame(kind = "salehi", set = sal$line[i], dataset_id = sal$id[i], design = "ends", from = 1L, to = 2L,
                     true_beta = NA_real_, true_bottleneck = NA_real_, replay = FALSE, stringsAsFactors = FALSE),
          noise_measure(col(x, 1), col(x, 2), diff(tt), rho, passages))
    }
  }
} else if (task$kind == "numbat") {
  nb <- utils::read.delim(file.path(dev_repo, "dev", "realdata", "beta", "tasks", "beta_numbat.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(nb))) {
    r <- rd_restore(nb$cancer[i], nb$sample[i])
    if (r$status != "ok") next
    x <- rd_ends(r$restored)
    add(data.frame(kind = "numbat", set = nb$cancer[i], dataset_id = nb$sample[i], design = "ends", from = 1L, to = 2L,
                   true_beta = NA_real_, true_bottleneck = NA_real_, replay = FALSE, stringsAsFactors = FALSE),
        noise_measure(col(x, 1), col(x, 2), diff(range(r$times)), 0.141, NA_real_))
  }
} else stop("unknown kind ", task$kind)
res <- do.call(rbind, rows)
dev_save_rds(res, out_file)
cat(task$kind, nrow(res), "rows\n")
