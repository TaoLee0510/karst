# Reproduction check of the release candidate (PLAN.md, part 0): refit one
# stored KARST fit with karst 50e5889 and compare.
#   Rscript rc_check.R TASK_FILE ID OUT_DIR
# TASK_FILE columns: type (paper | betagrid | own), path (stored fit, for
# paper / betagrid), dataset_id, round_dir (own).
args <- commandArgs(TRUE)
task <- utils::read.delim(args[[1]], stringsAsFactors = FALSE)[as.integer(args[[2]]), ]
out_dir <- args[[3]]; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("check_%03d.rds", as.integer(args[[2]])))
stopifnot(normalizePath(find.package("karst")) == normalizePath(file.path(Sys.getenv("KARST_LIB"), "karst")))
res <- list(task = task, karst = as.character(utils::packageVersion("karst")))
maxdiff <- function(a, b) { a <- as.numeric(a); b <- as.numeric(b); if (length(a) != length(b)) Inf else max(abs(a - b), 0, na.rm = FALSE) }
if (task$type %in% c("paper", "betagrid")) {
  old <- readRDS(task$path)
  # Refit from the round's input, as the paper's run_task.R / fit_karst() do
  # (a multi-time-point fit stores only the time points it kept, so its own
  # data are not the input): .../fits/salehi/<trajectory>/<design>/karst/...
  parts <- strsplit(task$path, "/", fixed = TRUE)[[1]]
  i <- which(parts == "salehi")[1]; traj <- parts[i + 1L]; design <- parts[i + 2L]
  yi <- readRDS(file.path(Sys.getenv("PAPER_DATA"), "processed", "salehi", "alfak_inputs", paste0(traj, ".Rds")))
  if (design == "ends" && ncol(yi$x) > 2L) yi$x <- yi$x[, c(1L, ncol(yi$x)), drop = FALSE]
  x <- as.matrix(yi$x); storage.mode(x) <- "double"
  times <- as.numeric(colnames(x)) * yi$dt
  new <- karst::fit_sparse_joint(x, dt = diff(times), beta = old$data$beta, min_cn = 0)
  num <- vapply(old$summary, is.numeric, logical(1))
  labs <- c(old$summary$karyotype, old$dropped)
  po <- stats::predict(old, labs); pn <- stats::predict(new, labs)
  res$diff <- c(objective = abs(new$diagnostics$objective - old$diagnostics$objective),
                summary = max(vapply(names(old$summary)[num], function(n) maxdiff(old$summary[[n]], new$summary[[n]]), numeric(1))),
                prediction = maxdiff(po$relative_fitness, pn$relative_fitness), prediction_sd = maxdiff(po$relative_sd, pn$relative_sd))
  res$same_karyotypes <- identical(old$summary$karyotype, new$summary$karyotype) && identical(old$dropped, new$dropped)
} else {
  clone <- Sys.getenv("KARST_CLONE")
  for (part in c("landscape.R", "score.R")) source(file.path(clone, "benchmark", "sparse_two_timepoint", part))
  dataset <- readRDS(file.path(task$round_dir, "datasets", paste0(task$dataset_id, ".rds")))
  nodes <- readRDS(file.path(task$round_dir, "nodes", paste0(task$dataset_id, ".rds")))
  counts <- dataset$counts; storage.mode(counts) <- "integer"
  set.seed(dataset$seed)
  fit <- karst::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                                 min_cn = 1L, max_cn = 8L, initial_state = "balance")
  st <- readRDS(file.path(task$round_dir, "fits", paste0(task$dataset_id, "__karst_sparse_joint_balance.rds")))
  pr <- stats::predict(fit, st$prediction$karyotype)
  res$diff <- c(prediction = maxdiff(st$prediction$fitness, pr$relative_fitness), prediction_sd = maxdiff(st$prediction$fitness_sd, pr$relative_sd))
  res$same_karyotypes <- TRUE
}
res$identical <- all(res$diff == 0) && res$same_karyotypes
saveRDS(res, out_file)
cat(task$type, if (is.na(task$path)) "" else basename(dirname(dirname(dirname(task$path)))), "identical:", res$identical, " max diffs:", paste(names(res$diff), signif(res$diff, 3), collapse = " "), "\n")
