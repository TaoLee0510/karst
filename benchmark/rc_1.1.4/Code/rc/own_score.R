# Part 2 of PLAN.md: score the own benchmark's held-out main table at every
# multiple, as run.R --stage=score does (stt_score_predictions() on the round's
# nodes and truth), for KARST default and ALFA-K minobs 10, and the round's
# zero prediction. x1 rows are the round's own fits.
#   Rscript own_score.R TASK_FILE ID OUT_DIR      (array of one; OR_CORES workers)
args <- commandArgs(TRUE)
# rc freeze: code Code/own, inputs STT_DATA_DIR, x1 fits OWN_ROUND_FITS (Code/hpc/rc_job.sh)
script_dir <- file.path(Sys.getenv("RC_HOME"), "Code", "own", "benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "generate.R", "methods.R", "score.R")) source(file.path(script_dir, part))
data_dir <- Sys.getenv("STT_DATA_DIR")
config <- stt_config(output_dir = args[[3]], data_dir = data_dir)
out <- args[[3]]
index <- utils::read.csv(file.path(data_dir, "dataset_index_heldout.csv"), stringsAsFactors = FALSE)
index <- index[index$scoreboard, , drop = FALSE]
methods <- c("karst_sparse_joint_balance", "alfak_minobs10")
cores <- as.integer(Sys.getenv("OR_CORES", "1"))
fit_dir <- function(m) if (m == 1) Sys.getenv("OWN_ROUND_FITS") else file.path(out, sprintf("x%02d", m), "fits")
rows <- parallel::mclapply(index$dataset_id, function(id) {
  dataset <- readRDS(stt_dataset_path(config, id))
  nodes <- readRDS(file.path(data_dir, "nodes", paste0(id, ".rds")))  # the round's cached nodes (run.R's stt_dataset_nodes)
  truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  r <- list()
  for (m in c(1, 2, 4, 8, 16, 32, 64)) for (meth in c(methods, if (m == 1) "zero")) {
    path <- file.path(fit_dir(m), sprintf("%s__%s.rds", id, meth))
    if (!file.exists(path)) next
    fit <- readRDS(path)
    pred <- if (is.null(fit$prediction)) data.frame(karyotype = character(), fitness = numeric()) else fit$prediction
    sc <- stt_score_predictions(nodes, truth_rate, pred, config$scoring$reference_min_count)
    r[[length(r) + 1L]] <- cbind(data.frame(dataset_id = id, multiple = m, method = meth, status = fit$status,
                                            seconds = fit$seconds, stringsAsFactors = FALSE), sc$metrics)
  }
  do.call(rbind, r)
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(rows, inherits, logical(1), "try-error")
if (any(bad)) { cat("errors:", sum(bad), "\n"); print(rows[which(bad)[1]]) }
metrics <- merge(index, do.call(rbind, rows[!bad]), by = "dataset_id")
utils::write.csv(metrics, file.path(out, "metrics_heldout_grid.csv"), row.names = FALSE)
cat("rows", nrow(metrics), "\n"); print(table(metrics$method, metrics$multiple))
