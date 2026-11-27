# Part 2 of PLAN.md: direction accuracy of one-copy changes on the own benchmark's
# held-out main table (distances 0-2) at every multiple, with the round report's
# definition (doc/report/build/own_benchmark_summary.R): the share of neighbouring
# pairs (one copy of one chromosome apart, true change != 0) whose predicted change
# has the true sign; ties count 1/2. x1 rows are the round's own fits.
#   Rscript own_direction.R TASK_FILE ID OUT_DIR      (array of one; OR_CORES workers)
args <- commandArgs(TRUE)
clone <- Sys.getenv("KARST_CLONE", "/share/lab_crd/taoli/Project/karst")
script_dir <- file.path(clone, "benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "score.R")) source(file.path(script_dir, part))
round_dir <- file.path(clone, "benchmark", "results", "sparse_two_timepoint_6ae2a24")
config <- stt_config(round_dir)
out <- args[[3]]
index <- utils::read.csv(file.path(round_dir, "dataset_index_heldout.csv"), stringsAsFactors = FALSE)
index <- index[index$scoreboard, , drop = FALSE]
methods <- c("karst_sparse_joint_balance", "alfak_minobs10")
cores <- as.integer(Sys.getenv("OR_CORES", "1"))
fit_dir <- function(m) if (m == 1) file.path(round_dir, "fits") else file.path(out, sprintf("x%02d", m), "fits")

direction <- function(nodes, truth, prediction) {
  table <- stt_score_predictions(nodes, truth, prediction, config$scoring$reference_min_count)$nodes
  cn <- stt_parse_karyotypes(table$karyotype)
  index_of <- stats::setNames(seq_len(nrow(table)), table$karyotype)
  edges <- do.call(rbind, lapply(seq_len(ncol(cn)), function(chr) {
    moved <- cn
    moved[, chr] <- moved[, chr] + 1
    child <- index_of[stt_format_karyotypes(moved)]
    keep <- which(!is.na(child))
    cbind(keep, unname(child[keep]))
  }))
  d <- pmax(table$distance[edges[, 1]], table$distance[edges[, 2]])
  dt <- table$truth_relative[edges[, 2]] - table$truth_relative[edges[, 1]]
  dp <- table$pred_relative[edges[, 2]] - table$pred_relative[edges[, 1]]
  dp[abs(dp) < 1e-12] <- 0
  score <- ifelse(dp == 0, 0.5, as.numeric(sign(dp) == sign(dt)))
  ok <- dt != 0
  do.call(rbind, lapply(0:2, function(k) {
    use <- ok & d == k
    data.frame(stratum = paste0("d", k), n_pairs = sum(use), direction = if (any(use)) mean(score[use]) else NA_real_)
  }))
}

rows <- parallel::mclapply(index$dataset_id, function(id) {
  dataset <- readRDS(file.path(round_dir, "datasets", paste0(id, ".rds")))
  nodes <- readRDS(file.path(round_dir, "nodes", paste0(id, ".rds")))
  truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  r <- list()
  for (m in c(1, 2, 4, 8, 16, 32, 64)) for (meth in methods) {
    path <- file.path(fit_dir(m), sprintf("%s__%s.rds", id, meth))
    if (!file.exists(path)) next
    fit <- readRDS(path)
    pred <- if (is.null(fit$prediction)) data.frame(karyotype = character(), fitness = numeric()) else fit$prediction
    r[[length(r) + 1L]] <- cbind(data.frame(dataset_id = id, multiple = m, method = meth, stringsAsFactors = FALSE),
                                 direction(nodes, truth_rate, pred))
  }
  do.call(rbind, r)
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(rows, inherits, logical(1), "try-error")
if (any(bad)) { print(rows[which(bad)[1]]); stop("errors: ", sum(bad)) }
x <- do.call(rbind, rows[!bad])
utils::write.csv(x, file.path(out, "direction_heldout_grid.csv"), row.names = FALSE)
cat("rows", nrow(x), "\n"); print(table(x$method, x$multiple))
