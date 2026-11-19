#!/usr/bin/env Rscript
# History choice study, development benchmark (two time points): one dataset,
# fitted as the benchmark's karst_sparse_joint_balance (the package default:
# balance prior, candidate histories), every candidate scored against the
# true growth rates.
#
#   Rscript validation/sparse_joint/history_choice/dev_task.R DATASET_ID OUT_DIR
#
# Run from the repository root with karst 1.1.1 or later on the library path.
args <- commandArgs(TRUE)
dataset_id <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(dataset_id, ".rds"))
if (file.exists(out_file)) quit(save = "no")
script_dir <- file.path("benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "generate.R", "methods.R", "score.R")) source(file.path(script_dir, part))
source(file.path("validation", "sparse_joint", "history_choice", "capture.R"))
cat("karst", as.character(utils::packageVersion("karst")), find.package("karst"), "\n")
stopifnot(utils::packageVersion("karst") >= "1.1.1")

config <- stt_config()
method <- config$methods$karst_sparse_joint_balance
dataset <- readRDS(stt_dataset_path(config, dataset_id))
nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(dataset_id, ".rds")))
truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
counts <- dataset$counts
storage.mode(counts) <- "integer"

hc_install()
set.seed(dataset$seed)
result <- hc_fit_with_candidates(function() {
  do.call(karst::fit_sparse_joint, c(
    list(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
         min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments))
})

score_fit <- function(fit) {
  prediction <- stats::predict(fit, nodes$karyotype)
  prediction <- data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                           fitness_sd = prediction$relative_sd, stringsAsFactors = FALSE)
  metrics <- stt_score_predictions(nodes, truth_rate, prediction, config$scoring$reference_min_count)$metrics
  metrics[metrics$stratum %in% c("d0", "d1", "d2"),
          c("stratum", "n_nodes", "rmse", "rmse_zero", "spearman", "pearson", "coverage95")]
}

if (identical(result$status, "ok")) {
  scores <- lapply(seq_len(nrow(result$candidates)), function(i) {
    fit <- result$fits[[i]]
    if (is.null(fit)) return(NULL)
    cbind(history = result$candidates$history[[i]], score_fit(fit))
  })
  result$scores <- do.call(rbind, scores)
  result$fits <- NULL
  result$fit <- NULL
}
result$dataset_id <- dataset_id
result$scenario <- dataset$scenario
result$karst_version <- as.character(utils::packageVersion("karst"))
temporary <- paste0(out_file, ".tmp", Sys.getpid())
saveRDS(result, temporary)
file.rename(temporary, out_file)
cat(dataset_id, result$status, nrow(result$candidates), "candidates",
    if (!is.null(result$identical_kept)) paste("kept identical:", result$identical_kept), round(result$seconds), "s\n")
