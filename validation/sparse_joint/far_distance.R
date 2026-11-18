# Extrapolation to distances 3 and 4 from the observed karyotypes.
#
#   Rscript validation/sparse_joint/far_distance.R [--split=dev|heldout|both] [--workers=2] [--ids=a,b]
#                                                  [--out=DIR] [--methods=a,b] [--fit-only=yes]
#
# For every core dataset of the benchmark (patient 1), karyotypes three and
# four steps from the observed ones are sampled (stt_far_nodes() in
# benchmark/sparse_two_timepoint/score.R), every method is fitted again with
# these as targets, and predictions are scored against the true landscape,
# centred on the same reference karyotypes as the scoreboard. The intermediate
# nodes of the sampling walks are passed to the methods too, but not scored.
# Results go to benchmark/results/sparse_two_timepoint/sweeps/far/; the
# scoreboard's own nodes, fits and metrics are not touched. Sparse joint fits
# are saved whole (objects/), so that further distances need no refit.
# Set STT_R_LIB to the library holding the karst build under test.
# --methods fits only those methods; --fit-only=yes saves the fits and writes
# no metrics (one fit per HPC task: --ids=<one dataset> --methods=<one>
# --fit-only=yes; a final run without these options scores every fit, all
# loaded from fits/). Scoring with --methods replaces those methods' rows in
# metrics_far_<split>.csv and keeps the others', as run.R --stage=score does,
# and fits only what is missing for those methods.

library_dir <- Sys.getenv("STT_R_LIB", unset = "")
if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R")) {
  source(file.path("benchmark/sparse_two_timepoint", part))
}
args <- commandArgs(TRUE)
option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else default
}
splits <- switch(option("split", "both"), both = c("dev", "heldout"), option("split", "both"))
workers <- as.integer(option("workers", "2"))
only_ids <- option("ids", "")
fit_only <- identical(option("fit-only", "no"), "yes")
config <- stt_config()
out <- option("out", file.path(config$output_dir, "sweeps", "far"))
for (sub_dir in c("nodes", "fits", "objects")) dir.create(file.path(out, sub_dir), recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

method_names <- c("zero", "karst_sparse_joint", "karst_sparse_joint_balance",
                  "alfak_minobs5", "alfak_minobs10", "alfak_minobs20")
only_methods <- option("methods", "")
if (nzchar(only_methods)) {
  keep <- strsplit(only_methods, ",", fixed = TRUE)[[1L]]
  stopifnot(all(keep %in% method_names))
  method_names <- keep
}

far_nodes <- function(dataset) {
  path <- file.path(out, "nodes", paste0(dataset$dataset_id, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  nodes <- stt_far_nodes(dataset$counts, config$scoring$min_cn, config$scoring$max_cn, 2000L, dataset$seed + 3L)
  # written whole and then renamed: tasks of one dataset may run at the same time
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(nodes, temporary)
  file.rename(temporary, path)
  nodes
}

# Sparse joint models are fitted here rather than through the benchmark
# wrapper so that the fit object can be kept.
fit_sparse_joint_kept <- function(dataset, nodes, method, object_path) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  started <- proc.time()[["elapsed"]]
  fit <- do.call(karst::fit_sparse_joint, c(
    list(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
         min_cn = method$min_cn, max_cn = method$max_cn),
    method$arguments
  ))
  saveRDS(fit, object_path)
  prediction <- stats::predict(fit, nodes$karyotype)
  list(status = "ok", detail = sprintf("states=%d", fit$diagnostics$n_states),
       seconds = proc.time()[["elapsed"]] - started,
       prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                               fitness_sd = prediction$relative_sd, shell = prediction$shell,
                               stringsAsFactors = FALSE))
}

run_dataset <- function(id) {
  dataset <- readRDS(stt_dataset_path(config, id))
  nodes <- far_nodes(dataset)
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  scored <- nodes$scored
  rows <- lapply(method_names, function(name) {
    method <- config$methods[[name]]
    path <- file.path(out, "fits", sprintf("%s__%s.rds", id, name))
    result <- if (file.exists(path)) readRDS(path) else {
      r <- if (identical(method$type, "karst_sparse_joint")) {
        tryCatch(fit_sparse_joint_kept(dataset, nodes, method, file.path(out, "objects", sprintf("%s__%s.rds", id, name))),
                 error = function(e) list(status = "failed", detail = conditionMessage(e), seconds = NA_real_, prediction = NULL))
      } else {
        work <- file.path(tempdir(), paste0("far_", Sys.getpid()))
        dir.create(work, showWarnings = FALSE)
        stt_run_method(method, dataset, nodes, work)
      }
      saveRDS(r, path)
      r
    }
    if (fit_only) return(NULL)
    prediction <- if (is.null(result$prediction)) {
      data.frame(karyotype = character(), fitness = numeric(), stringsAsFactors = FALSE)
    } else result$prediction
    score <- stt_score_predictions(nodes[scored, ], truth[scored], prediction, config$scoring$reference_min_count)
    cbind(data.frame(dataset_id = id, method = name, status = result$status,
                     seconds = if (is.null(result$seconds)) NA_real_ else result$seconds,
                     stringsAsFactors = FALSE), score$metrics)
  })
  do.call(rbind, rows)
}

for (split in splits) {
  index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
  index <- index[index$scoreboard & index$role == "core", ]
  ids <- if (nzchar(only_ids)) intersect(strsplit(only_ids, ",")[[1L]], index$dataset_id) else index$dataset_id
  if (!length(ids)) next
  cat(sprintf("%s: %d datasets x %d methods, %d workers\n", split, length(ids), length(method_names), workers))
  results <- parallel::mclapply(ids, function(id) tryCatch(run_dataset(id), error = function(e) conditionMessage(e)),
                                mc.cores = workers, mc.preschedule = FALSE)
  failed <- vapply(results, is.character, logical(1))
  if (any(failed)) cat("Failed datasets:", paste(ids[failed], unlist(results[failed]), sep = ": ", collapse = "; "), "\n")
  if (fit_only) {
    cat("Fitted", paste(method_names, collapse = ","), "for", length(ids) - sum(failed), "dataset(s) of", split, "\n")
    if (any(failed)) quit(save = "no", status = 1)
    next
  }
  metrics <- merge(index, do.call(rbind, results[!failed]), by = "dataset_id")
  file <- if (nzchar(only_ids)) sprintf("metrics_far_%s_subset.csv", split) else sprintf("metrics_far_%s.csv", split)
  # Scoring some methods replaces their rows and keeps those of the others.
  if (nzchar(only_methods) && file.exists(file.path(out, file))) {
    previous <- utils::read.csv(file.path(out, file), stringsAsFactors = FALSE)
    previous <- previous[!previous$method %in% method_names, , drop = FALSE]
    if (nrow(previous) && !setequal(names(previous), names(metrics))) stop("columns of ", file, " differ; not merged")
    metrics <- rbind(previous[names(metrics)], metrics)
  }
  utils::write.csv(metrics, file.path(out, file), row.names = FALSE)
  cat("Wrote", file.path(out, file), nrow(metrics), "rows\n")
}
