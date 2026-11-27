# Coverage of the nominal 95% intervals in KARST's own benchmark (held-out, core,
# default model), with the reported sd (calibrated) and with the model's own
# posterior sd (`model_sd`, before the two model-error terms), for observed
# karyotypes (cells at either time point) and unobserved ones, by wavelength and
# distance (0-2, the scoreboard distances). Aggregated as the benchmark does:
# mean over designs, then over landscapes.
#
#   Rscript build/own_coverage.R RESULTS_DIR OUT.json    (from doc/report)
#
# OWN_CORES workers (default 2); KARST_OWN_ROUND (default: the folder's name) is
# written into the output as `round`, as build/own_benchmark_summary.R does.

args <- commandArgs(TRUE)
if (!length(args) %in% 2:3) stop("usage: Rscript build/own_coverage.R RESULTS_DIR OUT.json [ROUND]")
res <- normalizePath(args[1])
# rc freeze: the inputs (datasets, nodes, index) are in STT_DATA_DIR
dat <- normalizePath(Sys.getenv("STT_DATA_DIR", res))
cores <- as.integer(Sys.getenv("OWN_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "2")))
own_round_label <- if (length(args) >= 3) args[3] else Sys.getenv("KARST_OWN_ROUND", basename(res))
repo <- normalizePath("../..")
old_wd <- setwd(repo)
for (part in c("landscape.R", "config.R", "score.R")) source(file.path("benchmark/sparse_two_timepoint", part))
config <- stt_config()
setwd(old_wd)

method <- "karst_sparse_joint_balance"
index <- read.csv(file.path(dat, "dataset_index_heldout.csv"), stringsAsFactors = FALSE)
index <- index[index$scoreboard & index$role == "core", ]

one_set <- function(id, dataset, nodes_dir, fits_dir, keep) {
  nodes <- readRDS(file.path(nodes_dir, paste0(id, ".rds")))
  fit <- readRDS(file.path(fits_dir, sprintf("%s__%s.rds", id, method)))
  p <- fit$prediction
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  out <- lapply(c(calibrated = "fitness_sd", model = "model_sd"), function(col) {
    q <- data.frame(karyotype = p$karyotype, fitness = p$fitness, fitness_sd = p[[col]])
    stt_score_predictions(nodes, truth, q, config$scoring$reference_min_count)$nodes
  })
  t <- out$calibrated
  inside <- function(tab) abs(tab$truth_relative - tab$pred_relative) <= stats::qnorm(0.975) * tab$pred_sd
  d <- data.frame(dataset_id = id, distance = t$distance, observed = (t$count_t0 + t$count_t1) > 0,
                  calibrated = inside(out$calibrated), model = inside(out$model))
  d[d$distance %in% keep & is.finite(out$calibrated$pred_sd) & is.finite(out$model$pred_sd), ]
}

rows <- parallel::mclapply(index$dataset_id, function(id) {
  dataset <- readRDS(file.path(dat, "datasets", paste0(id, ".rds")))
  # distances 0-2 only: the fits of the far-distance sweep do not keep `model_sd`
  one_set(id, dataset, file.path(dat, "nodes"), file.path(res, "fits"), 0:2)
}, mc.cores = cores, mc.preschedule = FALSE)
bad <- vapply(rows, inherits, logical(1), "try-error")
if (any(bad)) stop("coverage failed for: ", paste(index$dataset_id[bad], collapse = ", "))
x <- merge(do.call(rbind, rows), index[c("dataset_id", "landscape_id", "wavelength_class", "scenario")])

pooled <- function(group) {
  a <- aggregate(cbind(calibrated, model) ~ ., x[c(group, "scenario", "landscape_id", "calibrated", "model")], mean)
  a <- aggregate(cbind(calibrated, model) ~ ., a[setdiff(names(a), "scenario")], mean)
  aggregate(cbind(calibrated, model) ~ ., a[setdiff(names(a), "landscape_id")], mean)
}
out <- list(results = res, round = own_round_label, n_karyotypes = nrow(x),
            by_observed = pooled("observed"), by_cell = pooled(c("wavelength_class", "distance")))
jsonlite::write_json(out, args[2], auto_unbox = TRUE, pretty = TRUE, digits = NA, dataframe = "rows")
cat("wrote", args[2], "\n")
