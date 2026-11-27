# Copy of validation/sparse_joint/calibration_fit.R for the short-interval
# study: reads the fits of results/devbench (devbench.R) instead of the
# development benchmark's own results. Fits the constants that turn the
# model's posterior standard deviation of relative fitness into the reported
# one, with the package's formula (sparse_joint_calibrate() in
# R/sparse-joint.R):
#
#   reported sd = sqrt(model sd^2 + (deficit_k * relative fitness)^2 + (floor * birth rate * (1 + shell))^2)
#
# where relative fitness is the fit's own (relative to the count-weighted mean
# of its major karyotypes, not the scoring reference) and deficit_k is
# `deficit` for observed karyotypes and `deficit_unobserved` for karyotypes
# seen at neither time point (here: distance 0 or not).
#
#   Rscript validation/sparse_joint/short_interval/calibration_fit.R <method> [check: floor,deficit[,deficit_unobserved]]
#
# Without `check`, fits floor and deficit by Gaussian likelihood on the core
# scenarios of the split, with unobserved karyotypes taking the same deficit
# (the two-constant form), and prints the coverage they give. With `check`,
# only evaluates the given constants, for instance the package's
# 0.009,0.6,1.1 or another method; two constants give the two-constant form.
# Every dataset, and within it every kind of karyotype, has the same weight, so
# the 4000 unobserved karyotypes of a dataset do not outvote its 20 abundant
# ones.

si_script_dir <- "validation/sparse_joint/short_interval"
source(file.path(si_script_dir, "setup.R"))
for (part in c("methods.R", "report.R")) source(file.path(si_repo, "benchmark", "sparse_two_timepoint", part))
args <- commandArgs(TRUE)
method <- args[[1L]]
split <- "dev"
check <- if (length(args) > 1L) as.numeric(strsplit(args[[2L]], ",", fixed = TRUE)[[1L]]) else NULL
config <- stt_config(file.path(si_default_output, "devbench"))
birth_rate <- 1
index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
index <- index[index$scoreboard, ]

nodes <- do.call(rbind, lapply(index$dataset_id, function(id) {
  path <- file.path(config$output_dir, "fits", sprintf("%s__%s.rds", id, method))
  if (!file.exists(path)) return(NULL)
  prediction <- readRDS(path)$prediction
  if (is.null(prediction)) return(NULL)
  dataset <- readRDS(stt_dataset_path(config, id))
  scored <- readRDS(file.path(config$output_dir, "nodes", paste0(id, ".rds")))
  truth <- stt_true_growth_rate(scored$karyotype, dataset$truth)
  table <- stt_score_predictions(scored, truth, prediction, config$scoring$reference_min_count)$nodes
  at <- match(table$karyotype, prediction$karyotype)
  table$model_sd <- prediction$model_sd[at]
  table$shell <- prediction$shell[at]
  # The fit's relative fitness, which the package's deficit multiplies.
  table$relative <- prediction$fitness[at]
  cbind(dataset_id = id, scenario = dataset$scenario, table)
}))
nodes <- merge(nodes, index[c("dataset_id", "role")])
nodes <- nodes[is.finite(nodes$model_sd) & is.finite(nodes$relative), ]
nodes$error <- nodes$pred_relative - nodes$truth_relative
total <- nodes$count_t0 + nodes$count_t1
nodes$kind <- ifelse(nodes$distance == 0L, ifelse(total >= 5, "observed, 5+ cells", "observed, <5 cells"),
                     paste("unobserved, distance", nodes$distance))
cell <- paste(nodes$dataset_id, nodes$kind)
nodes$weight <- 1 / as.numeric(table(cell)[cell])

# `constants`: floor, deficit and optionally deficit_unobserved (otherwise
# `deficit` for all), as in the package.
reported_sd <- function(constants, part) {
  unobserved <- if (length(constants) > 2L) constants[[3L]] else constants[[2L]]
  deficit <- ifelse(part$distance == 0L, constants[[2L]], unobserved)
  sqrt(part$model_sd^2 + (deficit * part$relative)^2 +
         (constants[[1L]] * birth_rate * (1 + part$shell))^2)
}
core <- nodes[nodes$role == "core", ]
if (is.null(check)) {
  objective <- function(log_constants) {
    variance <- reported_sd(exp(log_constants), core)^2
    sum(core$weight * (log(variance) + core$error^2 / variance))
  }
  best <- stats::optim(log(c(0.01, 0.5)), objective, control = list(maxit = 2000, reltol = 1e-10))
  constants <- exp(best$par)
  cat(sprintf("Fitted on %d datasets (%s, core scenarios): floor = %.4f  deficit = %.3f  (convergence %d)\n",
              length(unique(core$dataset_id)), split, constants[1L], constants[2L], best$convergence))
} else {
  constants <- check
  cat(sprintf("Checking floor = %.4f  deficit = %.3f  deficit_unobserved = %.3f on %s, method %s\n",
              constants[1L], constants[2L], if (length(constants) > 2L) constants[3L] else constants[2L],
              split, method))
}

coverage <- function(part, by) {
  part$covered <- abs(part$error) <= stats::qnorm(0.975) * reported_sd(constants, part)
  part$covered_model <- abs(part$error) <= stats::qnorm(0.975) * part$model_sd
  part$reported <- reported_sd(constants, part)
  per_dataset <- stats::aggregate(part[c("covered", "covered_model", "reported", "model_sd")],
                                  part[c("dataset_id", by)], mean)
  per_dataset$abs_error <- stats::aggregate(abs(part$error), part[c("dataset_id", by)], mean)$x
  out <- stats::aggregate(per_dataset[c("covered", "covered_model", "reported", "model_sd", "abs_error")],
                          per_dataset[by], mean)
  out$worst_dataset <- stats::aggregate(per_dataset$covered, per_dataset[by], min)$x
  out$datasets <- stats::aggregate(per_dataset$covered, per_dataset[by], length)$x
  # after the key columns (one or two)
  names(out)[length(by) + 1:5] <- c("coverage95", "coverage95_model", "mean_reported_sd", "mean_model_sd",
                                    "mean_abs_error")
  out
}
options(width = 200)
cat("\n== core scenarios, by kind of karyotype ==\n")
print(format(coverage(core, "kind"), digits = 3), row.names = FALSE)
core$stratum <- paste0("d", core$distance)
cat("\n== core scenarios, by scenario and distance ==\n")
print(format(coverage(core, c("scenario", "stratum")), digits = 3), row.names = FALSE)
stress <- nodes[nodes$role == "stress", ]
if (nrow(stress)) {
  stress$stratum <- paste0("d", stress$distance)
  cat("\n== stress scenarios, by scenario and distance ==\n")
  print(format(coverage(stress, c("scenario", "stratum")), digits = 3), row.names = FALSE)
}
