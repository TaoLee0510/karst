# Calibration constants as validation/sparse_joint/calibration_fit.R fits them
# (dev/PLAN.md, clarification 2): reported sd = sqrt(model_sd^2 +
# (deficit_k * relative)^2 + (floor * birth_rate * (1 + shell))^2), deficit_k =
# deficit for distance 0 (observed) and deficit * 1.1 / 0.6 otherwise. Fitted by
# Gaussian likelihood on the core datasets, every dataset and kind of karyotype
# with the same weight.

dev_calib_kind <- function(nodes) {
  total <- if (all(c("count_t0", "count_t1") %in% names(nodes))) nodes$count_t0 + nodes$count_t1 else NA
  ifelse(nodes$distance == 0L, ifelse(!is.na(total) & total >= 5, "observed, 5+ cells", "observed, <5 cells"),
         paste("unobserved, distance", nodes$distance))
}

dev_calib_sd <- function(constants, nodes, birth_rate = 1, model_sd = nodes$model_sd) {
  deficit <- ifelse(nodes$distance == 0L, constants[[2L]], constants[[2L]] * 1.1 / 0.6)
  sqrt(model_sd^2 + (deficit * nodes$relative)^2 + (constants[[1L]] * birth_rate * (1 + nodes$shell))^2)
}

# `nodes`: rows of scored nodes with dataset_id, distance, error (pred - truth),
# model_sd, relative, shell (and counts for the kind).
dev_calib_fit <- function(nodes, birth_rate = 1, model_sd = nodes$model_sd) {
  ok <- is.finite(model_sd) & is.finite(nodes$relative) & is.finite(nodes$error)
  x <- nodes[ok, ]
  ms <- model_sd[ok]
  cell <- paste(x$dataset_id, dev_calib_kind(x))
  w <- 1 / as.numeric(table(cell)[cell])
  objective <- function(log_constants) {
    v <- dev_calib_sd(exp(log_constants), x, birth_rate, ms)^2
    sum(w * (log(v) + x$error^2 / v))
  }
  best <- stats::optim(log(c(0.01, 0.5)), objective, control = list(maxit = 2000, reltol = 1e-10))
  stats::setNames(exp(best$par), c("floor", "deficit"))
}

# Mean over datasets of each dataset's coverage, by distance.
dev_calib_coverage <- function(nodes, constants, model_sd = nodes$model_sd, by = "distance") {
  sd <- dev_calib_sd(constants, nodes, 1, model_sd)
  ok <- is.finite(sd) & is.finite(nodes$error)
  x <- nodes[ok, ]
  x$covered <- abs(x$error) <= stats::qnorm(0.975) * sd[ok]
  x$width <- 2 * stats::qnorm(0.975) * sd[ok]
  per <- stats::aggregate(cbind(covered, width) ~ ., x[c("dataset_id", by, "covered", "width")], mean)
  out <- stats::aggregate(cbind(covered, width) ~ ., per[c(by, "covered", "width")], mean)
  out$worst <- stats::aggregate(covered ~ ., per[c(by, "covered")], min)$covered
  out
}
