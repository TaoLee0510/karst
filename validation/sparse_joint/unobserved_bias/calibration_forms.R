#!/usr/bin/env Rscript
# Forms of the reported uncertainty, fitted to the pooled default fits of
# calibration.R (results/calibration_nodes.rds) and compared by coverage. Every
# form is the package's (sparse_joint_calibrate() in R/sparse-joint.R):
#
#   sd^2 = model^2 + (deficit_k * relative)^2 + (floor * (1 + shell))^2
#
# with `relative` the fit's own relative fitness (relative to the count-weighted
# mean of its major karyotypes, not the scoring reference) and deficit_k the
# deficit of observed karyotypes or, for karyotypes seen at neither time point,
# deficit_unobserved:
#
#   shipped   the package's constants: floor 0.009, deficit 0.6,
#             deficit_unobserved 1.1
#   before    floor 0.009 and deficit 0.6 for all karyotypes, the constants
#             before 2026-11-13 (and those the pooled fits were made with)
#   refit     floor and one deficit for all, refitted
#   split     floor, deficit and deficit_unobserved, refitted
#
# Then the candidates for deficit_unobserved (0.6, 1.0, 1.1, 1.2; floor 0.009
# and deficit 0.6), the table behind the package's 1.1 (README, decision of
# 2026-11-13): over the cells of unobserved karyotypes (set, wavelength and
# distance 1 or 2), the lowest coverage and the number below 0.93 outside the
# stress scenarios, the lowest in the stress scenarios, and the median width of
# unobserved karyotypes' intervals relative to deficit 0.6 outside them.
#
# The baseline birth rate is 1 in every fit here, so the floor's factor
# birth_rate is left out. Constants are fitted by Gaussian likelihood with
# every dataset, and within it every kind of karyotype, weighted equally, on
# the sets named by --fit (default: dev core, diagnosis, rugged; the stress
# scenarios, where the landscape changes, are reported only, as before).
# --leave-out=<set> refits without one set and reports it.
#
#   Rscript validation/sparse_joint/unobserved_bias/calibration_forms.R [--leave-out=dev_core|diagnosis|rugged]

if (!exists("ub_dir")) ub_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})

# Constants the pooled fits were made with (af362af): one deficit for all.
ub_fit_calibration <- c(floor = 0.009, deficit = 0.6)

# The fit's relative fitness, which the package's deficit multiplies. Tables
# collected with a `relative` column (the prediction's `fitness`) use it.
# Otherwise its size is recovered from the fit's calibrated and model standard
# deviations, which differ by exactly the deficit and floor terms of the
# constants the fits were made with; only its size enters the formula.
ub_relative <- function(x, calibration = ub_fit_calibration) {
  if (!is.null(x$relative)) return(x$relative)
  deficit_term <- x$fitness_sd^2 - x$model_sd^2 - (calibration[["floor"]] * (1 + x$shell))^2
  sqrt(pmax(deficit_term, 0)) / calibration[["deficit"]]
}

ub_load_calibration_nodes <- function(dir = ub_dir) {
  x <- readRDS(file.path(dir, "results", "calibration_nodes.rds"))
  x$error <- x$pred - x$truth
  x$set <- ifelse(x$source == "dev", paste0("dev_", x$role), x$source)
  x$observed <- x$distance == 0L
  x$relative <- ub_relative(x)
  cell <- paste(x$dataset_id, x$kind)
  x$weight <- 1 / as.numeric(table(cell)[cell])
  x
}

# `p`: floor, deficit, deficit_unobserved.
ub_package_sd <- function(p, x) {
  deficit <- ifelse(x$observed, p[[2L]], p[[3L]])
  sqrt(x$model_sd^2 + (deficit * x$relative)^2 + (p[[1L]] * (1 + x$shell))^2)
}

ub_forms <- list(
  shipped = list(start = numeric(0), sd = function(p, x) ub_package_sd(c(0.009, 0.6, 1.1), x),
                 names = character(0)),
  before = list(start = numeric(0), sd = function(p, x) ub_package_sd(c(0.009, 0.6, 0.6), x),
                names = character(0)),
  refit = list(start = log(c(0.009, 0.6)), sd = function(p, x) ub_package_sd(c(p[[1L]], p[[2L]], p[[2L]]), x),
               names = c("floor", "deficit")),
  split = list(start = log(c(0.009, 0.6, 0.6)), sd = function(p, x) ub_package_sd(p, x),
               names = c("floor", "deficit_observed", "deficit_unobserved"))
)

ub_fit_form <- function(form, x) {
  if (!length(form$start)) return(numeric(0))
  objective <- function(log_p) {
    v <- form$sd(exp(log_p), x)^2
    sum(x$weight * (log(v) + x$error^2 / v))
  }
  best <- stats::optim(form$start, objective, control = list(maxit = 4000, reltol = 1e-12))
  stats::setNames(exp(best$par), form$names)
}

# Dataset-weighted coverage and mean reported sd, by set, wavelength and kind.
ub_coverage <- function(form, p, x) {
  x$sd <- form$sd(p, x)
  x$covered <- abs(x$error) <= stats::qnorm(0.975) * x$sd
  per <- stats::aggregate(cbind(covered, sd) ~ set + wavelength_class + kind + dataset_id, x, mean)
  stats::aggregate(cbind(covered, sd) ~ set + wavelength_class + kind, per, mean)
}

ub_wide <- function(table, value) {
  w <- stats::reshape(table[c("set", "wavelength_class", "kind", value)], idvar = c("set", "wavelength_class"),
                      timevar = "kind", direction = "wide")
  names(w) <- sub(paste0(value, "."), "", names(w), fixed = TRUE)
  w[order(w$set, w$wavelength_class), ]
}

# One row per candidate deficit for unobserved karyotypes (see the header).
ub_unobserved_candidates <- function(x, candidates = c(0.6, 1.0, 1.1, 1.2), floor = 0.009, deficit = 0.6,
                                     stress_sets = "dev_stress") {
  reference <- ub_package_sd(c(floor, deficit, deficit), x)
  unobserved <- !x$observed & !x$set %in% stress_sets
  do.call(rbind, lapply(candidates, function(k) {
    p <- c(floor, deficit, k)
    form <- list(sd = function(p, x) ub_package_sd(p, x))
    table <- ub_coverage(form, p, x)
    table <- table[grepl("^unobserved", table$kind), ]
    stress <- table$set %in% stress_sets
    data.frame(deficit_unobserved = k,
               worst_non_stress = min(table$covered[!stress]),
               cells_below_0.93 = sum(table$covered[!stress] < 0.93), cells = sum(!stress),
               worst_stress = if (any(stress)) min(table$covered[stress]) else NA_real_,
               width_median = stats::median(ub_package_sd(p, x)[unobserved] / reference[unobserved]))
  }))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  leave_out <- sub("^--leave-out=", "", c(grep("^--leave-out=", args, value = TRUE), "--leave-out=")[[1L]])
  x <- ub_load_calibration_nodes()
  fit_sets <- setdiff(c("dev_core", "diagnosis", "rugged"), leave_out)
  train <- x[x$set %in% fit_sets, ]
  options(width = 220)
  for (name in names(ub_forms)) {
    form <- ub_forms[[name]]
    p <- ub_fit_form(form, train)
    cat(sprintf("\n==== %s %s (fitted on %s)\n", name,
                paste(sprintf("%s=%.4f", names(p), p), collapse = " "), paste(fit_sets, collapse = ", ")))
    table <- ub_coverage(form, p, x)
    table$covered <- round(table$covered, 3)
    print(ub_wide(table, "covered"), row.names = FALSE)
  }
  cat("\n==== deficit for unobserved karyotypes (floor 0.009, deficit 0.6)\n")
  print(format(ub_unobserved_candidates(x), digits = 3), row.names = FALSE)
}
