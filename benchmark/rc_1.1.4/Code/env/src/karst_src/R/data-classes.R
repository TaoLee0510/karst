# `times` (sampling times since the first) is stored for more than two
# timepoints only, so that two-timepoint objects are as they always were.
new_karst_data <- function(counts, dt, beta, metadata = list(), times = NULL) {
  out <- list(
    counts = counts,
    labels = rownames(counts),
    dt = dt,
    beta = beta,
    metadata = metadata
  )
  if (!is.null(times)) out$times <- times
  structure(out, class = "karst_data")
}

#' @export
print.karst_data <- function(x, ...) {
  cat("<karst_data>\n")
  cat("  karyotypes:", nrow(x$counts), "\n")
  cat("  timepoints:", paste(colnames(x$counts), collapse = ", "), "\n")
  if (!is.null(x$times)) cat("  times:", paste(format(x$times), collapse = ", "), "\n")
  cat("  dt:", x$dt, " beta:", x$beta, "\n")
  invisible(x)
}
