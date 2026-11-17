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

new_karst_graph <- function(x) {
  class(x) <- "karst_graph"
  x
}

new_karst_local_fit <- function(x) {
  class(x) <- "karst_local_fit"
  x
}

new_karst_global_fit <- function(x) {
  class(x) <- c("karst_global_fit", setdiff(class(x), "karst_global_fit"))
  x
}

new_karst_fit <- function(x) {
  class(x) <- "karst_fit"
  x
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

#' @export
print.karst_graph <- function(x, ...) {
  cat("<karst_graph>\n")
  cat("  nodes:", length(x$labels), "\n")
  cat("  chromosomes:", x$n_chr, "\n")
  cat("  shell_depth:", x$shell_depth, "\n")
  invisible(x)
}

#' @export
print.karst_fit <- function(x, ...) {
  cat("<karst_fit>\n")
  cat("  local nodes:", nrow(x$local$summary), "\n")
  cat("  global nodes:", nrow(x$global$summary), "\n")
  cat("  convergence:", x$local$diagnostics$convergence, "\n")
  if (!is.null(x$xval$R2R)) {
    cat("  xval R2R:", x$xval$R2R, "\n")
  }
  invisible(x)
}
