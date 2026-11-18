#' Prepare karyotype count data
#'
#' @param counts Matrix or data frame with one column per timepoint, at least
#'   two, in time order, and dot-separated karyotype labels as row names.
#' @param dt With two timepoints, the positive time between them. With more,
#'   either `times` or the positive time between each pair of consecutive
#'   timepoints (one value fewer than there are timepoints).
#' @param beta Missegregation rate per chromosome copy and division.
#' @param metadata Optional list stored on the returned object.
#' @param times Optional sampling times, one per column, strictly increasing,
#'   instead of `dt`. Only differences matter.
#'
#' @return An `karst_data` object. `dt` is the time from the first to the last
#'   timepoint. With more than two timepoints, `times` holds the time of each
#'   since the first; [karst_data_times()] returns it for any number.
#' @export
prepare_karst_data <- function(counts, dt = 1, beta = 0.00005, metadata = list(), times = NULL) {
  counts <- validate_count_matrix(counts)
  n_time <- ncol(counts)
  colnames(counts) <- paste0("t", seq_len(n_time) - 1L)
  times <- resolve_sample_times(n_time, dt, times, dt_given = !missing(dt))
  validate_scalar(beta, "beta", lower = 0, upper = 1)
  if (is.null(metadata)) metadata <- list()
  if (!is.list(metadata)) stop("`metadata` must be a list.", call. = FALSE)
  parse_karyotypes(rownames(counts))
  new_karst_data(counts = counts, dt = times[[n_time]], beta = beta, metadata = metadata,
                  times = if (n_time > 2L) times else NULL)
}

# Sampling times since the first timepoint, from `times` or from `dt`: one
# interval for two timepoints, one per pair of consecutive timepoints for
# more. Two timepoints give c(0, dt) exactly.
resolve_sample_times <- function(n_time, dt, times, dt_given) {
  if (!is.null(times)) {
    if (dt_given) stop("Give `dt` or `times`, not both.", call. = FALSE)
    if (!is.numeric(times) || length(times) != n_time || any(!is.finite(times))) {
      stop(sprintf("`times` must give %d finite times, one per timepoint.", n_time), call. = FALSE)
    }
    if (any(diff(times) <= 0)) stop("`times` must be strictly increasing.", call. = FALSE)
    return(as.numeric(times) - times[[1L]])
  }
  if (n_time == 2L) {
    validate_scalar(dt, "dt", lower = .Machine$double.eps)
    return(c(0, dt))
  }
  if (!is.numeric(dt) || length(dt) != n_time - 1L || any(!is.finite(dt)) || any(dt <= 0)) {
    stop(sprintf(paste0("With %d timepoints give `times`, or `dt` as the %d positive intervals between ",
                        "consecutive timepoints."), n_time, n_time - 1L), call. = FALSE)
  }
  c(0, cumsum(as.numeric(dt)))
}

#' Sampling times of count data
#'
#' @param data An `karst_data` object.
#' @return The time of each timepoint since the first: `c(0, dt)` for two
#'   timepoints.
#' @export
karst_data_times <- function(data) {
  if (!inherits(data, "karst_data")) stop("`data` must be an `karst_data` object.", call. = FALSE)
  if (is.null(data$times)) c(0, data$dt) else data$times
}
