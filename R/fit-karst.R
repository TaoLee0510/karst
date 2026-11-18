#' Fit a karst fitness landscape
#'
#' The main entry point: fits the sparse joint count model to sparse counts at
#' two or more timepoints with [fit_sparse_joint()], whose result is returned.
#' Its outputs can be converted to those of alfak with [as_alfak()], or
#' produced directly with alfak's arguments by [alfak_compat()].
#'
#' @param counts Count matrix with one column per timepoint, at least two, or
#'   an `karst_data` object.
#' @param dt Time interval if `counts` is not already prepared; with more than
#'   two timepoints, one per pair of consecutive timepoints (or give `times`
#'   in `sparse_joint_control`).
#' @param beta Missegregation rate per chromosome copy and division.
#' @param sparse_joint_control Named list of further arguments for
#'   [fit_sparse_joint()], such as `birth_rate`, `min_cn` and `max_cn`.
#' @param ... Not used; any argument given here is an error. The two-layer and
#'   joint-count pipelines, and the `pipeline` argument that selected them,
#'   were removed in version 1.0.0.
#'
#' @return An `karst_sparse_joint_fit` object.
#' @export
fit_karst <- function(counts,
                      dt = 1,
                      beta = 0.00005,
                      sparse_joint_control = list(),
                      ...) {
  supplied <- names(list(...))
  if ("pipeline" %in% supplied) {
    stop(paste0("`pipeline` was removed in karst 1.0.0 together with the two-layer and joint-count ",
                "pipelines; fit_karst() fits the sparse joint count model."), call. = FALSE)
  }
  if (...length()) {
    if (is.null(supplied)) supplied <- rep("", ...length())
    named <- nzchar(supplied)
    refused <- c(if (any(named)) paste("no argument", paste0("`", supplied[named], "`", collapse = ", ")),
                 if (any(!named)) "no further unnamed arguments")
    stop(sprintf(paste0("fit_karst() takes %s; pass arguments of fit_sparse_joint() ",
                        "in `sparse_joint_control`."), paste(refused, collapse = " and ")),
         call. = FALSE)
  }
  if (!is.list(sparse_joint_control) || (length(sparse_joint_control) &&
      (is.null(names(sparse_joint_control)) || any(!nzchar(names(sparse_joint_control)))))) {
    stop("`sparse_joint_control` must be a named list.", call. = FALSE)
  }
  do.call(fit_sparse_joint, c(list(counts, dt = dt, beta = beta), sparse_joint_control))
}
