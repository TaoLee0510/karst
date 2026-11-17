#' Summarize karst fits
#'
#' @param object An `karst_fit`, `karst_local_fit`, or graph posterior list.
#' @param layer For full fits, `"global"` or `"local"`.
#' @param ... Reserved.
#'
#' @return Data frame of posterior node summaries.
#' @export
summarize_karst <- function(object, layer = c("global", "local"), ...) {
  layer <- match.arg(layer)
  if (inherits(object, "karst_fit")) {
    if (layer == "global") return(object$global$summary)
    return(object$local$summary)
  }
  if (inherits(object, "karst_local_fit")) return(object$summary)
  if (is.list(object) && !is.null(object$summary)) return(object$summary)
  stop("Unsupported object type for summarize_karst().", call. = FALSE)
}

#' @export
summary.karst_fit <- function(object, ...) {
  summarize_karst(object, ...)
}

#' @export
summary.karst_local_fit <- function(object, ...) {
  summarize_karst(object, ...)
}

#' @export
summary.karst_global_fit <- function(object, ...) {
  summarize_karst(object, ...)
}
