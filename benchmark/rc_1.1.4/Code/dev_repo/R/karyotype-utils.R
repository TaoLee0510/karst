#' Parse karyotype labels
#'
#' Converts dot-separated copy-number labels such as `"2.2.3.1"` to an integer
#' matrix. Parsing is compiled and validates dimensional consistency.
#'
#' @param labels Character vector of karyotype labels.
#' @return Integer matrix with one row per label.
#' @keywords internal
parse_karyotypes <- function(labels) {
  karst_parse_karyotypes_cpp(as.character(labels))
}

#' Format karyotype matrix rows
#'
#' @param karyotypes Integer matrix.
#' @return Character vector of dot-separated labels.
#' @keywords internal
format_karyotypes <- function(karyotypes) {
  karst_stringify_karyotypes_cpp(as.matrix(karyotypes))
}

# Labels as format_karyotypes() writes them, so that labels can be matched as
# strings: "2.02.3" becomes "2.2.3". A label must be unsigned integers
# separated by single dots, as parse_karyotypes() requires; anything else
# ("2.2.", " 2.2", "+2.2") is an error. Duplicates are kept.
canonical_karyotype_labels <- function(labels) {
  labels <- as.character(labels)
  bad <- is.na(labels) | !grepl("^[0-9]+(\\.[0-9]+)*$", labels)
  if (any(bad)) {
    stop(sprintf("Malformed karyotype label '%s': copy numbers must be integers separated by single dots.",
                 labels[bad][[1L]]), call. = FALSE)
  }
  gsub("(^|\\.)0+([0-9])", "\\1\\2", labels)
}
