#' KARST: Karyotype Abundance-Based Reconstruction of Fitness Landscape from Sparse Time Points
#'
#' `karst` implements KARST (Karyotype Abundance-based Reconstruction of
#' fitness landscape from Sparse Time points). It infers fitness landscapes of
#' aneuploid karyotypes from sparse single-cell karyotype counts at two or more
#' time points with the sparse joint count model ([fit_karst()],
#' [fit_sparse_joint()]). KARST was released as alfak2 up to version 0.2.0.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib karst, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(".data"))
