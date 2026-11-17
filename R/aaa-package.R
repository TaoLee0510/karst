#' KARST: Karyotype Abundance-Based Reconstruction from Sparse Time Points
#'
#' `karst` implements KARST (Karyotype Abundance-based Reconstruction from
#' Sparse Time points). It infers fitness landscapes of aneuploid karyotypes
#' from sparse single-cell karyotype counts at two or more time points; the
#' default model is the sparse joint count model ([fit_sparse_joint()]). The
#' earlier two-layer pipeline, which builds bounded karyotype graphs, fits a
#' local hierarchical Bayesian model and propagates estimates across the graph
#' with a Gaussian-process prior, remains available. KARST was released as
#' alfak2 up to version 0.2.0.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib karst, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(".data"))
