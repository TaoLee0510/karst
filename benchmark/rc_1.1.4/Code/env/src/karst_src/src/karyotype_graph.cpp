// [[Rcpp::depends(RcppEigen)]]
#include "karst_core.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerMatrix karst_parse_karyotypes_cpp(Rcpp::CharacterVector labels) {
  return karst::labels_to_matrix(labels);
}

// [[Rcpp::export]]
Rcpp::CharacterVector karst_stringify_karyotypes_cpp(Rcpp::IntegerMatrix karyotypes) {
  return karst::matrix_to_labels(karyotypes);
}
