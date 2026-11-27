// [[Rcpp::depends(RcppEigen)]]
#include "karst_core.h"

using namespace Rcpp;

// Chromosome-level missegregation probability; the sparse joint kernel is
// built from it in R.
// [[Rcpp::export]]
double karst_pij_cpp(int parent_cn, int child_cn, double beta) {
  if (beta < 0.0 || beta > 1.0) Rcpp::stop("`beta` must be in [0, 1].");
  return karst::chromosome_transition_probability(parent_cn, child_cn, beta);
}
