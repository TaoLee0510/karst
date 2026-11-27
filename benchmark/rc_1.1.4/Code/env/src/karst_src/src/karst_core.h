#ifndef KARST_CORE_H
#define KARST_CORE_H

#include <RcppEigen.h>
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace karst {

inline std::string make_key(const std::vector<int>& x) {
  std::string out;
  out.reserve(x.size() * 3);
  for (size_t i = 0; i < x.size(); ++i) {
    if (i) out.push_back('.');
    out += std::to_string(x[i]);
  }
  return out;
}

// Copy numbers are unsigned decimal integers separated by single dots: no
// sign, no space, no empty field (also none at the end, which getline() would
// drop). Leading zeros are read as the number they write; callers that match
// labels as strings first rewrite them with canonical_karyotype_labels().
inline std::vector<int> parse_label(const std::string& s) {
  std::vector<int> out;
  if (s.empty()) Rcpp::stop("Malformed karyotype label: empty string.");
  if (s.back() == '.') Rcpp::stop("Malformed karyotype label '%s': empty copy-number field.", s);
  std::string token;
  std::stringstream ss(s);
  while (std::getline(ss, token, '.')) {
    if (token.empty()) Rcpp::stop("Malformed karyotype label '%s': empty copy-number field.", s);
    for (char c : token) {
      if (c < '0' || c > '9') {
        Rcpp::stop("Malformed karyotype label '%s': copy numbers must be integers.", s);
      }
    }
    char* end = nullptr;
    errno = 0;
    long v = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0' || errno == ERANGE) {
      Rcpp::stop("Malformed karyotype label '%s': copy numbers must be integers.", s);
    }
    if (v < 0 || v > std::numeric_limits<int>::max()) {
      Rcpp::stop("Malformed karyotype label '%s': copy number outside integer range.", s);
    }
    out.push_back(static_cast<int>(v));
  }
  if (out.empty()) Rcpp::stop("Malformed karyotype label: empty string.");
  return out;
}

inline Rcpp::IntegerMatrix labels_to_matrix(const Rcpp::CharacterVector& labels) {
  int n = labels.size();
  if (n == 0) Rcpp::stop("At least one karyotype label is required.");
  std::vector<int> first = parse_label(Rcpp::as<std::string>(labels[0]));
  int p = first.size();
  Rcpp::IntegerMatrix mat(n, p);
  for (int j = 0; j < p; ++j) mat(0, j) = first[j];
  std::unordered_set<std::string> seen;
  seen.insert(make_key(first));
  for (int i = 1; i < n; ++i) {
    std::vector<int> x = parse_label(Rcpp::as<std::string>(labels[i]));
    if (static_cast<int>(x.size()) != p) {
      Rcpp::stop("All karyotypes must have the same chromosome dimension.");
    }
    std::string key = make_key(x);
    if (!seen.insert(key).second) Rcpp::stop("Duplicate karyotype label: %s", key);
    for (int j = 0; j < p; ++j) mat(i, j) = x[j];
  }
  return mat;
}

inline Rcpp::CharacterVector matrix_to_labels(const Rcpp::IntegerMatrix& mat) {
  int n = mat.nrow(), p = mat.ncol();
  Rcpp::CharacterVector labels(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int> x(p);
    for (int j = 0; j < p; ++j) x[j] = mat(i, j);
    labels[i] = make_key(x);
  }
  return labels;
}

inline double choose_double(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  if (k == 0 || k == n) return 1.0;
  k = std::min(k, n - k);
  double out = 1.0;
  for (int i = 1; i <= k; ++i) {
    out *= static_cast<double>(n - k + i) / i;
  }
  return out;
}

inline double chromosome_transition_probability(int parent_cn, int child_cn, double beta) {
  if (parent_cn < 0 || child_cn < 0) return 0.0;
  int diff = std::abs(parent_cn - child_cn);
  if (diff > parent_cn) return 0.0;
  double total = 0.0;
  for (int m = diff; m <= parent_cn; m += 2) {
    int lost_to_child = (m + parent_cn - child_cn) / 2;
    if (lost_to_child < 0 || lost_to_child > m) continue;
    total += choose_double(parent_cn, m) *
      std::pow(beta, m) *
      std::pow(1.0 - beta, parent_cn - m) *
      std::pow(0.5, m) *
      choose_double(m, lost_to_child);
  }
  return total;
}

} // namespace karst

#endif
