// Paired variant of benchmark/sparse_two_timepoint/abm.cpp, written for
// doc/notes/nonpositive_fitness.md. Two changes, so that two runs with the
// same seed differ only by the rule for clones with fitness <= 0:
// - the population is an ordered map, and every clone draws from its own
//   random stream for each step and purpose (division, culling), so removing
//   one clone leaves the draws of all others unchanged;
// - `nonpositive_dies` switches between the benchmark's rule (no division,
//   the cells survive) and the ALFA-K simulator's (alfakR 4ad82b0: the whole
//   clone dies in that step); the cells removed are returned as
//   `died_nonpositive`.
// Everything else is the benchmark's simulator. The original header follows.
//
// Agent-based truth simulator for the sparse two-timepoint benchmark.
//
// Division and missegregation follow the ALFA-K simulations: in each step of
// length delta_t a cell divides with probability fitness * delta_t, every
// chromosome copy of a dividing cell missegregates independently with
// probability beta, and a missegregating copy goes to one daughter (+1) and is
// missing from the other (-1). When the population exceeds max_population it
// is thinned binomially, which models passaging.
//
// Two things differ from the ALFA-K simulator. Daughters with a copy number
// outside [min_cn, max_cn] die: real data contain no nullisomic karyotypes,
// and both fitting methods assume nullisomy is not viable. And the number of
// dividing cells is binomial rather than a capped Poisson, so the mean growth
// factor per step is exactly 1 + fitness * delta_t for a clone of any size.
//
// Faithful divisions are added in bulk; only divisions with at least one
// missegregation are simulated one by one. This is the same process, not an
// approximation.

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <random>
#include <string>
#include <unordered_map>
#include <map>
#include <cstdint>
#include <vector>

namespace {

typedef std::vector<int> Karyotype;

struct KaryotypeHash {
  std::size_t operator()(const Karyotype& k) const {
    std::size_t h = 1469598103934665603ULL;
    for (int v : k) {
      h ^= static_cast<std::size_t>(v + 1);
      h *= 1099511628211ULL;
    }
    return h;
  }
};

struct Clone {
  long long count;
  double fitness;
};

typedef std::map<Karyotype, Clone> Population;

double grf_fitness(const Karyotype& k, const Rcpp::NumericMatrix& centroids, double wavelength) {
  const int n_waves = centroids.nrow();
  double total = 0.0;
  for (int w = 0; w < n_waves; ++w) {
    double squared = 0.0;
    for (std::size_t c = 0; c < k.size(); ++c) {
      double diff = static_cast<double>(k[c]) - centroids(w, c);
      squared += diff * diff;
    }
    total += std::sin(std::sqrt(squared) / wavelength);
  }
  return total / (M_PI * std::sqrt(static_cast<double>(n_waves)));
}

bool viable(const Karyotype& k, int min_cn, int max_cn) {
  for (int v : k) if (v < min_cn || v > max_cn) return false;
  return true;
}

std::string label(const Karyotype& k) {
  std::string out;
  for (std::size_t c = 0; c < k.size(); ++c) {
    if (c) out.push_back('.');
    out += std::to_string(k[c]);
  }
  return out;
}


// Counter-based streams: every clone draws, at every step and for each purpose,
// from its own stream, keyed by (seed, step, karyotype, purpose). Removing one
// clone then leaves the draws of every other clone unchanged.
struct SplitMix {
  typedef std::uint64_t result_type;
  std::uint64_t state;
  explicit SplitMix(std::uint64_t s) : state(s) {}
  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return ~static_cast<result_type>(0); }
  result_type operator()() {
    std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
};
std::uint64_t mix(std::uint64_t a, std::uint64_t b) {
  SplitMix g(a ^ (b + 0x9E3779B97F4A7C15ULL + (a << 6) + (a >> 2)));
  g(); return g();
}
std::uint64_t stream_key(double seed, int step, const Karyotype& k, int purpose) {
  std::uint64_t h = static_cast<std::uint64_t>(KaryotypeHash()(k));
  return mix(mix(mix(static_cast<std::uint64_t>(seed), static_cast<std::uint64_t>(step)), h), static_cast<std::uint64_t>(purpose));
}

template <class Engine>
int truncated_error_count(int n_copies, double beta, double p_any, Engine& rng) {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double target = unif(rng) * p_any;
  double term = std::pow(1.0 - beta, n_copies);
  double cumulative = 0.0;
  for (int k = 1; k <= n_copies; ++k) {
    term *= (static_cast<double>(n_copies - k + 1) / k) * (beta / (1.0 - beta));
    cumulative += term;
    if (target <= cumulative) return k;
  }
  return n_copies;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List stt_abm_paired_cpp(Rcpp::IntegerMatrix initial_karyotypes, Rcpp::NumericVector initial_counts,
                              Rcpp::NumericMatrix centroids, double wavelength, double beta, double delta_t,
                              int n_steps, Rcpp::IntegerVector record_steps, double max_population,
                              double culling_survival, int min_cn, int max_cn, double seed, bool nonpositive_dies) {
  const int n_chr = initial_karyotypes.ncol();
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double died_nonpositive = 0.0;
  Population population;
  auto add_cells = [&](const Karyotype& k, long long n) {
    auto it = population.find(k);
    if (it == population.end()) population.emplace(k, Clone{n, grf_fitness(k, centroids, wavelength)});
    else it->second.count += n;
  };
  for (int i = 0; i < initial_karyotypes.nrow(); ++i) {
    Karyotype k(n_chr);
    for (int c = 0; c < n_chr; ++c) k[c] = initial_karyotypes(i, c);
    if (!viable(k, min_cn, max_cn)) Rcpp::stop("An initial karyotype is outside the viability bounds.");
    if (initial_counts[i] > 0) add_cells(k, static_cast<long long>(initial_counts[i]));
  }
  std::vector<int> out_step; std::vector<std::string> out_label; std::vector<double> out_count;
  auto record = [&](int step) {
    for (const auto& entry : population) {
      if (entry.second.count <= 0) continue;
      out_step.push_back(step); out_label.push_back(label(entry.first)); out_count.push_back(static_cast<double>(entry.second.count));
    }
  };
  std::size_t next_record = 0;
  if (next_record < static_cast<std::size_t>(record_steps.size()) && record_steps[next_record] == 0) { record(0); ++next_record; }
  std::vector<std::pair<Karyotype, long long>> births;
  for (int step = 1; step <= n_steps; ++step) {
    births.clear();
    for (auto& entry : population) {
      Clone& clone = entry.second;
      if (clone.count <= 0) continue;
      if (clone.fitness <= 0.0) {
        if (nonpositive_dies) { died_nonpositive += static_cast<double>(clone.count); clone.count = 0; }
        continue;
      }
      const Karyotype& parent = entry.first;
      SplitMix rng(stream_key(seed, step, parent, 1));
      double p_divide = std::min(1.0, clone.fitness * delta_t);
      long long n_divide = std::binomial_distribution<long long>(clone.count, p_divide)(rng);
      if (n_divide <= 0) continue;
      int n_copies = 0;
      for (int v : parent) n_copies += v;
      double p_any = 1.0 - std::pow(1.0 - beta, n_copies);
      long long n_erroneous = (beta > 0.0 && n_copies > 0) ? std::binomial_distribution<long long>(n_divide, p_any)(rng) : 0;
      clone.count += n_divide - n_erroneous;
      clone.count -= n_erroneous;
      for (long long d = 0; d < n_erroneous; ++d) {
        int n_errors = truncated_error_count(n_copies, beta, p_any, rng);
        Karyotype first = parent, second = parent;
        std::vector<int> chosen;
        while (static_cast<int>(chosen.size()) < n_errors) {
          int position = static_cast<int>(unif(rng) * n_copies);
          if (position >= n_copies) position = n_copies - 1;
          bool repeated = false;
          for (int previous : chosen) if (previous == position) { repeated = true; break; }
          if (!repeated) chosen.push_back(position);
        }
        for (int position : chosen) {
          int chr = 0, seen = 0;
          for (; chr < n_chr; ++chr) { seen += parent[chr]; if (position < seen) break; }
          if (unif(rng) < 0.5) { ++first[chr]; --second[chr]; } else { --first[chr]; ++second[chr]; }
        }
        if (viable(first, min_cn, max_cn)) births.emplace_back(first, 1);
        if (viable(second, min_cn, max_cn)) births.emplace_back(second, 1);
      }
    }
    for (const auto& birth : births) add_cells(birth.first, birth.second);
    double total = 0.0;
    for (auto it = population.begin(); it != population.end();) {
      if (it->second.count <= 0) it = population.erase(it);
      else { total += static_cast<double>(it->second.count); ++it; }
    }
    if (max_population > 0 && total > max_population) {
      double keep = std::min(culling_survival, max_population / total);
      for (auto it = population.begin(); it != population.end();) {
        SplitMix rng(stream_key(seed, step, it->first, 2));
        it->second.count = std::binomial_distribution<long long>(it->second.count, keep)(rng);
        if (it->second.count <= 0) it = population.erase(it); else ++it;
      }
    }
    if (next_record < static_cast<std::size_t>(record_steps.size()) && record_steps[next_record] == step) { record(step); ++next_record; }
    if (population.empty()) break;
    if (step % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(Rcpp::Named("step") = Rcpp::wrap(out_step), Rcpp::Named("karyotype") = Rcpp::wrap(out_label),
                            Rcpp::Named("count") = Rcpp::wrap(out_count), Rcpp::Named("died_nonpositive") = died_nonpositive);
}
