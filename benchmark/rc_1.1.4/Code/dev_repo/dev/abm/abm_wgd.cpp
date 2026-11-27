// dev copy of validation/sparse_joint/series_methods/abm_constant.cpp (itself a
// copy of benchmark/sparse_two_timepoint/abm.cpp). Two additions, neither of
// which draws a random number, so that a run with landscape_type 0 replays
// the original simulator exactly:
//  - every passage (binomial thinning) is logged: step, cells before, the
//    survival probability, cells after;
//  - landscape_type 1 (rough Mount Fuji) and 2 (NK on the copy-number lattice)
//    give other fitness landscapes, from a hash of the karyotype, so that
//    they need no table. dev_landscape_fitness() gives the same values in R.
// `constant_fitness` > 0 gives every viable karyotype that division rate.
// Exported as dev_abm_cpp so that it does not replace stt_abm_cpp or sm_abm_cpp.
//
// R10 copy (dev/phase3/r10_spec.md): exported as wgd_abm_cpp and
// wgd_landscape_fitness. Additions: whole-genome doubling (a division becomes,
// with probability `omega`, one cell with karyotype 2k; it dies if 2k leaves
// [min_cn, max_cn]); cells of ploidy 4 (mean copy number >= 3) missegregate
// at `beta_factor` x beta per copy; landscape_type 0 is ploidy-aware: the GRF
// at k * 2 / ploidy(k), plus `wgd_gain` for ploidy 4.
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

typedef std::unordered_map<Karyotype, Clone, KaryotypeHash> Population;

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

// splitmix64: a deterministic 64-bit mix, for landscape noise without tables.
unsigned long long mix64(unsigned long long x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

double hash_unit(unsigned long long h) {
  return (static_cast<double>(h >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

// Standard normal from a hash of the integers in `key` and a seed (Box-Muller).
double hash_normal(const std::vector<long long>& key, unsigned long long seed) {
  unsigned long long h = mix64(seed);
  for (long long v : key) h = mix64(h ^ static_cast<unsigned long long>(v + 1000003LL));
  double u1 = hash_unit(h), u2 = hash_unit(mix64(h ^ 0xA5A5A5A5ULL));
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

// landscape_type 1, rough Mount Fuji: par = (top, slope per copy, noise sd,
// noise seed); ref = the summit karyotype. f = top - slope * |k - ref|_1 +
// sd * z(k).
// landscape_type 2, NK: par = (mean, sd, K, seed). Chromosome c contributes
// z(c, k_c, k_{c+1}, ..., k_{c+K}) (indices mod n), and f = mean + sd *
// sum_c z / sqrt(n).
double other_fitness(const Karyotype& k, int type, const Rcpp::NumericVector& par, const Rcpp::IntegerVector& ref) {
  const int n = static_cast<int>(k.size());
  if (type == 1) {
    double distance = 0.0;
    for (int c = 0; c < n; ++c) distance += std::abs(k[c] - ref[c]);
    std::vector<long long> key(k.begin(), k.end());
    return par[0] - par[1] * distance + par[2] * hash_normal(key, static_cast<unsigned long long>(par[3]));
  }
  if (type == 2) {
    const int K = static_cast<int>(par[2]);
    double total = 0.0;
    for (int c = 0; c < n; ++c) {
      std::vector<long long> key;
      key.push_back(c);
      for (int j = 0; j <= K; ++j) key.push_back(k[(c + j) % n]);
      total += hash_normal(key, static_cast<unsigned long long>(par[3]));
    }
    return par[0] + par[1] * total / std::sqrt(static_cast<double>(n));
  }
  Rcpp::stop("Unknown landscape_type.");
  return 0.0;
}

int ploidy(const Karyotype& k) {
  double total = 0.0;
  for (int v : k) total += v;
  return (total / static_cast<double>(k.size()) >= 3.0) ? 4 : 2;
}

double grf_fitness_scaled(const Karyotype& k, const Rcpp::NumericMatrix& centroids, double wavelength, double gain) {
  const int n_waves = centroids.nrow();
  const double factor = 2.0 / ploidy(k);
  double total = 0.0;
  for (int w = 0; w < n_waves; ++w) {
    double squared = 0.0;
    for (std::size_t c = 0; c < k.size(); ++c) {
      double diff = static_cast<double>(k[c]) * factor - centroids(w, c);
      squared += diff * diff;
    }
    total += std::sin(std::sqrt(squared) / wavelength);
  }
  return total / (M_PI * std::sqrt(static_cast<double>(n_waves))) + (ploidy(k) == 4 ? gain : 0.0);
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

// Number of missegregating copies in a division, given that there is at least
// one: a zero-truncated Binomial(n_copies, beta) drawn by inverse CDF.
int truncated_error_count(int n_copies, double beta, double p_any, std::mt19937_64& rng) {
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  double target = unif(rng) * p_any;
  double term = std::pow(1.0 - beta, n_copies);  // P(k = 0)
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
Rcpp::List wgd_abm_cpp(Rcpp::IntegerMatrix initial_karyotypes,
                       Rcpp::NumericVector initial_counts,
                       Rcpp::NumericMatrix centroids,
                       double wavelength,
                       double beta,
                       double delta_t,
                       int n_steps,
                       Rcpp::IntegerVector record_steps,
                       double max_population,
                       double culling_survival,
                       int min_cn,
                       int max_cn,
                       double seed,
                       double constant_fitness,
                       int landscape_type = 0,
                       Rcpp::NumericVector landscape_par = Rcpp::NumericVector::create(),
                       Rcpp::IntegerVector landscape_ref = Rcpp::IntegerVector::create(),
                       double omega = 0.0,
                       double beta_factor = 2.5,
                       double wgd_gain = 0.0) {
  const int n_chr = initial_karyotypes.ncol();
  if (centroids.ncol() != n_chr) Rcpp::stop("Centroids and karyotypes differ in dimension.");
  if (beta < 0.0 || beta >= 1.0) Rcpp::stop("`beta` must be in [0, 1).");
  std::mt19937_64 rng(static_cast<unsigned long long>(seed));
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  Population population;
  auto add_cells = [&](const Karyotype& k, long long n) {
    auto it = population.find(k);
    if (it == population.end()) {
      double f = constant_fitness > 0.0 ? constant_fitness :
        (landscape_type == 0 ? grf_fitness_scaled(k, centroids, wavelength, wgd_gain) : other_fitness(k, landscape_type, landscape_par, landscape_ref));
      population.emplace(k, Clone{n, f});
    } else {
      it->second.count += n;
    }
  };
  for (int i = 0; i < initial_karyotypes.nrow(); ++i) {
    Karyotype k(n_chr);
    for (int c = 0; c < n_chr; ++c) k[c] = initial_karyotypes(i, c);
    if (!viable(k, min_cn, max_cn)) Rcpp::stop("An initial karyotype is outside the viability bounds.");
    if (initial_counts[i] > 0) add_cells(k, static_cast<long long>(initial_counts[i]));
  }

  std::vector<int> passage_step;
  std::vector<double> passage_before, passage_keep, passage_after;
  std::vector<int> out_step;
  std::vector<std::string> out_label;
  std::vector<double> out_count;
  auto record = [&](int step) {
    for (const auto& entry : population) {
      if (entry.second.count <= 0) continue;
      out_step.push_back(step);
      out_label.push_back(label(entry.first));
      out_count.push_back(static_cast<double>(entry.second.count));
    }
  };
  std::size_t next_record = 0;
  if (next_record < static_cast<std::size_t>(record_steps.size()) && record_steps[next_record] == 0) {
    record(0);
    ++next_record;
  }

  std::vector<std::pair<Karyotype, long long>> births;
  for (int step = 1; step <= n_steps; ++step) {
    births.clear();
    for (auto& entry : population) {
      Clone& clone = entry.second;
      if (clone.count <= 0 || clone.fitness <= 0.0) continue;  // no division at non-positive fitness
      const Karyotype& parent = entry.first;
      double p_divide = std::min(1.0, clone.fitness * delta_t);
      long long n_divide = std::binomial_distribution<long long>(clone.count, p_divide)(rng);
      if (n_divide <= 0) continue;
      int n_copies = 0;
      for (int v : parent) n_copies += v;
      // R10: whole-genome doublings among the divisions (one cell 2k each)
      if (omega > 0.0) {
        long long n_wgd = std::binomial_distribution<long long>(n_divide, omega)(rng);
        if (n_wgd > 0) {
          n_divide -= n_wgd;
          clone.count -= n_wgd;
          Karyotype doubled = parent;
          for (int& v : doubled) v *= 2;
          if (viable(doubled, min_cn, max_cn)) births.emplace_back(doubled, n_wgd);
        }
      }
      double beta_k = (ploidy(parent) == 4) ? beta * beta_factor : beta;
      if (beta_k >= 1.0) beta_k = 0.999;
      double p_any = 1.0 - std::pow(1.0 - beta_k, n_copies);
      long long n_erroneous = (beta_k > 0.0 && n_copies > 0)
        ? std::binomial_distribution<long long>(n_divide, p_any)(rng) : 0;
      // A faithful division replaces one parent by two: net +1.
      clone.count += n_divide - n_erroneous;
      // An erroneous division removes the parent and adds its viable daughters.
      clone.count -= n_erroneous;
      for (long long d = 0; d < n_erroneous; ++d) {
        int n_errors = truncated_error_count(n_copies, beta_k, p_any, rng);
        Karyotype first = parent, second = parent;
        std::vector<int> chosen;
        chosen.reserve(n_errors);
        while (static_cast<int>(chosen.size()) < n_errors) {
          int position = static_cast<int>(unif(rng) * n_copies);
          if (position >= n_copies) position = n_copies - 1;
          bool repeated = false;
          for (int previous : chosen) if (previous == position) { repeated = true; break; }
          if (!repeated) chosen.push_back(position);
        }
        for (int position : chosen) {
          int chr = 0, seen = 0;
          for (; chr < n_chr; ++chr) {
            seen += parent[chr];
            if (position < seen) break;
          }
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
      double after = 0.0;
      for (auto it = population.begin(); it != population.end();) {
        it->second.count = std::binomial_distribution<long long>(it->second.count, keep)(rng);
        if (it->second.count <= 0) it = population.erase(it); else { after += static_cast<double>(it->second.count); ++it; }
      }
      passage_step.push_back(step);
      passage_before.push_back(total);
      passage_keep.push_back(keep);
      passage_after.push_back(after);
    }
    if (next_record < static_cast<std::size_t>(record_steps.size()) && record_steps[next_record] == step) {
      record(step);
      ++next_record;
    }
    if (population.empty()) break;
    if (step % 256 == 0) Rcpp::checkUserInterrupt();
  }
  return Rcpp::List::create(
    Rcpp::Named("step") = Rcpp::wrap(out_step),
    Rcpp::Named("karyotype") = Rcpp::wrap(out_label),
    Rcpp::Named("count") = Rcpp::wrap(out_count),
    Rcpp::Named("passages") = Rcpp::DataFrame::create(
      Rcpp::Named("step") = Rcpp::wrap(passage_step),
      Rcpp::Named("before") = Rcpp::wrap(passage_before),
      Rcpp::Named("keep") = Rcpp::wrap(passage_keep),
      Rcpp::Named("after") = Rcpp::wrap(passage_after))
  );
}

// Fitness of karyotypes (one row each) on a landscape of type 1 or 2, as
// dev_abm_cpp() assigns it.
// [[Rcpp::export]]
Rcpp::NumericVector wgd_landscape_fitness_other(Rcpp::IntegerMatrix karyotypes, int landscape_type,
                                          Rcpp::NumericVector landscape_par, Rcpp::IntegerVector landscape_ref) {
  Rcpp::NumericVector out(karyotypes.nrow());
  for (int i = 0; i < karyotypes.nrow(); ++i) {
    Karyotype k(karyotypes.ncol());
    for (int c = 0; c < karyotypes.ncol(); ++c) k[c] = karyotypes(i, c);
    out[i] = other_fitness(k, landscape_type, landscape_par, landscape_ref);
  }
  return out;
}


// Fitness of karyotypes (one row each) on the ploidy-aware GRF of landscape_type 0.
// [[Rcpp::export]]
Rcpp::NumericVector wgd_landscape_fitness(Rcpp::IntegerMatrix karyotypes, Rcpp::NumericMatrix centroids,
                                          double wavelength, double wgd_gain) {
  Rcpp::NumericVector out(karyotypes.nrow());
  for (int i = 0; i < karyotypes.nrow(); ++i) {
    Karyotype k(karyotypes.ncol());
    for (int c = 0; c < karyotypes.ncol(); ++c) k[c] = karyotypes(i, c);
    out[i] = grf_fitness_scaled(k, centroids, wavelength, wgd_gain);
  }
  return out;
}
