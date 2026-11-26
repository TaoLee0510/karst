# Fitness of a karyotype from its own counts: the targets of the
# cross-validation with as_alfak(xval_target = "own").
#
# alfak's cross-validation (alfakR dev2, xval()) compares the Kriging
# prediction for a held-out karyotype with an estimate of its fitness from its
# own counts: one bootstrap draw of its first-stage fitness. The default
# cross-validation of as_alfak() compares the held-out prediction with the
# model's estimate from all the data instead, which comes from the same
# smooth model and lies close to it. The counterpart of alfak's target here:
# every latent variable of the fit except the karyotype's own stays at the
# fit's mode, and the karyotype's own terms of the joint objective are
# minimised over its log birth rate and its own latent abundances, without the
# landscape prior on its birth rate. Its own terms are
#   - its counts at every time point, as the fit uses them (weights included);
#   - its initial abundance term: the balance prior of a non-major state fed
#     by majors, otherwise the prior of a free initial abundance (pairwise
#     method: one per pair it belongs to);
#   - in the state-space method, its own transitions from time point to time
#     point.
# Its neighbours enter only through the influx they send it, at the mode; the
# reference state's growth over each interval and the depth offsets stay as
# fitted (also for the reference state itself), and the karyotype's outflow to
# its neighbours and its part in their balance prior are left out. With
# `resample`, the counts are first drawn once again cell by cell, at each
# time point Binomial(cells at the time point, count / cells), as alfak's
# bootstrap resamples cells; the draw is made for every state of the domain
# at once, so a state's draw does not depend on which states are asked for.
#
# The birth rate is searched within [min(g) - 6, max(g) + 3] of the fit's log
# birth rates (a factor of about 400 below the lowest and 20 above the
# highest) and latent log abundances within [-30, log(cells) + 5]. A state
# that the draw leaves without cells ends at the lower bound of its birth
# rate: no growth, the least fit it can be.

# log phi1(a, b; T) as the TMB templates compute it (sparse_joint_log_phi1():
# the rate difference smoothed at zero), so that a state's terms here equal
# those of the fit's objective.
sparse_joint_own_log_phi1 <- function(a, b, dt) {
  gap <- (a - b) * dt
  x <- sqrt(gap * gap + 1e-12)
  0.5 * (a + b) * dt + 0.5 * x + log(dt) + log(-expm1(-x) / x)
}

sparse_joint_own_lse <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

# Poisson log density of `y` with log intensity `l` (y integer; 0 counts
# at intensity 0 have density 1).
sparse_joint_own_pois <- function(y, l) (if (y > 0) y * l else 0) - exp(l) - lgamma(y + 1)

# Everything a state's own terms need, from a fit: the TMB data and the
# latent variables at the fit's mode (as sparse_joint_rebuild() assembles
# them), the edges into each state, the counts as states x time points and
# the cells at each time point.
sparse_joint_own_setup <- function(fit) {
  if (inherits(fit, "karst_sparse_joint_rate_fit")) fit <- fit$fits[[fit$point]]
  if (is.null(fit$settings) || is.null(fit$posterior)) {
    stop("This fit predates stored settings; refit it with the current karst.", call. = FALSE)
  }
  problem <- sparse_joint_problem(fit$data, fit$settings)
  td <- problem$tmb_data
  models <- c("sparse_joint", "sparse_joint_series", "sparse_joint_pairwise", "sparse_joint_state_space")
  if (!td$model %in% models) stop("Own-count targets are not available for model ", td$model, ".", call. = FALSE)
  par <- problem$parameters
  stored <- intersect(names(par), names(fit$posterior))
  par[stored] <- fit$posterior[stored]
  n <- length(td$weight)
  # Matrices of latent abundances (pairwise: states x pairs; state space:
  # states x later time points) may be stored as vectors.
  if (td$model == "sparse_joint_pairwise") par$zeta <- matrix(par$zeta, nrow = n)
  if (td$model == "sparse_joint_state_space") par$zeta_later <- matrix(par$zeta_later, nrow = n)
  y <- if (td$model == "sparse_joint") cbind(td$y0, td$y1) else as.matrix(td$y)
  from <- td$edge_from + 1L
  to <- td$edge_to + 1L
  self <- 2 * exp(td$log_self) - 1
  cells <- colSums(as.matrix(fit$data$counts))
  if (length(cells) != ncol(y)) stop("the counts of the fit and of its domain have different time points", call. = FALSE)
  # Counts of each state that enter the likelihood (sample_weight of the
  # anchored method: time points in between that count for the majors only).
  used <- matrix(td$weight, n, ncol(y))
  if (!is.null(td$sample_weight)) used[, -1L] <- used[, -1L] * td$sample_weight
  list(model = td$model, td = td, par = par, n = n, y = y, cells = cells, used = used, labels = problem$domain$labels,
       from = from, into = split(seq_along(to), factor(to, levels = seq_len(n))), self = self,
       net = self * exp(par$g), n_major = td$n_major,
       g_bounds = c(min(par$g) - 6, max(par$g) + 3), z_bounds = c(-30, log(max(cells)) + 5))
}

# Initial abundance term of state s with log abundance z and net rate ns:
# the balance prior when it is fed by majors (two time points, anchored,
# state space), otherwise the prior of a free abundance.
sparse_joint_own_init <- function(st, s, z, ns) {
  td <- st$td
  if (td$balance == 1L && s > st$n_major) {
    e <- st$into[[s]]
    e <- e[st$from[e] <= st$n_major]
    if (length(e)) {
      p <- st$from[e]
      start <- sparse_joint_own_lse(log(2) + td$edge_log_prob[e] + st$par$g[p] + st$par$zeta[p] -
                                      st$net[p] * td$history + sparse_joint_own_log_phi1(ns, st$net[p], td$history))
      level <- if (st$model == "sparse_joint" && !is.null(st$par$balance_level)) st$par$balance_level else 0
      residual <- z - start - level
      return(-stats::dnorm(residual - td$balance_offset[s], 0, td$balance_sd, log = TRUE))
    }
  }
  -stats::dnorm(z, td$zeta_prior_mean[s], td$zeta_prior_sd[s], log = TRUE)
}

# log of what state s has at time t after starting from log abundance `z0`,
# with net rate ns: its own growth plus the influx from its neighbours, which
# start from `from0` (their log abundances; only the edges `edges` count).
sparse_joint_own_log_n <- function(st, s, z0, ns, t, from0, edges = st$into[[s]]) {
  own <- z0 + ns * t
  if (!length(edges)) return(own)
  i <- st$from[edges]
  sparse_joint_own_lse(c(own, log(2) + st$td$edge_log_prob[edges] + st$par$g[i] + from0[i] +
                            sparse_joint_own_log_phi1(ns, st$net[i], t)))
}

# The own objective of state s: a function of c(g, its latent abundances),
# for the counts `ys` (one per time point), with the starting values and
# bounds of the search.
sparse_joint_own_objective <- function(st, s, ys) {
  td <- st$td
  par <- st$par
  w <- st$used[s, ]
  ref <- td$reference + 1L
  if (st$model == "sparse_joint") {
    centre <- st$net[ref] * td$dt
    fn <- function(p) {
      ns <- st$self[s] * exp(p[1L])
      l1 <- sparse_joint_own_log_n(st, s, p[2L], ns, td$dt, par$zeta) - centre + par$theta1
      sparse_joint_own_init(st, s, p[2L], ns) - w[1L] * sparse_joint_own_pois(ys[1L], p[2L]) -
        w[2L] * sparse_joint_own_pois(ys[2L], l1)
    }
    start <- c(par$g[s], par$zeta[s])
    nz <- 1L
  } else if (st$model == "sparse_joint_series") {
    tau <- td$tau
    centre <- st$net[ref] * tau
    fn <- function(p) {
      ns <- st$self[s] * exp(p[1L])
      out <- sparse_joint_own_init(st, s, p[2L], ns) - w[1L] * sparse_joint_own_pois(ys[1L], p[2L])
      for (k in seq_along(tau)) {
        if (w[k + 1L] == 0) next
        l <- sparse_joint_own_log_n(st, s, p[2L], ns, tau[k], par$zeta) - centre[k] + par$theta[k]
        out <- out - w[k + 1L] * sparse_joint_own_pois(ys[k + 1L], l)
      }
      out
    }
    start <- c(par$g[s], par$zeta[s])
    nz <- 1L
  } else if (st$model == "sparse_joint_pairwise") {
    pairs <- which(td$member[s, ] > 0)
    zeta <- as.matrix(par$zeta)
    centre <- st$net[ref] * td$interval
    fn <- function(p) {
      ns <- st$self[s] * exp(p[1L])
      out <- 0
      for (m in seq_along(pairs)) {
        k <- pairs[m]
        z <- p[1L + m]
        e <- st$into[[s]]
        e <- e[td$member[st$from[e], k] > 0]
        init <- NULL
        if (td$balance == 1L && td$pair_major[s, k] == 0) {
          fe <- e[td$pair_major[st$from[e], k] > 0]
          if (length(fe)) {
            q <- st$from[fe]
            start_k <- sparse_joint_own_lse(log(2) + td$edge_log_prob[fe] + par$g[q] + zeta[q, k] -
                                              st$net[q] * td$pair_history[k] +
                                              sparse_joint_own_log_phi1(ns, st$net[q], td$pair_history[k]))
            init <- -stats::dnorm(z - start_k - td$balance_offset[s, k], 0, td$balance_sd, log = TRUE)
          }
        }
        if (is.null(init)) init <- -stats::dnorm(z, td$zeta_prior_mean[s, k], td$zeta_prior_sd[s, k], log = TRUE)
        l <- sparse_joint_own_log_n(st, s, z, ns, td$interval[k], zeta[, k], e) - centre[k] + par$theta[k]
        out <- out + init - w[1L] * (sparse_joint_own_pois(ys[k], z) + sparse_joint_own_pois(ys[k + 1L], l))
      }
      out
    }
    start <- c(par$g[s], zeta[s, pairs])
    nz <- length(pairs)
  } else {
    later <- as.matrix(par$zeta_later)
    interval <- td$interval
    centre <- st$net[ref] * interval
    drift <- exp(par$log_scale[length(par$log_scale)])
    v <- drift * interval * td$drift_weight[s, ]
    v <- v * td$max_drift_var / (v + td$max_drift_var)
    fn <- function(p) {
      ns <- st$self[s] * exp(p[1L])
      out <- sparse_joint_own_init(st, s, p[2L], ns) - w[1L] * sparse_joint_own_pois(ys[1L], p[2L])
      for (k in seq_along(interval)) {
        from_s <- if (k == 1L) p[2L] else p[1L + k]
        from0 <- if (k == 1L) par$zeta else later[, k - 1L]
        mean <- sparse_joint_own_log_n(st, s, from_s, ns, interval[k], from0) - centre[k] + par$theta[k]
        zk <- p[2L + k]
        out <- out - stats::dnorm(zk - mean - td$transition_offset[s, k], 0, sqrt(v[k]), log = TRUE) -
          w[1L] * sparse_joint_own_pois(ys[k + 1L], zk)
      }
      out
    }
    start <- c(par$g[s], par$zeta[s], later[s, ])
    nz <- 1L + ncol(later)
  }
  list(fn = fn, start = start, lower = c(st$g_bounds[1L], rep(st$z_bounds[1L], nz)),
       upper = c(st$g_bounds[2L], rep(st$z_bounds[2L], nz)))
}

# The log intensities of state s at the fit's mode, as its own objective
# computes them: the check that the objective here reproduces the fit's
# (sparse_joint_own_check()).
sparse_joint_own_intensity <- function(st, s) {
  td <- st$td
  par <- st$par
  ns <- st$net[s]
  ref <- td$reference + 1L
  if (st$model == "sparse_joint") {
    return(sparse_joint_own_log_n(st, s, par$zeta[s], ns, td$dt, par$zeta) - st$net[ref] * td$dt + par$theta1)
  }
  if (st$model == "sparse_joint_series") {
    return(vapply(seq_along(td$tau), function(k) {
      sparse_joint_own_log_n(st, s, par$zeta[s], ns, td$tau[k], par$zeta) - st$net[ref] * td$tau[k] + par$theta[k]
    }, numeric(1)))
  }
  if (st$model == "sparse_joint_pairwise") {
    zeta <- as.matrix(par$zeta)
    return(vapply(seq_len(ncol(zeta)), function(k) {
      if (td$member[s, k] == 0) return(NA_real_)
      e <- st$into[[s]]
      e <- e[td$member[st$from[e], k] > 0]
      sparse_joint_own_log_n(st, s, zeta[s, k], ns, td$interval[k], zeta[, k], e) -
        st$net[ref[k]] * td$interval[k] + par$theta[k]
    }, numeric(1)))
  }
  later <- as.matrix(par$zeta_later)
  vapply(seq_along(td$interval), function(k) {
    from0 <- if (k == 1L) par$zeta else later[, k - 1L]
    mean <- sparse_joint_own_log_n(st, s, from0[s], ns, td$interval[k], from0) - st$net[ref] * td$interval[k] +
      par$theta[k]
    later[s, k] - mean
  }, numeric(1))
}

# Largest difference between the log intensities (state space: transition
# residuals) of the own objectives at the fit's mode and those the fit's TMB
# objective reports, over every state of the domain.
sparse_joint_own_check <- function(fit) {
  st <- sparse_joint_own_setup(fit)
  rb <- sparse_joint_rebuild(fit)
  rep <- rb$objective$report(rb$objective$env$par)
  tmb <- switch(st$model,
                sparse_joint = matrix(rep$log_intensity1, ncol = 1L),
                sparse_joint_state_space = as.matrix(rep$transition_residual),
                as.matrix(rep$log_intensity))
  own <- t(vapply(seq_len(st$n), function(s) sparse_joint_own_intensity(st, s), numeric(ncol(tmb))))
  if (ncol(tmb) == 1L) own <- matrix(own, ncol = 1L)
  if (st$model == "sparse_joint_pairwise") tmb[st$td$member == 0] <- NA
  max(abs(own - tmb), na.rm = TRUE)
}

# Net growth rates (2 P_ii - 1) exp(g) of `karyotypes`, each from its own
# counts (see the top of this file). NA for a karyotype outside the domain or
# without cells in the counts the fit uses. Attributes: `at_bound` (the
# search ended at a bound of the birth rate) and `cells` (cells of each
# karyotype in those counts, before any resampling).
sparse_joint_own_fitness <- function(fit, karyotypes, resample = TRUE, seed = 1L) {
  st <- sparse_joint_own_setup(fit)
  state <- match(karyotypes, st$labels)
  cells <- vapply(seq_along(karyotypes), function(i) {
    if (is.na(state[i])) NA_real_ else sum(st$y[state[i], ] * (st$used[state[i], ] > 0))
  }, numeric(1))
  y <- st$y
  if (resample) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()),
            add = TRUE)
    set.seed(seed)
    for (t in seq_len(ncol(y))) y[, t] <- stats::rbinom(st$n, st$cells[t], pmin(1, st$y[, t] / st$cells[t]))
  }
  net <- rep(NA_real_, length(karyotypes))
  at_bound <- rep(NA, length(karyotypes))
  for (i in which(!is.na(state) & cells > 0)) {
    s <- state[i]
    obj <- sparse_joint_own_objective(st, s, y[s, ])
    start <- pmin(pmax(obj$start, obj$lower), obj$upper)
    # Far from the data (near a bound) an intensity can approach overflow;
    # the search sees at most 1e100 there, so that its finite differences
    # stay finite.
    fn <- function(p) {
      v <- obj$fn(p)
      if (is.finite(v) && v < 1e100) v else 1e100
    }
    best <- stats::optim(start, fn, method = "L-BFGS-B", lower = obj$lower, upper = obj$upper,
                         control = list(maxit = 500L, factr = 1e5))
    net[i] <- st$self[s] * exp(best$par[1L])
    at_bound[i] <- best$par[1L] <= obj$lower[1L] + 1e-6 || best$par[1L] >= obj$upper[1L] - 1e-6
  }
  structure(stats::setNames(net, karyotypes), at_bound = at_bound, cells = cells)
}
