# P1-P3 of dev/realdata/noise/PLAN_v2.md: the neutral missegregation envelope
# (beta in [5e-5, 1e-2]), standardised excess z of each karyotype, the Moran's I
# of z over karyotypes one copy apart, and the verdict.

env_betas <- exp(seq(log(5e-5), log(1e-2), length.out = 13))

# One-division transition of one chromosome type, copy number i -> j, as
# karst's chromosome_transition_probability(): m ~ Binomial(i, beta) copies
# missegregate, each going to this daughter with probability 1/2.
env_kernel <- function(beta, cmax) {
  K <- matrix(0, cmax + 1L, cmax + 1L)
  for (i in 0:cmax) for (j in 0:cmax) {
    d <- abs(i - j)
    if (d > i) next
    tot <- 0
    for (m in seq(d, i, by = 2)) {
      lost <- (m + i - j) / 2
      if (lost < 0 || lost > m) next
      tot <- tot + choose(i, m) * beta^m * (1 - beta)^(i - m) * 0.5^m * choose(m, lost)
    }
    K[i + 1L, j + 1L] <- tot
  }
  K
}

# Sum_m Pois(m; lambda) K^m (divisions along a lineage over the interval).
env_transition <- function(K, lambda) {
  if (lambda <= 0) return(diag(nrow(K)))
  mmax <- max(20L, ceiling(lambda + 10 * sqrt(lambda) + 10))
  out <- matrix(0, nrow(K), ncol(K)); Km <- diag(nrow(K))
  w <- stats::dpois(0:mmax, lambda)
  for (m in 0:mmax) { out <- out + w[m + 1L] * Km; Km <- Km %*% K }
  out
}

env_parse <- function(labels) do.call(rbind, lapply(strsplit(labels, ".", fixed = TRUE), as.integer))

# x0, x1: named count vectors; dt: interval; rho: divisions per unit time.
env_analyse <- function(x0, x1, dt, rho, perms = 999L, seed = 1L) {
  x0 <- x0[x0 > 0]; x1 <- x1[x1 > 0]
  n0 <- sum(x0); n1 <- sum(x1)
  J <- names(x0); U <- union(J, names(x1))
  kJ <- env_parse(J); kU <- env_parse(U)
  cmax <- max(c(kJ, kU)) + 4L
  p0 <- as.numeric(x0) / n0
  lambda <- 2 * rho * dt
  y1 <- stats::setNames(numeric(length(U)), U); y1[names(x1)] <- x1
  qs <- matrix(NA_real_, length(U), length(env_betas)); vs <- qs
  for (b in seq_along(env_betas)) {
    E <- env_transition(env_kernel(env_betas[b], cmax), lambda)
    P <- matrix(1, length(J), length(U))
    for (ch in seq_len(ncol(kJ))) P <- P * E[cbind(rep(kJ[, ch] + 1L, times = length(U)), rep(kU[, ch] + 1L, each = length(J)))]
    q <- as.numeric(crossprod(p0, P))
    qs[, b] <- q
    vs[, b] <- pmax(as.numeric(crossprod(p0, P^2)) - q^2, 0) / n0
  }
  ib_max <- max.col(qs, ties.method = "first"); ib_min <- max.col(-qs, ties.method = "first")
  q_max <- qs[cbind(seq_along(U), ib_max)]; v_max <- vs[cbind(seq_along(U), ib_max)]
  q_min <- qs[cbind(seq_along(U), ib_min)]; v_min <- vs[cbind(seq_along(U), ib_min)]
  unseen <- !(U %in% J)
  q_max[unseen] <- q_max[unseen] + 3 / n0
  in0 <- U %in% J
  rise <- y1 >= 5 & y1 > n1 * q_max
  fall <- in0 & n1 * q_min >= 5 & y1 < n1 * q_min
  z <- numeric(length(U))
  z[rise] <- (y1[rise] - n1 * q_max[rise]) / sqrt(n1 * q_max[rise] * (1 - q_max[rise]) + n1^2 * v_max[rise])
  z[fall] <- (y1[fall] - n1 * q_min[fall]) / sqrt(n1 * q_min[fall] * (1 - q_min[fall]) + n1^2 * v_min[fall])
  abnormal <- abs(z) > 3
  excess <- (sum((y1 - n1 * q_max)[rise & abnormal]) + sum((n1 * q_min - y1)[fall & abnormal])) / n1
  testable <- y1 >= 5 | (in0 & n1 * q_min >= 5)
  T <- which(testable)
  # edges: pairs of testable karyotypes one copy apart
  kt <- kU[T, , drop = FALSE]
  edges <- if (length(T) >= 2) {
    d <- as.matrix(stats::dist(kt, method = "manhattan"))
    which(d == 1 & upper.tri(d), arr.ind = TRUE)
  } else matrix(integer(0), 0, 2)
  zt <- z[T]
  moran <- function(v) if (nrow(edges) == 0 || sum(v^2) == 0) NA_real_ else
    (length(v) / nrow(edges)) * sum(v[edges[, 1]] * v[edges[, 2]]) / sum(v^2)
  I_obs <- moran(zt); p <- NA_real_
  if (is.finite(I_obs)) {
    set.seed(seed)
    I_perm <- replicate(perms, moran(sample(zt)))
    p <- (1 + sum(I_perm >= I_obs)) / (perms + 1)
  }
  # report: abnormal karyotypes with a same-sign abnormal one within 2 steps
  ab <- which(abnormal)
  near_share <- NA_real_
  if (length(ab) >= 2) {
    da <- as.matrix(stats::dist(kU[ab, , drop = FALSE], method = "manhattan")); diag(da) <- Inf
    same <- outer(sign(z[ab]), sign(z[ab])) > 0
    near_share <- mean(apply(da <= 2 & same, 1, any))
  }
  n_ab <- sum(abnormal)
  verdict <- if (n_ab < 3 || excess < 0.05) "no_anomaly" else if (nrow(edges) < 10 || !is.finite(p)) "undetermined" else
    if (p < 0.05) "selection" else if (p >= 0.20 && n_ab >= 5) "sampling" else "undetermined"
  list(summary = data.frame(n0 = n0, n1 = n1, dt = dt, rho = rho, lambda = lambda, n_karyotypes = length(U),
                            n_testable = length(T), n_edges = nrow(edges), n_abnormal = n_ab,
                            n_rise = sum(rise & abnormal), n_fall = sum(fall & abnormal), excess = excess,
                            moran = I_obs, p = p, near_share = near_share, verdict = verdict),
       karyotypes = data.frame(karyotype = U, x0 = ifelse(in0, x0[U], 0), x1 = y1, q_max = q_max, q_min = q_min,
                               z = z, abnormal = abnormal, testable = testable, row.names = NULL))
}

# Noise injection of P4: independent overdispersion (Dirichlet-multinomial,
# variance inflation phi) or region sampling (families: nearest of the 5 most
# frequent karyotypes of both time points; one log-normal factor per family).
env_inject <- function(x, kind, level, families = NULL) {
  n <- sum(x); p <- as.numeric(x) / n
  if (kind == "dm") {
    alpha <- (n - 1) / (level - 1) - 1
    g <- stats::rgamma(length(p), shape = pmax(alpha * p, 1e-12))
    pp <- g / sum(g)
  } else if (kind == "region") {
    fam <- families[names(x)]
    f <- exp(stats::rnorm(length(unique(fam)), 0, level)); names(f) <- unique(fam)
    pp <- p * f[fam]; pp <- pp / sum(pp)
  } else stop("unknown injection")
  y <- as.numeric(stats::rmultinom(1, n, pp)); names(y) <- names(x)
  y[y > 0]
}

env_families <- function(x0, x1) {
  tot <- tapply(c(x0, x1), names(c(x0, x1)), sum)
  top <- names(sort(tot, decreasing = TRUE))[seq_len(min(5L, length(tot)))]
  all <- names(tot)
  ka <- env_parse(all); kt <- env_parse(top)
  near <- vapply(seq_len(nrow(ka)), function(i) which.min(rowSums(abs(sweep(kt, 2, ka[i, ])))), integer(1))
  stats::setNames(top[near], all)
}
