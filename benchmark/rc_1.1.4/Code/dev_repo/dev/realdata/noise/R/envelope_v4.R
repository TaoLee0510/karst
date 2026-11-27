# v4 of dev/realdata/noise (PLAN_v4.md, step 1): anomalies at each rate of the
# envelope grid separately, and the four classes with beta_hat and phi_hat.
# Sourced after R/envelope.R, which is unchanged.

# Anomalies under one rate at a time: the same expected frequencies, variance
# terms, unseen allowance (3 / n0), rise/fall rules, |z| > 3 and excess as
# env_analyse(), but with q(beta) instead of the maximum or minimum over rates.
env_by_beta <- function(x0, x1, dt, rho) {
  x0 <- x0[x0 > 0]; x1 <- x1[x1 > 0]
  n0 <- sum(x0); n1 <- sum(x1)
  J <- names(x0); U <- union(J, names(x1))
  kJ <- env_parse(J); kU <- env_parse(U)
  cmax <- max(c(kJ, kU)) + 4L
  p0 <- as.numeric(x0) / n0
  lambda <- 2 * rho * dt
  y1 <- stats::setNames(numeric(length(U)), U); y1[names(x1)] <- x1
  in0 <- U %in% J
  unseen <- !in0
  do.call(rbind, lapply(seq_along(env_betas), function(b) {
    E <- env_transition(env_kernel(env_betas[b], cmax), lambda)
    P <- matrix(1, length(J), length(U))
    for (ch in seq_len(ncol(kJ))) P <- P * E[cbind(rep(kJ[, ch] + 1L, times = length(U)), rep(kU[, ch] + 1L, each = length(J)))]
    q <- as.numeric(crossprod(p0, P))
    v <- pmax(as.numeric(crossprod(p0, P^2)) - q^2, 0) / n0
    q_hi <- q; q_hi[unseen] <- q_hi[unseen] + 3 / n0
    rise <- y1 >= 5 & y1 > n1 * q_hi
    fall <- in0 & n1 * q >= 5 & y1 < n1 * q
    z <- numeric(length(U))
    z[rise] <- (y1[rise] - n1 * q_hi[rise]) / sqrt(n1 * q_hi[rise] * (1 - q_hi[rise]) + n1^2 * v[rise])
    z[fall] <- (y1[fall] - n1 * q[fall]) / sqrt(n1 * q[fall] * (1 - q[fall]) + n1^2 * v[fall])
    ab <- abs(z) > 3
    excess <- (sum((y1 - n1 * q_hi)[rise & ab]) + sum((n1 * q - y1)[fall & ab])) / n1
    data.frame(beta = env_betas[b], n_abnormal = sum(ab), excess = excess, anomalous = sum(ab) >= 3 && excess >= 0.05)
  }))
}

# The class of PLAN_v4.md step 1 from env_analyse() and env_by_beta().
env_v4_class <- function(env, by_beta) {
  verdict <- env$summary$verdict
  fewest <- by_beta$beta[which.min(by_beta$n_abnormal)]   # which.min: first, i.e. the smaller rate on ties
  k <- env$karyotypes
  phi_hat <- max(1, if (any(k$testable)) mean(k$z[k$testable]^2) else 1)
  if (verdict == "selection") {
    list(class = "selection", beta_hat = NA_real_, phi_hat = phi_hat)
  } else if (verdict %in% c("sampling", "undetermined")) {
    list(class = "noise", beta_hat = fewest, phi_hat = phi_hat)
  } else if (isTRUE(by_beta$anomalous[1])) {
    ok <- by_beta$beta[!by_beta$anomalous]
    list(class = "higher_beta", beta_hat = if (length(ok)) min(ok) else fewest, phi_hat = phi_hat)
  } else list(class = "normal", beta_hat = NA_real_, phi_hat = phi_hat)
}

# Development and confirmation halves (PLAN_v4.md, fixed before any fit).
v4_half <- function(landscape_id) {
  dev <- c(sprintf("nz_L%02d", c(1, 2, 5, 7, 8, 10)), sprintf("nz_c2_L%02d", c(1, 2, 3, 7, 8, 9)))
  ifelse(landscape_id %in% dev, "dev", "confirm")
}

# Robust variance inflation (PLAN_v4.md correction 16:24, candidate T2r): at
# the single rate beta, the full signed z of every testable karyotype (y1 >= 5,
# or seen at the first time point with n1 q >= 5); unseen karyotypes get the
# 3 / n0 allowance; phi_r = max(1, median(z^2) / qchisq(0.5, 1)); 1 with fewer
# than 5 testable karyotypes.
env_phi_robust <- function(x0, x1, dt, rho, beta) {
  x0 <- x0[x0 > 0]; x1 <- x1[x1 > 0]
  n0 <- sum(x0); n1 <- sum(x1)
  J <- names(x0); U <- union(J, names(x1))
  kJ <- env_parse(J); kU <- env_parse(U)
  cmax <- max(c(kJ, kU)) + 4L
  p0 <- as.numeric(x0) / n0
  E <- env_transition(env_kernel(beta, cmax), 2 * rho * dt)
  P <- matrix(1, length(J), length(U))
  for (ch in seq_len(ncol(kJ))) P <- P * E[cbind(rep(kJ[, ch] + 1L, times = length(U)), rep(kU[, ch] + 1L, each = length(J)))]
  q <- as.numeric(crossprod(p0, P))
  v <- pmax(as.numeric(crossprod(p0, P^2)) - q^2, 0) / n0
  y1 <- stats::setNames(numeric(length(U)), U); y1[names(x1)] <- x1
  in0 <- U %in% J
  qt <- q; qt[!in0] <- qt[!in0] + 3 / n0
  testable <- y1 >= 5 | (in0 & n1 * q >= 5)
  z <- (y1 - n1 * qt) / sqrt(n1 * qt * (1 - qt) + n1^2 * v)
  z2 <- z[testable]^2
  list(phi_r = if (sum(testable) < 5) 1 else max(1, stats::median(z2) / stats::qchisq(0.5, 1)),
       n_testable = sum(testable), median_z2 = if (length(z2)) stats::median(z2) else NA_real_)
}
