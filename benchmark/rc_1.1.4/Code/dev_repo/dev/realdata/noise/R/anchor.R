# v5 of dev/realdata/noise/PLAN_v5.md: high-frequency karyotypes, the
# reference karyotypes and their one-step neighbour cloud, and the lower bound
# on the missegregation rate from the cloud's growth relative to the
# references. Needs envelope.R (env_analyse, env_parse).

anc_min_cells <- 10
anc_beta_default <- 5e-5
anc_beta_cap <- 5e-4

# Per karyotype of env_analyse(): high frequency, the two-proportion z without
# missegregation and the "falling" flag.
anc_karyotypes <- function(env) {
  k <- env$karyotypes
  a <- k$x0; b <- k$x1; n0 <- sum(a); n1 <- sum(b)
  pb <- (a + b) / (n0 + n1)
  k$z2 <- (b / n1 - a / n0) / sqrt(pb * (1 - pb) * (1 / n0 + 1 / n1))
  k$high <- a >= anc_min_cells | b >= anc_min_cells
  k$falling <- k$high & k$z2 < -1.96
  k
}

# L1 distance of every karyotype to the nearest one of rows idx.
anc_distance <- function(K, idx) {
  out <- rep(Inf, nrow(K))
  for (i in idx) out <- pmin(out, rowSums(abs(sweep(K, 2, K[i, ]))))
  out
}

# The bound of one dataset for one version ("a" or "b"); k from anc_karyotypes().
anc_bound <- function(k, dt, rho, version, seed, draws = 4000L) {
  K <- env_parse(k$karyotype)
  a <- k$x0; b <- k$x1
  ref <- which(k$high & !k$falling)
  other <- which(k$high & k$falling)
  empty <- data.frame(version = version, n_high = sum(k$high), n_ref = length(ref), n_falling = length(other),
                      cloud_k = 0L, C0 = 0, C1 = 0, A0 = sum(a[ref]), A1 = sum(b[ref]), Gbar = NA_real_,
                      r0 = NA_real_, r1 = NA_real_, point = NA_real_, L = NA_real_, triggered = FALSE,
                      beta_used = anc_beta_default)
  if (!length(ref)) return(empty)
  d_ref <- anc_distance(K, ref)
  d_other <- anc_distance(K, other)
  cloud <- !k$high & d_ref == 1 & d_other != 1
  cloud <- cloud & !(k$abnormal & k$z > 0)
  if (version == "b" && sum(cloud) > 0) {
    idx <- which(cloud)
    drop <- idx[order(-(a[idx] + b[idx]))][seq_len(min(2L, length(idx)))]
    cloud[drop] <- FALSE
  }
  C0 <- sum(a[cloud]); C1 <- sum(b[cloud]); A0 <- sum(a[ref]); A1 <- sum(b[ref])
  cn <- K[ref, , drop = FALSE]
  G <- rowSums(cn * (1 + (cn >= 2)))
  w <- a[ref] / sum(a) + b[ref] / sum(b)
  Gbar <- sum(w * G) / sum(w)
  denom <- rho * Gbar * dt
  out <- empty
  out$cloud_k <- sum(cloud); out$C0 <- C0; out$C1 <- C1; out$Gbar <- Gbar
  out$r0 <- if (A0 > 0) C0 / A0 else Inf
  out$r1 <- if (A1 > 0) C1 / A1 else Inf
  out$point <- if (is.finite(out$r0) && is.finite(out$r1)) (out$r1 - out$r0) / denom else NA_real_
  if (C0 + C1 == 0) return(out)
  set.seed(seed)
  p0 <- stats::rbeta(draws, C0 + 0.5, A0 + 0.5); p1 <- stats::rbeta(draws, C1 + 0.5, A1 + 0.5)
  out$L <- as.numeric(stats::quantile((p1 / (1 - p1) - p0 / (1 - p0)) / denom, 0.05))
  out$triggered <- out$L > anc_beta_default
  out$beta_used <- if (out$triggered) min(anc_beta_cap, out$L) else anc_beta_default
  out
}

anc_seed <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L

# Everything of one dataset: the envelope (ceiling 1e-2), the karyotype table
# of the high-frequency ones and the two bounds.
anc_dataset <- function(x0, x1, dt, rho, id) {
  env <- env_analyse(x0, x1, dt, rho, perms = 1L, seed = anc_seed(id))
  k <- anc_karyotypes(env)
  list(n0 = sum(k$x0), n1 = sum(k$x1), n_karyotypes = nrow(k), dt = dt, rho = rho,
       high = k[k$high, c("karyotype", "x0", "x1", "z", "abnormal", "testable", "z2", "falling")],
       bounds = rbind(anc_bound(k, dt, rho, "a", anc_seed(id)), anc_bound(k, dt, rho, "b", anc_seed(id))),
       k_all = k)
}
