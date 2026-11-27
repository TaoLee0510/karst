# Measures of N1 of dev/realdata/noise/PLAN.md (and its corrections of
# 2026-11-20 11:40): from two samples of one population, the Good-Turing
# unseen share, the novelty, the apparent diversification rate r_hat and the
# temporal-method drift lower bound N_e.

noise_beta_hi <- 1e-2
noise_beta_default <- 5e-5

# x0, x1: named count vectors (karyotype labels "a.b.c"); dt: time between the
# samples; rho: divisions per unit time; passages: number of passages between
# the samples (NA if unknown).
noise_measure <- function(x0, x1, dt, rho, passages = NA_real_) {
  x0 <- x0[x0 > 0]; x1 <- x1[x1 > 0]
  n0 <- sum(x0); n1 <- sum(x1)
  u <- min(0.5, sum(x0 == 1) / n0)
  new <- !(names(x1) %in% names(x0))
  n_new <- sum(x1[new]) / n1
  k0 <- do.call(rbind, lapply(strsplit(names(x0), ".", fixed = TRUE), as.numeric))
  c_bar <- sum(rowSums(k0) * x0) / n0
  d_new <- if (any(new)) {
    k1 <- do.call(rbind, lapply(strsplit(names(x1)[new], ".", fixed = TRUE), as.numeric))
    vapply(seq_len(nrow(k1)), function(i) min(rowSums(abs(sweep(k0, 2, k1[i, ])))), numeric(1))
  } else numeric(0)
  share_d <- function(lo, hi) sum(x1[new][d_new >= lo & d_new <= hi]) / n1
  divisions <- dt * rho
  r_hat <- if (n_new <= u) 0 else if (n_new >= 1) Inf else -log((1 - n_new) / (1 - u)) / (c_bar * divisions)
  # temporal method (Nei-Tajima F_c) on the t0 karyotypes: those with >= 5 cells
  # each a class, the rest of t0 pooled; t1 cells outside t0 left out
  big <- names(x0)[x0 >= 5]
  cls0 <- c(x0[big], other = sum(x0[!(names(x0) %in% big)]))
  in0 <- x1[names(x1) %in% names(x0)]
  cls1 <- c(in0[big], other = sum(in0[!(names(in0) %in% big)]))
  cls1[is.na(cls1)] <- 0
  keep <- cls0 > 0 | cls1 > 0
  s0 <- sum(cls0); s1 <- sum(cls1)
  fc <- NA_real_; ne <- NA_real_
  if (sum(keep) >= 2 && s1 > 0) {
    xf <- cls0[keep] / s0; yf <- cls1[keep] / s1
    den <- (xf + yf) / 2 - xf * yf
    fc <- sum(((xf - yf)^2 / den)[den > 0]) / (sum(keep) - 1)
    excess <- fc - 1 / s0 - 1 / s1
    ne <- if (is.finite(passages) && passages > 0) (if (excess > 0) passages / excess else Inf) else NA_real_
  }
  data.frame(n0 = n0, n1 = n1, u = u, n_new = n_new, new_d1 = share_d(1, 1), new_d2 = share_d(2, 2),
             new_d3plus = share_d(3, Inf), c_bar = c_bar, dt = dt, rho = rho, divisions = divisions,
             r_hat = r_hat, multiple = r_hat / noise_beta_default, flagged = r_hat > noise_beta_hi,
             nu_hat = max(0, r_hat - noise_beta_default), fc = fc, s0 = s0, s1 = s1, passages = passages,
             ne = ne, classes = sum(keep))
}
