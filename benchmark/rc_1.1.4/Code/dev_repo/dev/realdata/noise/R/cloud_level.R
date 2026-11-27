# v6 of dev/realdata/noise/PLAN_v6.md: the level of the one-step neighbour
# cloud. S_A: cloud cells around the high-frequency karyotypes per cell of
# them; S_B: cell pairs one copy apart per pair of identical cells; both per
# route of missegregation. Needs envelope.R (env_analyse, env_parse).

cl_min_cells <- 10
cl_beta_default <- 5e-5
cl_cap_multiple <- 10

# Ordered pairs (i, j) of rows of K one copy apart.
cl_adjacent_pairs <- function(K) {
  key <- do.call(paste, c(as.data.frame(K), sep = "."))
  idx <- stats::setNames(seq_len(nrow(K)), key)
  out <- vector("list", 2L * ncol(K))
  for (ch in seq_len(ncol(K))) for (s in c(-1L, 1L)) {
    Kn <- K; Kn[, ch] <- Kn[, ch] + s
    j <- idx[do.call(paste, c(as.data.frame(Kn), sep = "."))]
    ok <- !is.na(j)
    out[[2L * (ch - 1L) + (s > 0) + 1L]] <- cbind(which(ok), as.integer(j[ok]))
  }
  do.call(rbind, out)
}

cl_routes <- function(K) rowSums(K * (1 + (K >= 2)))

# Both statistics of one dataset; x0, x1 named count vectors.
cl_stats <- function(x0, x1, dt, rho, seed) {
  env <- env_analyse(x0, x1, dt, rho, perms = 1L, seed = seed)
  k <- env$karyotypes
  a <- k$x0; b <- k$x1
  K <- env_parse(k$karyotype)
  G <- cl_routes(K)
  pr <- cl_adjacent_pairs(K)
  high <- a >= cl_min_cells | b >= cl_min_cells
  near_high <- rep(FALSE, length(a))
  if (nrow(pr)) near_high[unique(pr[high[pr[, 2]], 1])] <- TRUE
  cloud <- !high & near_high & !(k$abnormal & k$z > 0)
  H <- sum(a[high] + b[high]); C <- sum(a[cloud] + b[cloud])
  GA <- if (H > 0) sum((a + b)[high] * G[high]) / H else NA_real_
  S_A <- if (H > 0) (C + 0.5) / H / GA else NA_real_
  N0 <- if (nrow(pr)) sum(a[pr[, 1]] * a[pr[, 2]]) else 0
  N1 <- if (nrow(pr)) sum(b[pr[, 1]] * b[pr[, 2]]) else 0
  I0 <- sum(a * (a - 1)); I1 <- sum(b * (b - 1))
  GB <- sum((a + b) * G) / sum(a + b)
  S_B <- (N0 + N1 + 0.5) / (I0 + I1 + 0.5) / GB
  data.frame(n0 = sum(a), n1 = sum(b), n_karyotypes = length(a), n_high = sum(high), H = H, C = C,
             cloud_k = sum(cloud), GA = GA, S_A = S_A, N0 = N0, N1 = N1, I0 = I0, I1 = I1, GB = GB, S_B = S_B)
}

# beta_used from a statistic and the calibration (c, d, sigma) of step 2.
cl_beta_used <- function(S, cal) {
  x <- (log2(S) - cal$c - stats::qnorm(0.75) * cal$sigma) / cal$d
  x[!is.finite(x) & !is.na(x) & x < 0] <- 0
  out <- cl_beta_default * 2^pmin(pmax(x, 0), log2(cl_cap_multiple))
  out[is.na(S)] <- cl_beta_default
  out
}
