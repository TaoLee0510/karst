# Copy of ../outside_domain/R/engine.R (unchanged), so this study is self-contained.
# Deterministic expected-value version of the benchmark ABM (alfakR dev2
# f93b2d2, src/abm_core.cpp, run_karyotype_abm) on a fixed karyotype list.
#
# The ABM, per step of length dt: a karyotype j with n_j cells and rate f_j > 0
# divides Poisson(n_j f_j dt) times (capped at n_j); rates <= 0 do not divide.
# A division draws Binomial(N_j, p) errors (N_j = total copies); each error
# moves one copy of a chromosome, chosen uniformly among the N_j copies, from
# one daughter to the other. Daughters with a negative copy number, or off the
# fitness map (the list), are dropped. There is no death; culling above the
# population limit keeps a uniform fraction, so frequencies are unchanged in
# expectation.
#
# Here, per step: the expected divisions d_j = n_j max(f_j, 0) dt; with no
# error (probability (1 - p)^N_j) a division gives two copies of j; with one
# error (N_j p (1 - p)^(N_j - 1)) on chromosome i (probability c_ij / N_j) it
# gives j + e_i and j - e_i; two or more errors (the rest) are dropped, as are
# daughters off the list. The cap at n_j and stochastic loss of rare types are
# not modelled. Frequencies are returned normalised.

# Index of each one-copy gain and loss of each karyotype within the list.
engine_build <- function(labels) {
  cn <- do.call(rbind, lapply(strsplit(labels, ".", fixed = TRUE), as.integer))
  storage.mode(cn) <- "integer"
  key <- function(m) do.call(paste, c(as.data.frame(m), sep = "."))
  plus <- minus <- matrix(NA_integer_, nrow(cn), ncol(cn))
  for (i in seq_len(ncol(cn))) {
    up <- cn; up[, i] <- up[, i] + 1L
    plus[, i] <- match(key(up), labels)
    dn <- cn; dn[, i] <- dn[, i] - 1L
    ok <- dn[, i] >= 0L
    minus[ok, i] <- match(key(dn[ok, , drop = FALSE]), labels)
  }
  list(labels = labels, cn = cn, copies = rowSums(cn), plus = plus, minus = minus)
}

# Propagate counts `n0` (aligned with the list) for `steps` steps of `dt` at
# division rates `rates` and missegregation probability `p` per copy and
# division. `record`: the steps whose frequencies are returned (columns).
engine_run <- function(eng, n0, rates, p, dt = 0.1, steps, record = steps) {
  stopifnot(length(n0) == length(eng$labels), length(rates) == length(eng$labels), all(is.finite(rates)))
  n <- as.numeric(n0)
  N <- eng$copies
  p0 <- (1 - p)^N
  p1 <- N * p * (1 - p)^pmax(N - 1, 0)
  share <- eng$cn / pmax(N, 1)  # c_ij / N_j
  r <- pmax(rates, 0) * dt
  out <- matrix(NA_real_, length(n), length(record))
  lost <- 0
  for (s in seq_len(steps)) {
    d <- n * r
    nn <- n - d + 2 * p0 * d
    e1 <- d * p1
    for (i in seq_len(ncol(eng$cn))) {
      w <- e1 * share[, i]
      for (idx in list(eng$plus[, i], eng$minus[, i])) {
        ok <- !is.na(idx) & w > 0
        nn[idx[ok]] <- nn[idx[ok]] + w[ok]
        lost <- lost + sum(w[is.na(idx) & w > 0])
      }
    }
    lost <- lost + sum(d * (1 - p0 - p1))
    n <- nn
    hit <- match(s, record)
    if (!is.na(hit)) out[, hit] <- n / sum(n)
  }
  rownames(out) <- eng$labels
  attr(out, "lost_cells") <- lost
  attr(out, "final_total") <- sum(n)
  out
}
