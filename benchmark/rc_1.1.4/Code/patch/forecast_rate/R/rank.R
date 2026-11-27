# Ranking rules for the observed karyotypes left out of a sparse joint fit's
# domain (PLAN.md). Wraps karst's predict(); the package is not changed. `k`
# is karst's namespace or the local stand-in (od_karst() in
# ../outside_domain/R/common.R).

# L1 (lattice step) distance between the rows of two copy-number matrices.
or_l1 <- function(a, b) {
  d <- matrix(0, nrow(a), nrow(b))
  for (c in seq_len(ncol(a))) d <- d + abs(outer(a[, c], b[, c], "-"))
  d
}

or_counts <- function(fit) {
  x <- as.matrix(fit$data$counts)
  if (is.null(rownames(x))) rownames(x) <- fit$data$labels
  storage.mode(x) <- "double"
  x
}

# Time between the first and the last time point (days).
or_span <- function(fit) {
  if (!is.null(fit$data$times)) return(diff(range(fit$data$times)))
  sum(fit$data$dt)
}

# Features of each karyotype in O = fit$dropped (PLAN.md, Definitions), and the
# fit-level quantities the rules need. NULL when the fit has no O or no
# observed level (then predict() already uses the trend for them).
or_features <- function(k, fit) {
  O <- fit$dropped
  trend <- fit$trend
  if (!length(O) || is.null(trend$observed_level)) return(NULL)
  tr_of <- function(cn) as.numeric(k$sparse_joint_design(cn, trend$centre, trend$major, trend$shell_power)$design %*%
                                     trend$coefficients)
  cnO <- od_parse(O)
  s <- fit$summary
  total <- rowSums(as.matrix(s[grep("^count_", names(s))]))
  J <- which(total > 0)
  cnJ <- od_parse(s$karyotype[J])
  rJ <- s$log_birth_mean[J] - tr_of(cnJ)
  # residual correlation at lattice distance 1 and 2 among observed in-domain karyotypes
  dJJ <- or_l1(cnJ, cnJ)
  rho <- c(0, 0); rho_raw <- c(NA_real_, NA_real_); pairs <- c(0L, 0L)
  for (d in 1:2) {
    ij <- which(dJJ == d & upper.tri(dJJ), arr.ind = TRUE)
    pairs[d] <- nrow(ij)
    if (nrow(ij) >= 10L) {
      r <- suppressWarnings(stats::cor(c(rJ[ij[, 1]], rJ[ij[, 2]]), c(rJ[ij[, 2]], rJ[ij[, 1]])))
      rho_raw[d] <- r
      rho[d] <- if (is.finite(r)) min(max(r, 0), 1) else 0
    }
  }
  dOJ <- or_l1(cnO, cnJ)
  dstar <- apply(dOJ, 1L, min)
  nn <- vapply(seq_along(O), function(i) {
    d <- dstar[i]
    if (d > 2) return(0)
    rho[d] * mean(rJ[dOJ[i, ] == d])
  }, numeric(1))
  # count features
  x <- or_counts(fit)
  depth <- colSums(x)
  iO <- match(O, rownames(x))
  if (anyNA(iO)) stop("dropped karyotypes missing from the counts")
  y <- x[iO, , drop = FALSE]
  nb <- or_l1(cnO, od_parse(rownames(x))) == 1
  lam <- nb %*% x
  a <- log(rowSums(sweep(y + 0.5, 2L, depth, "/"))) - log(rowSums(sweep(lam + 0.5, 2L, depth, "/")))
  last <- ncol(x)
  sfe <- (log((y[, last] + 0.5) / depth[last]) - log((y[, 1L] + 0.5) / depth[1L])) / or_span(fit)
  minor <- !s$major & total > 0
  sigma_m <- if (sum(minor) >= 2L) stats::sd(s$log_birth_mean[minor]) else trend$observed_level[["sd"]]
  self <- 2 * k$sparse_joint_self_probability(cnO, fit$data$beta) - 1
  list(info = list(L = trend$observed_level[["log_birth"]], rho = rho, rho_raw = rho_raw, pairs = pairs, sigma_m = sigma_m,
                   ref = trend$reference$fitness, n_O = length(O), n_J = length(J)),
       O = data.frame(karyotype = O, self = self, tr = tr_of(cnO), dstar = dstar, n = nn, a = as.numeric(a),
                      s = as.numeric(sfe), y_first = y[, 1L], y_last = y[, last], stringsAsFactors = FALSE))
}

or_z <- function(x) {
  if (length(x) < 3L) return(0 * x)
  v <- sqrt(mean((x - mean(x))^2))
  if (!is.finite(v) || v == 0) 0 * x else (x - mean(x)) / v
}

# Log division rate of each karyotype in O under a rule. `w`: c(a = w_a, s = w_s) for R3.
or_lb <- function(feat, rule = c("current", "R1", "R2", "R3"), w = c(a = 0, s = 0)) {
  rule <- match.arg(rule)
  # Refinement of 2026-11-25 (PLAN.md, before any result): every rule is
  # L + sigma_m * z(score). Centring the raw trend on its mean put near
  # karyotypes far above the level when O also held far ones (local Salehi
  # fit: +7.7 per day). The ordering within O is the same as with the raw sums.
  o <- feat$O
  L <- feat$info$L
  if (rule == "current") return(rep(L, nrow(o)))
  score <- switch(rule, R1 = o$tr, R2 = o$tr + o$n,
                  R3 = or_z(o$tr + o$n) + w[["a"]] * or_z(o$a) + w[["s"]] * or_z(o$s))
  L + feat$info$sigma_m * or_z(score)
}

# Relative fitness of each karyotype in O under a rule (as predict() computes it).
or_relative <- function(feat, rule, w = c(a = 0, s = 0)) feat$O$self * exp(or_lb(feat, rule, w)) - feat$info$ref

# predict()'s table `base` for `labels`, with the rows in O replaced under a
# rule: fitness_mean, relative_fitness and the division rate (fitness / self).
# sd columns are unchanged (PLAN.md, D3a).
or_apply <- function(base, feat, rule, w = c(a = 0, s = 0)) {
  if (is.null(feat) || rule == "current") return(base)
  i <- match(base$karyotype, feat$O$karyotype)
  hit <- which(!is.na(i))
  lb <- or_lb(feat, rule, w)
  base$fitness_mean[hit] <- feat$O$self[i[hit]] * exp(lb[i[hit]])
  base$relative_fitness[hit] <- base$fitness_mean[hit] - feat$info$ref
  base
}

# Local runs only (no installed karst): an R version of karst's compiled
# chromosome transition probability (src/karst_core.h), added to the stand-in
# environment from od_karst(use_package = FALSE).
or_local_karst <- function() {
  k <- od_karst(use_package = FALSE)
  k$karst_pij_cpp <- function(parent_cn, child_cn, beta) {
    diff <- abs(parent_cn - child_cn); if (parent_cn < 0 || child_cn < 0 || diff > parent_cn) return(0)
    total <- 0
    for (m in seq(diff, parent_cn, by = 2)) {
      lost <- (m + parent_cn - child_cn) / 2
      if (lost < 0 || lost > m) next
      total <- total + choose(parent_cn, m) * beta^m * (1 - beta)^(parent_cn - m) * 0.5^m * choose(m, lost)
    }
    total
  }
  k
}
