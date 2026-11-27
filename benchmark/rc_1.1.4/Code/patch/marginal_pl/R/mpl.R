# Per-chromosome marginal pseudo-likelihood (PLAN.md). After a KARST fit, the
# trend's dosage and curvature coefficients are updated so that the observed
# karyotypes of the first time point, grown for dt at their fitness, give the
# copy-number distribution of every chromosome at the last time point:
#   PL(delta) = sum_i sum_c y1[i, c] log p1[i, c](delta)
#   maximise  -1/2 delta' S^-1 delta + w PL(delta)
# where delta = theta - theta_hat (intercept fixed: a common shift of fitness
# leaves p1 unchanged) and S is the fit's covariance of the coefficients.
# D1a: only karyotypes in the domain move (log birth + Dtrend); karyotypes of O
# keep KARST's common level; unseen karyotypes follow the new trend.
# D3a: selection only (no missegregation in the propagation).
# `k` is karst's namespace (od_karst()).

mpl_eps <- 1e-3  # mixing with a uniform over the copy states seen at t1

# Everything the PL needs from a fit (its own counts: first and last column).
mpl_setup <- function(k, fit) {
  x <- as.matrix(fit$data$counts)
  if (is.null(rownames(x))) rownames(x) <- fit$data$labels
  storage.mode(x) <- "double"
  last <- ncol(x)
  dt <- if (!is.null(fit$data$times)) diff(range(fit$data$times)) else sum(fit$data$dt)
  v0 <- x[x[, 1] > 0, 1]; v1 <- x[x[, last] > 0, last]
  cn0 <- od_parse(names(v0)); cn1 <- od_parse(names(v1))
  states <- 0:max(cn0, cn1)
  # y1[i, c]: cells at t1 with c copies of chromosome i
  y1 <- t(vapply(seq_len(ncol(cn1)), function(i) tapply(v1, factor(cn1[, i], levels = states), sum, default = 0), numeric(length(states))))
  y1[is.na(y1)] <- 0
  # indicator of each t0 karyotype's copy state on chromosome i: index into states
  idx0 <- cn0 - min(states) + 1L
  pr <- k$predict.karst_sparse_joint_fit(fit, names(v0))
  self <- 2 * k$sparse_joint_self_probability(cn0, fit$data$beta) - 1
  lb <- log(pr$fitness_mean / self)
  tr <- fit$trend
  des <- k$sparse_joint_design(cn0, tr$centre, tr$major, tr$shell_power)
  free <- which(des$group != 0)
  inside <- names(v0) %in% fit$summary$karyotype
  S <- tr$covariance[free, free, drop = FALSE]
  Sinv <- tryCatch(solve(S + diag(1e-10 * mean(diag(S)), nrow(S))), error = function(e) MASS::ginv(S))
  list(n0 = as.numeric(v0), lab0 = names(v0), self = self, lb = lb, D = des$design[, free, drop = FALSE] * inside,
       idx0 = idx0, y1 = y1, n1 = sum(v1), dt = dt, free = free, Sinv = Sinv, inside = inside,
       u = (y1 > 0) / pmax(rowSums(y1 > 0), 1))
}

# Predicted marginals at t1 from t0 growth rates `f` (net, per day).
mpl_marg <- function(S, f) {
  a <- S$n0 * exp((f - max(f)) * S$dt)
  A <- sum(a)
  p <- matrix(0, nrow(S$y1), ncol(S$y1))
  lev <- seq_len(ncol(p))
  for (i in seq_len(nrow(p))) p[i, ] <- as.numeric(tapply(a, factor(S$idx0[, i], levels = lev), sum, default = 0))
  list(p = p / A, a = a, A = A)
}

# PL per t1 cell for given t0 net growth rates `f` (any map: truth, KARST, neutral).
mpl_pl_rates <- function(S, f) {
  m <- mpl_marg(S, f)
  pt <- (1 - mpl_eps) * m$p + mpl_eps * S$u
  sum(S$y1 * log(pmax(pt, 1e-300))) / S$n1
}

# PL and its gradient in delta (D1a: only karyotypes in the domain move).
mpl_obj <- function(delta, S, w) {
  lb <- S$lb + drop(S$D %*% delta)
  f <- S$self * exp(lb)
  m <- mpl_marg(S, f)
  pt <- (1 - mpl_eps) * m$p + mpl_eps * S$u
  pl <- sum(S$y1 * log(pmax(pt, 1e-300)))
  # d PL / d a_k = (1 - eps) / A * sum_i [ y1[i, c_ki] / pt[i, c_ki] - sum_c y1[i, c] p[i, c] / pt[i, c] ]
  r <- S$y1 / pmax(pt, 1e-300)
  base <- rowSums(r * m$p)
  g_a <- numeric(length(m$a))
  for (i in seq_len(nrow(r))) g_a <- g_a + r[i, S$idx0[, i]] - base[i]
  g_a <- g_a * (1 - mpl_eps) / m$A
  g_delta <- drop(crossprod(S$D, g_a * m$a * S$dt * f))
  prior <- -0.5 * sum(delta * drop(S$Sinv %*% delta))
  list(value = prior + w * pl, grad = -drop(S$Sinv %*% delta) + w * g_delta, pl = pl / S$n1)
}

# Fit delta for weight w (BFGS on the negative objective).
mpl_fit <- function(S, w) {
  p <- length(S$free)
  if (!any(S$inside) || w == 0) return(list(delta = numeric(p), convergence = NA, pl = mpl_obj(numeric(p), S, 0)$pl))
  o <- stats::optim(numeric(p), function(d) -mpl_obj(d, S, w)$value, function(d) -mpl_obj(d, S, w)$grad,
                    method = "BFGS", control = list(maxit = 500, reltol = 1e-10))
  list(delta = o$par, convergence = o$convergence, pl = mpl_obj(o$par, S, 0)$pl, pl0 = mpl_obj(numeric(p), S, 0)$pl)
}

# predict() for `labels` after the update (D1a): every karyotype not in O gets
# log birth + Dtrend; O keeps KARST's level; the reference (majors weighted by
# their cells) is recomputed. sd columns unchanged. `rule_o = TRUE` (report
# only, D1b): O also get Dtrend.
mpl_predict <- function(k, fit, labels, delta, free, rule_o = FALSE) {
  base <- k$predict.karst_sparse_joint_fit(fit, labels)
  if (!any(delta != 0)) return(base)
  cn <- od_parse(labels); tr <- fit$trend
  des <- k$sparse_joint_design(cn, tr$centre, tr$major, tr$shell_power)$design[, free, drop = FALSE]
  move <- if (rule_o) rep(TRUE, length(labels)) else !labels %in% fit$dropped
  self <- 2 * k$sparse_joint_self_probability(cn, fit$data$beta) - 1
  lb <- log(base$fitness_mean / self) + ifelse(move, drop(des %*% delta), 0)
  base$fitness_mean <- self * exp(lb)
  s <- fit$summary; maj <- s$karyotype[s$major]
  wm <- rowSums(as.matrix(s[s$major, grep("^count_", names(s)), drop = FALSE]))
  cm <- od_parse(maj)
  dm <- k$sparse_joint_design(cm, tr$centre, tr$major, tr$shell_power)$design[, free, drop = FALSE]
  fm <- (2 * k$sparse_joint_self_probability(cm, fit$data$beta) - 1) * exp(s$log_birth_mean[s$major] + drop(dm %*% delta))
  base$relative_fitness <- base$fitness_mean - sum(wm * fm) / sum(wm)
  base
}
