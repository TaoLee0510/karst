# Version B of the paper's tests: alfak aligned to KARST (ALIGN=karst;
# docs/align_karst_6ae2a24.md).
#
# The default tests align KARST to alfak: for each minobs both methods report
# a landscape on alfak's karyotypes (the frequent karyotypes, with at least
# minobs cells, their one-step neighbours and those neighbours' neighbours),
# and every later step (forecasts, scores, Fig. 3-6) uses those karyotypes.
# Here the karyotypes come from KARST's side and do not depend on minobs: the
# extended list E of an input, as dev/realdata/extend builds it: every
# karyotype observed at a time point the fit uses, and their neighbours
# within two steps (one chromosome copy more or less per step, no copy number
# below 1). E holds alfak's lists of every minobs. On E:
#   alfak, for each minobs (its fits are kept): its landscape's mean on its
#     own list, and elsewhere the prediction of its stable Kriging model
#     (alfakR fitKrig()'s krig_stable_mean: fields::Krig on the bootstrap
#     means of the frequent karyotypes and their neighbours), rebuilt from
#     its bootstrap_res.Rds;
#   KARST, one landscape for every minobs: predict() on E, on the scale of
#     its minobs-10 landscape (as_alfak(): the level read back from it, the
#     efflux divisor of each karyotype).
# The sets that depended on minobs follow KARST: the frequent karyotypes
# (`fq`) are KARST's majors and `nn` their one-step neighbours, for both
# methods; `dist` is each karyotype's distance class to the observed ones
# (0, 1, 2), by which the simulations' landscape accuracy is reported
# (abm_fit_summary()). KARST's cross-validation score, which the paper's rule
# selects by, has folds centred on its majors and targets from each
# karyotype's own counts (as_alfak(xval_target = "own")); alfak keeps its
# own. Forecasts run alfak's ABM on E, as predict_evo() does, but are stored
# sparsely ("align_sparse_forecast": the frequencies above zero at each
# reported time), and are expanded for scoring over every karyotype of E
# that the forecasts or the observed counts hold (align_densify()), which
# gives the scores of the whole of E.
# With ALIGN=karst, bench_paths() puts the results in <round>/align_karst/,
# in the layout of a round, so that the collecting and analysis scripts run
# unchanged; `results_base` is the round itself, whose fits are read.

align_mode <- function() {
  a <- Sys.getenv("ALIGN", "")
  if (a %in% c("", "none")) return("none")
  if (!identical(a, "karst")) stop("ALIGN must be unset, \"none\" or \"karst\", not \"", a, "\"")
  a
}

# The extended list of an input count matrix (karyotypes x the time points
# fitted): E, and each member's distance class to the observed karyotypes.
align_list <- function(x) {
  nbr <- karst:::alfak_neighbours
  canon <- function(v) karst:::format_karyotypes(karst:::parse_karyotypes(v))
  observed <- unique(canon(rownames(x)[rowSums(x) > 0]))
  one <- setdiff(nbr(observed), observed)
  two <- setdiff(nbr(one), c(observed, one))
  E <- c(observed, one, two)
  list(E = E, dist = rep(c(0L, 1L, 2L), c(length(observed), length(one), length(two))))
}

# alfak's landscape of one fit folder on E: its mean on its own list, the
# stable Kriging prediction elsewhere. `agreement`: Pearson's r of the
# Kriging prediction and the mean on alfak's own list.
align_alfak_landscape <- function(dir, E, chunk = 20000L) {
  lsc <- readRDS(file.path(dir, "landscape.Rds"))
  if (!all(lsc$k %in% E)) stop("alfak's list is not inside E: ", dir)
  boot <- readRDS(file.path(dir, "bootstrap_res.Rds"))
  fboot <- cbind(boot$final_fitness, boot$nn_fitness)
  strs <- c(colnames(boot$final_fitness), colnames(boot$nn_fitness))
  fmean <- colMeans(fboot, na.rm = TRUE)[strs]
  valid <- is.finite(fmean)
  ktrain <- karst:::parse_karyotypes(strs)
  storage.mode(ktrain) <- "double"
  krig <- suppressWarnings(fields::Krig(ktrain[valid, , drop = FALSE], fmean[valid], cov.function = "stationary.cov",
                                        cov.args = alfakR:::krig_covariance_args(),
                                        nstep.cv = alfakR:::ALFAK_KRIG_NSTEP_CV, give.warnings = FALSE))
  predict_k <- function(k) {
    km <- karst:::parse_karyotypes(k)
    storage.mode(km) <- "double"
    idx <- split(seq_along(k), ceiling(seq_along(k) / chunk))
    unlist(lapply(idx, function(i) as.numeric(stats::predict(krig, km[i, , drop = FALSE]))), use.names = FALSE)
  }
  mean <- rep(NA_real_, length(E))
  inside <- match(lsc$k, E)
  mean[inside] <- lsc$mean
  outside <- setdiff(seq_along(E), inside)
  if (length(outside)) mean[outside] <- predict_k(E[outside])
  list(mean = mean, source = ifelse(seq_along(E) %in% inside, "landscape", "kriging"),
       agreement = stats::cor(predict_k(lsc$k), lsc$mean), krig = krig)
}

# KARST's landscape of a fit on E, on the scale of its stored minobs-10
# landscape (`dir10`): the level read back from it (one value, `level_spread`
# its largest deviation) and the efflux divisor.
align_karst_landscape <- function(fit, dir10, E, correct_efflux, chunk = 20000L) {
  lsc <- readRDS(file.path(dir10, "landscape.Rds"))
  divisor_of <- function(k) {
    if (isTRUE(correct_efflux)) 2 * karst:::sparse_joint_self_probability(karst:::parse_karyotypes(k), fit$data$beta) - 1
    else rep(1, length(k))
  }
  relative_of <- function(k) {
    idx <- split(seq_along(k), ceiling(seq_along(k) / chunk))
    unlist(lapply(idx, function(i) stats::predict(fit, k[i])$relative_fitness), use.names = FALSE)
  }
  probe <- lsc$k[seq_len(min(500L, nrow(lsc)))]
  lv <- lsc$mean[seq_along(probe)] * divisor_of(probe) - relative_of(probe)
  level <- stats::median(lv)
  mean <- (relative_of(E) + level) / divisor_of(E)
  m <- match(lsc$k, E)
  list(mean = mean, level = level, level_spread = max(abs(lv - level)),
       check_original = max(abs(mean[m] - lsc$mean)))
}

# A forecast on a list, as predict_evo(prediction_type = "ABM") runs it
# (same initial allocation, step grid, record interval and seed), read
# sparsely: for each reported time the frequencies above zero.
align_abm_forecast <- function(k, mean, x0, times, p, pop_size, delta_t, max_pop, culling, record_interval, seed) {
  init <- stats::setNames(numeric(length(k)), k)
  init[names(x0)] <- x0
  steps <- alfakR:::abm_times_to_steps(times, delta_t)
  n_steps <- max(steps)
  record <- alfakR:::resolve_abm_record_interval(steps, n_steps, record_interval)
  res <- alfakR:::run_karyotype_abm(initial_population_r = alfakR:::prepare_abm_initial_population(init, pop_size),
                                    fitness_map_r = stats::setNames(as.list(mean), k), p_missegregation = p,
                                    dt = delta_t, n_steps = as.integer(n_steps), max_population_size = max_pop,
                                    culling_survival_fraction = culling, record_interval = record,
                                    seed = as.integer(seed), grf_centroids = matrix(numeric(0), nrow = 0, ncol = 0),
                                    grf_lambda = NA_real_)
  freq <- lapply(steps, function(s) {
    v <- res[[as.character(s)]]
    if (!length(v) || sum(v) <= 0) return(stats::setNames(numeric(0), character(0)))
    stats::setNames(as.numeric(v) / sum(v), names(v))
  })
  structure(list(time = times, freq = freq), class = "align_sparse_forecast")
}

# Sparse forecast replicates as the dense data frames of predict_evo() (a
# time column, one column per karyotype), over the karyotypes any of them
# holds and those of `keep` (the observed karyotypes on the list): enough for
# every score of the paper, which are the same over the whole list.
align_densify <- function(preds, keep = character(0)) {
  support <- unique(c(unlist(lapply(preds, function(p) unlist(lapply(p$freq, names)))), keep))
  lapply(preds, function(p) {
    m <- matrix(0, length(p$time), length(support), dimnames = list(NULL, support))
    for (i in seq_along(p$time)) if (length(p$freq[[i]])) m[i, names(p$freq[[i]])] <- p$freq[[i]]
    cbind(data.frame(time = p$time), as.data.frame(m, check.names = FALSE))
  })
}
is_sparse_forecast <- function(x) inherits(x, "align_sparse_forecast")

# alfakR::find_steady_state(lscape, p, Nmax = 1), which Fig. 6 runs, for
# lists of any length. alfakR builds its matrix W (get_A_inputs) by visiting
# every pair of karyotypes and reserving n^2 entries, which fails beyond some
# 45,000 karyotypes, and E holds up to 1.6 million. With Nmax = 1, W links a
# karyotype only to itself and to its one-step neighbours on the list, so
# here those are found by lookup (keys: the counts as two base-B numbers, one
# complex number). Every entry is the same product, in the same order, of
# alfakR's per-chromosome probabilities (pij_cpp), and the eigenvector step
# is alfakR's, so W and the steady state are alfakR's
# (scripts/align_check_steady.R).
align_W <- function(k, p) {
  km <- karst:::parse_karyotypes(k)
  storage.mode(km) <- "integer"
  n <- nrow(km)
  nc <- ncol(km)
  top <- max(km)
  P <- matrix(0, top + 1L, top + 2L)
  for (a in 0:top) for (b in 0:(top + 1L)) P[a + 1L, b + 1L] <- alfakR:::pij_cpp(a, b, p)
  D <- matrix(P[cbind(as.vector(km) + 1L, as.vector(km) + 1L)], n, nc)
  # prefix products in alfakR's order: L[, c] = D[, 1] * ... * D[, c]
  L <- D
  if (nc > 1) for (m in 2:nc) L[, m] <- L[, m - 1L] * D[, m]
  B <- top + 2
  half <- ceiling(nc / 2)
  first <- seq_len(nc) <= half
  w <- B^(c(seq_len(half), seq_len(nc - half)) - 1)
  if (B^half >= 2^52) stop("copy numbers too high for exact keys")
  key <- complex(real = as.vector(km[, first, drop = FALSE] %*% w[first]),
                 imaginary = as.vector(km[, !first, drop = FALSE] %*% w[!first]))
  ii <- list(seq_len(n))
  jj <- list(seq_len(n))
  xx <- list(2 * L[, nc] - 1)
  for (c in seq_len(nc)) for (s in c(-1L, 1L)) {
    child <- km[, c] + s
    step <- if (first[c]) complex(real = s * w[c], imaginary = 0) else complex(real = 0, imaginary = s * w[c])
    j <- match(key + step, key)
    ok <- which(child >= 0L & !is.na(j))
    if (!length(ok)) next
    off <- P[cbind(km[ok, c] + 1L, child[ok] + 1L)]
    q <- if (c == 1L) off else L[ok, c - 1L] * off
    if (c < nc) for (m in (c + 1L):nc) q <- q * D[ok, m]
    ii[[length(ii) + 1L]] <- ok
    jj[[length(jj) + 1L]] <- j[ok]
    xx[[length(xx) + 1L]] <- 2 * q
  }
  ii <- unlist(ii)
  jj <- unlist(jj)
  xx <- unlist(xx)
  keep <- xx != 0
  Matrix::sparseMatrix(i = ii[keep], j = jj[keep], x = xx[keep], dims = c(n, n), dimnames = list(k, k), repr = "C")
}
align_find_steady_state <- function(lscape, p) {
  k <- lscape$k
  if (anyDuplicated(k)) stop("duplicate karyotypes")
  r_values <- lscape$mean
  names(r_values) <- k
  W <- align_W(k, p)
  M_eigs <- Matrix::t(W * r_values)
  dense <- function() {
    d <- eigen(as.matrix(M_eigs))
    i <- which.max(Re(d$values))
    list(values = d$values[i], vectors = d$vectors[, i, drop = FALSE])
  }
  eig <- if (nrow(M_eigs) < 3) dense() else tryCatch(RSpectra::eigs(M_eigs, k = 1, which = "LR"), error = function(e) {
    warning("RSpectra::eigs with 'LR' failed; falling back to dense eigen decomposition: ", conditionMessage(e),
            call. = FALSE)
    dense()
  })
  ss <- Re(eig$vectors[, 1])
  ss <- ss * sign(ss[which.max(abs(ss))[1]])
  ss[ss < 1e-10] <- 0
  if (is.na(sum(ss)) || sum(ss) <= 1e-9) stop("steady state sums to ", sum(ss))
  ss <- ss / sum(ss)
  names(ss) <- k
  ss
}
