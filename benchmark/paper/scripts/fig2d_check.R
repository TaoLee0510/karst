#!/usr/bin/env Rscript
# Why does karst with two time points rank the size of change badly in the
# paper's Fig. 2d? (Salehi, round 05a8855: Spearman 0.25 for karst against
# 0.45 for alfak, both on two time points.) A diagnosis only: it reads the
# fits and forecasts of the round and refits nothing.
#
#   Rscript scripts/fig2d_check.R --stage=salehi --index N   chunk N of 10 trajectories
#   Rscript scripts/fig2d_check.R --stage=abm --index N      chunk N of 4 simulations
#
# For every two-time-point fit ("ends" design) of karst and alfak, every
# minobs, and every daughter passage (Salehi) or the next sampling time (ABM),
# one row with:
#   base, pred0      Fig. 2d's two axes as the paper computes them: the
#                    Wasserstein distance from the parent to the daughter
#                    (observed change) and from the parent to the forecast
#                    (predicted change), both restricted to the landscape's
#                    karyotypes and renormalised. Salehi: at the paper's
#                    evaluation day (15, or 25 for the p53 lines).
#   base_full,       the same on every observed karyotype; the parent's cells
#   pred0_full       outside the landscape are kept unchanged in the forecast.
#   sel0             the change selection alone would make: the parent's
#                    frequencies times exp(fitness x time), renormalised (no
#                    missegregation, no drift).
#   var_f, price     fitness variance among the parent's cells and the norm of
#                    the Price vector (covariance of fitness and copy number),
#                    which set how fast the population and its mean karyotype
#                    move under selection.
# ABM rows also carry the same quantities from the true landscape.
#
# Output: BENCH_RESULTS/analysis/fig2d_check/{salehi,abm}/chunk_NNN.Rds. The
# analysis is scripts/fig2d_check_analyse.R.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
stage <- arg("stage")
index <- as.integer(arg("index"))
# FIG2D_OUT: another output folder, to read a results folder without writing
# into it.
out <- Sys.getenv("FIG2D_OUT")
if (!nzchar(out)) out <- file.path(paths$results, "analysis", "fig2d_check", stage)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

s2m <- function(k) do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
# Wasserstein distance (as the paper: Euclidean ground distance, p = 1)
# between two distributions given as masses on karyotype labels.
wdist <- function(ka, ma, kb, mb) {
  a <- is.finite(ma) & ma > 0
  b <- is.finite(mb) & mb > 0
  if (!any(a) || !any(b)) return(NaN)
  tryCatch(transport::wasserstein(transport::wpp(s2m(ka[a]), mass = ma[a] / sum(ma[a])),
                                  transport::wpp(s2m(kb[b]), mass = mb[b] / sum(mb[b]))),
           error = function(e) NaN)
}
entropy <- function(p) {
  p <- p[p > 0] / sum(p)
  -sum(p * log(p))
}
# Fitness variance and Price vector among cells with weights w.
selection_stats <- function(f, w, kmat, prefix) {
  w <- w / sum(w)
  fbar <- sum(w * f)
  kbar <- colSums(w * kmat)
  price <- colSums(w * (f - fbar) * sweep(kmat, 2L, kbar))
  r <- list(var_f = sum(w * (f - fbar)^2), price = sqrt(sum(price^2)))
  names(r) <- paste0(prefix, names(r))
  list(stats = r, price = price)
}
angle <- function(a, b) {
  m <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  if (!is.finite(m) || m == 0) return(NaN)
  180 * acos(max(-1, min(1, sum(a * b) / m))) / pi
}
mean_preds <- function(dir) {
  paths <- sort(list.files(file.path(dir, "abm_preds"), full.names = TRUE))
  if (!length(paths)) return(NULL)
  preds <- lapply(paths, readRDS)
  mat <- Reduce(`+`, lapply(preds, function(x) as.matrix(x[, -1]))) / length(preds)
  list(time = preds[[1]]$time, mat = mat)
}
fit_ok <- function(dir) {
  run <- file.path(dir, "run.Rds")
  file.exists(run) && identical(readRDS(run)$status, "ok") && file.exists(file.path(dir, "landscape.Rds"))
}

# One row: the parent `p` and daughter `d` (counts named by karyotype, every
# observed karyotype), the landscape, the forecast `fc` at the evaluation time
# (on the landscape's karyotypes) and the elapsed time `dt`.
change_row <- function(p, d, lscape, fc, dt) {
  k_all <- union(names(p), names(d))
  pa <- stats::setNames(numeric(length(k_all)), k_all)
  da <- pa
  pa[names(p)] <- p
  da[names(d)] <- d
  inl <- k_all %in% lscape$k
  yb <- pa * inl / sum(pa[inl])
  yt <- da * inl / sum(da[inl])
  kb <- names(yb)
  ix <- match(kb[inl & yb > 0], lscape$k)
  f <- lscape$mean[ix]
  wpar <- yb[inl & yb > 0]
  kmat <- s2m(kb[inl & yb > 0])
  sel <- wpar * exp(f * dt)
  st <- selection_stats(f, wpar, kmat, "")
  # mean karyotypes (Euclidean axes) on the landscape's karyotypes
  k_l <- s2m(lscape$k)
  kbar_b <- colSums(wpar * kmat) / sum(wpar)
  kt_k <- kb[inl & yt > 0]
  kbar_t <- colSums(yt[inl & yt > 0] * s2m(kt_k)) / sum(yt[inl & yt > 0])
  kbar_p <- colSums(fc * k_l) / sum(fc)
  # the parent's cells outside the landscape kept unchanged
  cov_p <- sum(pa[inl]) / sum(pa)
  aug_k <- c(lscape$k, k_all[!inl])
  aug_m <- c(cov_p * fc / sum(fc), pa[!inl] / sum(pa))
  data.frame(base = wdist(kb, yt, kb, yb), pred0 = wdist(kb, yb, lscape$k, fc), pred = wdist(kb, yt, lscape$k, fc),
             sel0 = wdist(kb[inl & yb > 0], wpar, kb[inl & yb > 0], sel),
             base_full = wdist(k_all, da, k_all, pa), pred0_full = wdist(k_all, pa, aug_k, aug_m),
             base_euc = sqrt(sum((kbar_t - kbar_b)^2)), pred0_euc = sqrt(sum((kbar_p - kbar_b)^2)),
             angle_pred = angle(kbar_t - kbar_b, kbar_p - kbar_b), angle_price = angle(kbar_t - kbar_b, st$price),
             var_f = st$stats$var_f, price = st$stats$price, sd_f_range = diff(range(f)),
             cov_parent = cov_p, cov_daughter = sum(da[inl]) / sum(da),
             cells_parent = sum(pa), cells_daughter = sum(da), kary_parent = sum(pa > 0),
             kary_parent_in = sum(pa[inl] > 0), entropy_parent = entropy(yb[inl]),
             landscape_size = nrow(lscape), n_fq = sum(lscape$fq), dt = dt)
}

karst_extras <- function(dir, p) {
  out <- data.frame(n_major = NA_integer_, history = NA_real_, domain_parent = NA_real_, fit_dt = NA_real_,
                    cells_t0 = NA_real_, cells_t1 = NA_real_, birth_rate = NA_real_, balance_level = NA_real_,
                    shell_coef = NA_real_, scale_residual = NA_real_, intermediate_weight = NA_real_,
                    samples_fitted = NA_real_)
  path <- file.path(dir, "karst_fit.Rds")
  if (!file.exists(path)) return(out)
  fit <- readRDS(path)
  counts <- fit$data$counts
  out$n_major <- fit$diagnostics$n_major
  out$history <- fit$diagnostics$history
  out$domain_parent <- sum(p[names(p) %in% fit$summary$karyotype]) / sum(p)
  out$fit_dt <- fit$data$dt
  out$cells_t0 <- sum(counts[, 1L])
  out$cells_t1 <- sum(counts[, ncol(counts)])
  out$birth_rate <- fit$diagnostics$birth_rate
  # A candidate build's level of the balance prior (karst's
  # validation/sparse_joint/cloud_level), and what the history acts through.
  if (!is.null(fit$diagnostics$balance_level)) out$balance_level <- fit$diagnostics$balance_level
  out$shell_coef <- fit$trend$coefficients[[length(fit$trend$coefficients)]]
  out$scale_residual <- fit$trend$scales[["residual"]]
  out$intermediate_weight <- if (is.null(fit$diagnostics$intermediate_weight)) NA_real_ else fit$diagnostics$intermediate_weight
  out$samples_fitted <- ncol(counts)
  out
}

if (stage == "salehi") {
  meta <- utils::read.csv(paths$salehi_meta)
  lineages <- readRDS(paths$salehi_lineages)
  treat_lut <- stats::setNames(meta$on_treatment, paste0("u", meta$uid))
  inputs <- list.files(paths$salehi_inputs)
  all_traj <- sort(salehi_trajectories(paths)$traj_id)
  chunk <- all_traj[seq((index - 1L) * 10L + 1L, min(length(all_traj), index * 10L))]
  rows <- list()
  for (id in chunk) {
    li <- lineages[[id]]
    type <- salehi_lineage_fig2(id, meta)
    maxp <- if (isTRUE(type %in% c("p53 k.o", "p53 w.t"))) 25 else 15
    # FIG2D_METHODS: other fitted methods to read instead, e.g. the candidate
    # arms of karst's cloud-level study in their own results folder.
    for (method in strsplit(Sys.getenv("FIG2D_METHODS", "karst,alfak"), ",", fixed = TRUE)[[1L]]) for (mo in c(5, 10, 20)) {
      # FIG2D_DESIGN: "ends" (default) or "all" (every passage of the trajectory).
      dir <- file.path(paths$results, "fits", "salehi", id, Sys.getenv("FIG2D_DESIGN", "ends"), method,
                       paste0("minobs_", mo))
      if (!fit_ok(dir)) next
      pr <- mean_preds(dir)
      if (is.null(pr) || !any(pr$time == maxp) || !any(pr$time == 1)) next
      lscape <- readRDS(file.path(dir, "landscape.Rds"))
      fc <- as.numeric(pr$mat[which(pr$time == maxp)[1], ])
      if (!identical(colnames(pr$mat), lscape$k)) fc <- fc[match(lscape$k, colnames(pr$mat))]
      for (dj in li$dec1) {
        fj <- paste(meta[meta$uid == dj, c("datasetname", "timepoint")], collapse = "_")
        if (sum(grepl(fj, inputs)) < 1) next
        tru <- readRDS(file.path(paths$salehi_inputs, inputs[grepl(fj, inputs)][1]))$x
        p <- stats::setNames(tru[, ncol(tru) - 1L], rownames(tru))
        d <- stats::setNames(tru[, ncol(tru)], rownames(tru))
        p <- p[p > 0]
        d <- d[d > 0]
        if (!length(p) || !length(d) || !any(names(p) %in% lscape$k) || !any(names(d) %in% lscape$k)) next
        r <- change_row(p, d, lscape, fc, maxp)
        # The trajectory's own recent and total change, which a forecast could
        # follow: from the passage before the parent to the parent (on the
        # landscape's karyotypes, renormalised, as `base`; and on all cells),
        # and from the trajectory's first passage to the parent (all cells).
        traj <- salehi_input(id, ends = FALSE, paths)$x
        r$passages <- ncol(traj)
        r$base_prev <- r$base_prev_all <- r$span_all <- NA_real_
        if (ncol(traj) >= 2L) {
          q <- stats::setNames(traj[, ncol(traj) - 1L], rownames(traj))
          q <- q[q > 0]
          pp <- stats::setNames(traj[, ncol(traj)], rownames(traj))
          pp <- pp[pp > 0]
          ql <- q[names(q) %in% lscape$k]
          pl <- pp[names(pp) %in% lscape$k]
          r$base_prev_all <- wdist(names(q), q, names(pp), pp)
          if (length(ql) && length(pl)) r$base_prev <- wdist(names(ql), ql, names(pl), pl)
          f1 <- stats::setNames(traj[, 1L], rownames(traj))
          f1 <- f1[f1 > 0]
          r$span_all <- wdist(names(f1), f1, names(pp), pp)
        }
        extra <- if (startsWith(method, "karst")) karst_extras(dir, p) else karst_extras(file.path(dirname(dir), "none"), p)
        rows[[length(rows) + 1L]] <- cbind(data.frame(fi = id, method = method, min_obs = as.character(mo),
                                                      dj = dj, type = type, tt = maxp,
                                                      test_treat = as.character(treat_lut[paste0("u", dj)]),
                                                      stringsAsFactors = FALSE), r, extra)
      }
    }
    cat(id, "done\n")
  }
  saveRDS(do.call(rbind, rows), file.path(out, sprintf("chunk_%03d.Rds", index)))
}

if (stage == "abm") {
  sims <- sort(sub("\\.Rds$", "", list.files(paths$abm_raw, pattern = "\\.Rds$")))
  chunk <- sims[seq((index - 1L) * 4L + 1L, min(length(sims), index * 4L))]
  rows <- list()
  for (sim_id in chunk) {
    info <- abm_sim_info(sim_id)
    sim <- read_abm_sim(sim_id, paths$abm_raw)
    x <- as.matrix(sim$abm_output$x)  # a data frame in some simulations; keeps the karyotype row names
    times <- as.numeric(colnames(x))
    t0 <- max(times[times < 120])
    t1 <- min(times[times > t0])
    p <- x[, times == t0]
    d <- x[, times == t1]
    p <- p[p > 0]
    d <- d[d > 0]
    f_true_p <- stats::setNames(grf_fitness(names(p), info$w, sim$true_landscape), names(p))
    for (design in c("w2_ends", "w4_ends", "w8_ends")) for (method in c("karst", "alfak")) for (mo in c(5, 10, 20)) {
      dir <- file.path(paths$results, "fits", "abm", sim_id, design, method, paste0("minobs_", mo))
      if (!fit_ok(dir)) next
      pr <- mean_preds(dir)
      if (is.null(pr) || nrow(pr$mat) < 2L) next
      lscape <- readRDS(file.path(dir, "landscape.Rds"))
      fc <- as.numeric(pr$mat[2L, ])
      if (!identical(colnames(pr$mat), lscape$k)) fc <- fc[match(lscape$k, colnames(pr$mat))]
      if (!any(names(p) %in% lscape$k) || !any(names(d) %in% lscape$k)) next
      r <- change_row(p, d, lscape, fc, t1 - t0)
      # the same from the true landscape, on the parent's cells in the landscape
      inl <- names(p) %in% lscape$k
      w <- p[inl] / sum(p[inl])
      kmat <- s2m(names(w))
      f_pred <- lscape$mean[match(names(w), lscape$k)]
      f_true <- f_true_p[names(w)]
      st_true <- selection_stats(f_true, w, kmat, "true_")
      st_pred <- selection_stats(f_pred, w, kmat, "")
      sel_true <- w * exp(f_true * (t1 - t0))
      truth <- data.frame(true_var_f = st_true$stats$true_var_f, true_price = st_true$stats$true_price,
                          true_sel0 = wdist(names(w), w, names(w), sel_true),
                          angle_price_truth = angle(st_pred$price, st_true$price),
                          spearman_f_parent = suppressWarnings(stats::cor(f_pred, f_true, method = "spearman")))
      rows[[length(rows) + 1L]] <- cbind(data.frame(sim_id = sim_id, w = info$w, design = design, method = method,
                                                    min_obs = as.character(mo), t0 = t0, t1 = t1,
                                                    stringsAsFactors = FALSE), r, truth,
                                         if (method == "karst") karst_extras(dir, p) else karst_extras("none", p))
    }
    cat(sim_id, "done\n")
  }
  saveRDS(do.call(rbind, rows), file.path(out, sprintf("chunk_%03d.Rds", index)))
}
