# Scores of a fitted ABM landscape, as S04_process_abm_results.R of the paper
# computes them. abm_fit_summary() returns the same list that S04 saves as
# abm_preds_summary.Rds (info, summary, df_so, df_sc, df_de, df_dw, df_a), so
# the paper's own outputs can check this code.

# Rescaled R-squared of S04 (both vectors centred first).
r2r <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

# Cross-validation score from xval.Rds: a number in the paper's alfak, a
# list with R2R in alfak dev and karst.
xval_score <- function(xval) {
  if (is.list(xval)) xval$R2R else as.numeric(xval)
}

# Accuracy of the landscape against the true fitness `f` of its karyotypes.
landscape_scores <- function(lscape, f, xval) {
  d2 <- !(lscape$fq | lscape$nn)
  fq <- lscape$fq
  nn <- lscape$nn
  m <- lscape$mean
  data.frame(r = cor(m, f), rfq = cor(m[fq], f[fq]), rnn = cor(m[nn], f[nn]), rd2 = cor(m[d2], f[d2]),
             rho = cor(m, f, method = "spearman"), rhofq = cor(m[fq], f[fq], method = "spearman"),
             rhonn = cor(m[nn], f[nn], method = "spearman"), rhod2 = cor(m[d2], f[d2], method = "spearman"),
             R = r2r(f, m), Rfq = r2r(f[fq], m[fq]), Rnn = r2r(f[nn], m[nn]), Rd2 = r2r(f[d2], m[d2]),
             Rxv = xval_score(xval), nfq = sum(fq))
}

# Forecast scores against the true frequencies at the sampling times from the
# forecast start onwards. `preds` is the list of forecast replicates (data
# frames with a time column and one column per landscape karyotype, rows at
# the start and at each later sampling time); `sim_x` the full count matrix of
# the simulation.
forecast_scores <- function(lscape, preds, sim_x, last_time = 120) {
  preds <- as.data.frame(Reduce(`+`, lapply(preds, function(df) as.matrix(df[-1]))) / length(preds))
  tt <- as.numeric(colnames(sim_x))
  t0 <- max(tt[tt < last_time])
  measure_tps <- utils::head(tt[tt >= t0], 11)
  yi <- sim_x[, colnames(sim_x) %in% measure_tps]
  y <- preds * 0
  rx <- rownames(yi)[rownames(yi) %in% lscape$k]
  for (i in seq_along(measure_tps)) y[i, rx] <- yi[rx, i] / sum(yi[rx, i])
  y0 <- as.numeric(y[1, ])
  k <- do.call(rbind, lapply(strsplit(lscape$k, ".", fixed = TRUE), as.numeric))
  k0 <- colSums(y0 * k)
  mag <- function(v) sqrt(sum(v^2))
  angle <- function(a, b) 180 * acos(sum(a * b) / (mag(a) * mag(b))) / pi
  rows <- seq_len(nrow(y))
  yr <- function(i) as.numeric(y[i, ])
  pr <- function(i) as.numeric(preds[i, ])

  df_a <- data.frame(angle = vapply(rows, function(i) angle(colSums(yr(i) * k) - k0, colSums(pr(i) * k) - k0), 0))
  df_de <- do.call(rbind, lapply(rows, function(i) {
    ky <- colSums(yr(i) * k)
    kp <- colSums(pr(i) * k)
    data.frame(de_p = sqrt(sum((kp - ky)^2)), de_0 = sqrt(sum((k0 - ky)^2)))
  }))
  s0 <- transport::wpp(k, mass = y0)
  df_dw <- do.call(rbind, lapply(rows, function(i) {
    tryCatch({
      sy <- transport::wpp(k, mass = yr(i))
      sp <- transport::wpp(k, mass = pr(i))
      data.frame(dw_p = transport::wasserstein(sy, sp), dw_0 = transport::wasserstein(sy, s0))
    }, error = function(e) data.frame(dw_p = NaN, dw_0 = NaN))
  }))
  df_sc <- do.call(rbind, lapply(rows, function(i) {
    a <- yr(i)
    p <- pr(i)
    data.frame(sc_p = sum(p * a) / prod(sqrt(sum(p^2)), sqrt(sum(a^2))),
               sc_0 = sum(y0 * a) / prod(sqrt(sum(y0^2)), sqrt(sum(a^2))))
  }))
  df_so <- do.call(rbind, lapply(rows, function(i) {
    a <- yr(i)
    p <- pr(i)
    data.frame(so_p = sum(pmin(p, a)) / min(sum(p), sum(a)), so_0 = sum(pmin(y0, a)) / min(sum(y0), sum(a)))
  }))
  rownames(df_a) <- rownames(df_dw) <- rownames(df_sc) <- rownames(df_so) <- measure_tps
  list(df_so = df_so[c("so_p", "so_0")], df_sc = df_sc[c("sc_p", "sc_0")], df_de = df_de,
       df_dw = df_dw, df_a = df_a)
}

# Everything S04 saves for one fit directory. `info` identifies the fit.
# `with_forecast`: the forecasts are complete (run_task.R passes
# forecast_complete(), R/forecast.R).
abm_fit_summary <- function(fit_dir, sim, info, last_time = 120, with_forecast = TRUE) {
  lscape <- readRDS(file.path(fit_dir, "landscape.Rds"))
  xval <- readRDS(file.path(fit_dir, "xval.Rds"))
  f <- grf_fitness(lscape$k, info$w, sim$true_landscape)
  out <- list(info = info, summary = landscape_scores(lscape, f, xval))
  if (with_forecast) {
    paths <- list.files(file.path(fit_dir, "abm_preds"), full.names = TRUE)
    if (!length(paths)) stop("No forecasts in ", fit_dir)
    out <- c(out, forecast_scores(lscape, lapply(paths, readRDS), sim$abm_output$x, last_time))
  }
  out
}
