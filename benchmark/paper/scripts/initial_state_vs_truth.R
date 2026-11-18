# karst's two-time-point fits of this round (initial_state = "auto", 827c4b5
# on) against the same inputs' fits of a reference round with the balance
# prior (05a8855), each against the truth, on the ABM and missegregation
# simulations, by the state "auto" chose. Read only; refits nothing.
#
#   Rscript scripts/initial_state_vs_truth.R REFERENCE_RESULTS
#
# Per fit (input, design, minobs; both rounds' fits ok):
#   landscape  r, rho, R (every karyotype of the landscape), rnn and Rnn (one
#              copy-number step from a frequent karyotype), rd2 (two steps), from
#              score.Rds
#   CV         r_heldout_true: the cross-validation's held-out predictions
#              against the true fitness of those karyotypes (xval.Rds; as
#              scripts/cv_score_analysis.R); Rxv: the CV score as computed
#   forecast   ABM only, the paper's Fig. 1g scores against the no-change
#              baseline at passage 10 after the last time point (and the mean
#              over passages 1-10): overlap and cosine (higher is better),
#              Euclidean and Wasserstein distance (lower is better), and the
#              angle of the predicted change to the true one (degrees)
# and the check of initial_state = "auto" (tables/initial_state.tsv of this
# round: level, choice, the check fit's convergence) and both fits'
# convergence.
#
# Writes BENCH_RESULTS/analysis/initial_state/vs_truth_fits.tsv (one row per
# fit) and vs_truth.tsv (paired summaries by source, design, wavelength or
# true rate, and choice), and prints them with the check fits' convergence.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
reference <- normalizePath(commandArgs(TRUE)[[1]])
out <- file.path(paths$results, "analysis", "initial_state")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 220)
centred_r <- function(a, b) if (length(a) >= 3 && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_
st <- utils::read.delim(file.path(paths$results, "tables", "initial_state.tsv"), stringsAsFactors = FALSE)
st <- st[st$source %in% c("abm", "misseg"), ]

fit_row <- function(root, d, truth_of) {
  run <- file.path(d, "run.Rds")
  if (!file.exists(run)) return(NULL)
  r <- readRDS(run)
  if (!identical(r$status, "ok") || !file.exists(file.path(d, "score.Rds"))) return(NULL)
  sc <- readRDS(file.path(d, "score.Rds"))
  s <- sc$summary
  row <- data.frame(r = s$r, rho = s$rho, R = s$R, rnn = s$rnn, Rnn = s$Rnn, rd2 = s$rd2, Rxv = s$Rxv,
                    convergence = first_or_na(r$convergence), history = first_or_na(r$history))
  x <- readRDS(file.path(d, "xval.Rds"))
  p <- unname(x$predictions)
  k <- x$xval_data$k
  ok <- is.finite(p)
  row$r_heldout_true <- if (sum(ok) >= 3) centred_r(truth_of(k[ok]), p[ok]) else NA_real_
  if (!is.null(sc$df_de)) {
    i10 <- nrow(sc$df_de)  # passage 10 after the last time point (rows: passages 0-10)
    pass <- seq(2L, i10)   # passages 1-10
    row$overlap_win10 <- as.numeric(sc$df_so$so_p[i10] > sc$df_so$so_0[i10])
    row$cosine_win10 <- as.numeric(sc$df_sc$sc_p[i10] > sc$df_sc$sc_0[i10])
    row$euclid_win10 <- as.numeric(sc$df_de$de_p[i10] < sc$df_de$de_0[i10])
    row$wasser_win10 <- as.numeric(sc$df_dw$dw_p[i10] < sc$df_dw$dw_0[i10])
    row$euclid_win_mean <- mean(sc$df_de$de_p[pass] < sc$df_de$de_0[pass])
    row$wasser_win_mean <- mean(sc$df_dw$dw_p[pass] < sc$df_dw$dw_0[pass])
    row$euclid_ratio10 <- sc$df_de$de_p[i10] / sc$df_de$de_0[i10]
    row$angle10 <- sc$df_a$angle[i10]
  }
  row
}

rows <- list()
for (i in seq_len(nrow(st))) {
  s <- st[i, ]
  raw <- if (s$source == "abm") paths$abm_raw else paths$misseg_raw
  sim <- read_abm_sim(s$dataset, raw)
  info <- abm_sim_info(s$dataset)
  truth_of <- function(k) grf_fitness(k, info$w, sim$true_landscape)
  for (m in Sys.glob(file.path(paths$results, "fits", s$source, s$dataset, s$design, "karst", "minobs_*"))) {
    mo <- sub("^minobs_", "", basename(m))
    new <- fit_row(paths$results, m, truth_of)
    old <- fit_row(reference, file.path(reference, "fits", s$source, s$dataset, s$design, "karst", basename(m)), truth_of)
    if (is.null(new) || is.null(old)) next
    names(new) <- paste0(names(new), ".new")
    names(old) <- paste0(names(old), ".old")
    rows[[length(rows) + 1L]] <- cbind(data.frame(source = s$source, dataset = s$dataset, design = s$design,
                                                  minobs = mo, lambda = info$w, rate = info$pmis, chosen = s$chosen,
                                                  level = s$level, check_convergence = s$check_convergence,
                                                  stringsAsFactors = FALSE), new, old)
  }
}
fits <- bind_rows_fill(rows)
utils::write.table(fits, file.path(out, "vs_truth_fits.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

higher <- c(r = TRUE, rho = TRUE, R = TRUE, rnn = TRUE, Rnn = TRUE, rd2 = TRUE, r_heldout_true = TRUE, Rxv = TRUE,
            overlap_win10 = TRUE, cosine_win10 = TRUE, euclid_win10 = TRUE, wasser_win10 = TRUE,
            euclid_win_mean = TRUE, wasser_win_mean = TRUE, euclid_ratio10 = FALSE, angle10 = FALSE)
summarise <- function(d) {
  do.call(rbind, lapply(names(higher), function(v) {
    a <- d[[paste0(v, ".new")]]
    b <- d[[paste0(v, ".old")]]
    if (is.null(a) || is.null(b)) return(NULL)
    ok <- is.finite(a) & is.finite(b)
    if (!any(ok)) return(NULL)
    a <- a[ok]
    b <- b[ok]
    diff <- if (higher[[v]]) a - b else b - a
    data.frame(score = v, pairs = sum(ok), new = round(mean(a), 3), old = round(mean(b), 3),
               median_new = round(stats::median(a), 3), median_old = round(stats::median(b), 3),
               median_diff = round(stats::median(a - b), 3), better = round(mean(diff > 1e-12), 2),
               worse = round(mean(diff < -1e-12), 2),
               p = if (sum(ok) > 5 && any(abs(diff) > 1e-12)) signif(suppressWarnings(
                 stats::wilcox.test(a, b, paired = TRUE)$p.value), 2) else NA_real_)
  }))
}
fits$group <- ifelse(fits$source == "abm", paste("lambda", fits$lambda), paste("rate", fits$rate))
summ <- do.call(rbind, lapply(split(fits, list(fits$source, fits$design, fits$group, fits$chosen), drop = TRUE), function(d) {
  s <- summarise(d)
  if (is.null(s)) return(NULL)
  cbind(data.frame(source = d$source[1], design = d$design[1], group = d$group[1], chosen = d$chosen[1],
                   inputs = length(unique(d$dataset)), stringsAsFactors = FALSE), s)
}))
pooled <- do.call(rbind, lapply(split(fits, list(fits$source, fits$design, fits$chosen), drop = TRUE), function(d) {
  s <- summarise(d)
  if (is.null(s)) return(NULL)
  cbind(data.frame(source = d$source[1], design = d$design[1], group = "all", chosen = d$chosen[1],
                   inputs = length(unique(d$dataset)), stringsAsFactors = FALSE), s)
}))
summ <- rbind(pooled, summ)
summ$reference <- basename(reference)
utils::write.table(summ, file.path(out, "vs_truth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("fits paired:", nrow(fits), " (this round against", basename(reference), ")\n")
show <- function(sel, title) {
  cat("\n==", title, "\n")
  x <- summ[sel, setdiff(names(summ), c("reference", "new", "old"))]
  print(x, row.names = FALSE)
}
key_scores <- c("r", "rnn", "rd2", "r_heldout_true", "euclid_win10", "wasser_win10", "euclid_win_mean", "angle10")
show(summ$design == "w2_ends" & summ$chosen == "free" & summ$score %in% key_scores,
     "w2_ends, inputs switched to free: this round (free) against the reference (balance)")
show(summ$design == "w2_ends" & summ$chosen == "balance" & summ$group == "all" & summ$score %in% key_scores,
     "w2_ends, inputs kept on balance (should not change)")
show(summ$design != "w2_ends" & summ$chosen == "free" & summ$score %in% key_scores,
     "w4/w8_ends, inputs switched to free")

cat("\n== The check fit's convergence (1 = failed the convergence check) by choice\n")
st$group <- ifelse(st$source == "abm", paste("lambda", st$lambda), paste("rate", st$rate))
cv <- stats::aggregate(cbind(inputs = 1, check_unconverged = check_convergence %in% 1) ~ source + design + chosen,
                       data = st[!is.na(st$chosen), ], FUN = sum)
cv$share <- round(cv$check_unconverged / cv$inputs, 3)
print(cv[order(cv$source, cv$design, cv$chosen), ], row.names = FALSE)
cv2 <- stats::aggregate(cbind(inputs = 1, check_unconverged = check_convergence %in% 1) ~ source + group + chosen,
                        data = st[!is.na(st$chosen) & st$design == "w2_ends", ], FUN = sum)
cv2$share <- round(cv2$check_unconverged / cv2$inputs, 3)
cat("\nw2_ends by wavelength / true rate:\n")
print(cv2[order(cv2$source, cv2$group, cv2$chosen), ], row.names = FALSE)
cat("\nLevel of switched inputs whose check did / did not converge (w2_ends):\n")
w2 <- st[!is.na(st$chosen) & st$design == "w2_ends" & st$chosen == "free", ]
if (nrow(w2)) {
  print(stats::aggregate(level ~ source + check_convergence, w2, function(v) round(stats::quantile(v, c(.1, .5, .9)), 2)))
}
cat("\nConvergence of the fits themselves (1 = failed), w2_ends, this round vs reference, by choice:\n")
fc <- fits[fits$design == "w2_ends", ]
print(stats::aggregate(cbind(new = convergence.new %in% 1, old = convergence.old %in% 1, fits = 1) ~ source + chosen,
                       fc, sum), row.names = FALSE)
