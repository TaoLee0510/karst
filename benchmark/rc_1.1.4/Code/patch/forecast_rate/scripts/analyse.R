# Study B (PLAN.md): Q1-Q4 on the simulations and the Salehi description.
#   Rscript scripts/analyse.R OUT_DIR RES_DIR
args <- commandArgs(TRUE); out_dir <- args[[1]]; res_dir <- args[[2]]
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
options(width = 220)
B <- 2000L; seed <- 20261125L

## simulations -----------------------------------------------------------------------
L <- list(); status <- list(); A <- list()
for (f in list.files(out_dir, "[.]rds$", recursive = TRUE, full.names = TRUE)) {
  if (!grepl("/(sim|alfak)_(abm|misseg)/", f) || grepl("/pilot/", f)) next
  x <- readRDS(f); t <- x$task
  if (grepl("/alfak_", f)) { if (identical(x$status, "ok")) A[[length(A) + 1L]] <- cbind(data.frame(source = t$source, sim_id = t$sim_id, design = t$design, rep = t$rep, half = t$half,
    unit = paste(t$source, t$sim_id), w = x$info$w, true_mult = round(x$true_mult), horizon = NA_real_, stringsAsFactors = FALSE), x$rows); next }
  status[[length(status) + 1L]] <- data.frame(source = t$source, sim_id = t$sim_id, design = t$design, status = x$status,
    same_input = if (is.null(x$same_input)) NA else all(x$same_input), fits = if (is.null(x$fits_found)) "" else paste(x$fits_found, collapse = ","))
  if (!identical(x$status, "ok")) next
  L[[length(L) + 1L]] <- cbind(data.frame(source = t$source, sim_id = t$sim_id, design = t$design, rep = t$rep, half = t$half,
    unit = paste(t$source, t$sim_id), w = x$info$w, true_mult = round(x$true_mult), horizon = x$horizon, stringsAsFactors = FALSE), x$rows)
}
status <- do.call(rbind, status); L <- do.call(rbind, L)
A <- if (length(A)) do.call(rbind, A) else NULL
cat("ALFA-K arm outputs:", if (is.null(A)) 0 else length(unique(paste(A$source, A$sim_id, A$design))), "inputs\n")
if (!is.null(A)) { A$horizon <- L$horizon[match(paste(A$source, A$sim_id, A$design), paste(L$source, L$sim_id, L$design))]; L <- rbind(L, A[, names(L)]) }
cat("simulated inputs:", nrow(status), " ok", sum(status$status == "ok"), " inputs whose fits differ from the rebuilt input:",
    sum(status$same_input %in% FALSE), "\n"); print(table(status$source, status$fits))
L$input <- paste(L$source, L$sim_id, L$design)
L$set <- ifelse(L$source == "abm", "main", paste0("vM x", L$true_mult))
utils::write.table(L, file.path(res_dir, "sim_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# paired comparison of arm a against arm b on the rows `x` (one per input and arm)
pair <- function(x, a, b) {
  A <- x[x$arm == a, ]; Bm <- x[x$arm == b, ]
  m <- merge(A, Bm, by = c("input", "unit", "set", "half", "true_mult"), suffixes = c("_a", "_b"))
  m[is.finite(m$w1_a) & is.finite(m$w1_b), ]
}
# arm names that depend on each input's true rate
arm_for <- function(x, pattern) if (grepl("%d", pattern)) sprintf(pattern, x$true_mult) else rep(pattern, nrow(x))
pair_t <- function(x, pa, pb) {
  x <- x; x$arm_a <- arm_for(x, pa); x$arm_b <- arm_for(x, pb)
  A <- x[x$arm == x$arm_a, ]; Bm <- x[x$arm == x$arm_b, ]
  m <- merge(A, Bm, by = c("input", "unit", "set", "half", "true_mult"), suffixes = c("_a", "_b"))
  m[is.finite(m$w1_a) & is.finite(m$w1_b), ]
}
stat_boot <- function(m, fun) {
  est <- fun(m)
  ids <- unique(m$unit); rows <- split(seq_len(nrow(m)), factor(m$unit, levels = ids))
  set.seed(seed)
  b <- replicate(B, fun(m[unlist(rows[sample.int(length(ids), replace = TRUE)]), , drop = FALSE]))
  c(est = est, low = unname(stats::quantile(b, 0.025, na.rm = TRUE)), high = unname(stats::quantile(b, 0.975, na.rm = TRUE)), n = nrow(m), sims = length(ids))
}
w1r <- function(m) stats::median(m$w1_a / m$w1_b)
angd <- function(m) stats::median(m$angle_a - m$angle_b, na.rm = TRUE)
newr <- function(m) { ok <- m$new_err_b > 1e-9; stats::median(m$new_err_a[ok] / m$new_err_b[ok]) }
halr <- function(m) { ok <- m$halo_err_b > 1e-9; stats::median(m$halo_err_a[ok] / m$halo_err_b[ok]) }
summ <- function(m, label) {
  data.frame(comparison = label, inputs = nrow(m), sims = length(unique(m$unit)),
             t(stat_boot(m, w1r)[c("est", "low", "high")]), angle_change = angd(m), new_err_ratio = newr(m), halo_err_ratio = halr(m),
             check.names = FALSE)
}
vm <- L[L$source == "misseg", ]

## Q1 --------------------------------------------------------------------------------------
hi <- vm[vm$true_mult >= 4, ]
fo <- pair_t(hi, "f1_g%d", "f1_g1")
# "f%d_g%d" needs the rate twice
hi$arm_both <- sprintf("f%d_g%d", hi$true_mult, hi$true_mult)
pair_both <- function(x, other) {
  A <- x[x$arm == x$arm_both, ]; Bm <- x[x$arm == other(x), ]
  m <- merge(A, Bm, by = c("input", "unit", "set", "half", "true_mult"), suffixes = c("_a", "_b"))
  m[is.finite(m$w1_a) & is.finite(m$w1_b), ]
}
both_fo <- pair_both(hi, function(x) sprintf("f1_g%d", x$true_mult))
both <- pair_both(hi, function(x) rep("f1_g1", nrow(x)))
q1 <- rbind(summ(fo, "FO (f1, g true) vs base (f1, g1)"), summ(both_fo, "BOTH (f = g = true) vs FO"), summ(both, "BOTH vs base"))
cat("\n== Q1 (varyMisseg, true rate >= x4; W1 = median paired ratio, 95% interval by simulation)\n"); print(q1, row.names = FALSE, digits = 3)
m_fo <- q1$est[1]; m_b <- q1$est[2]; m_bb <- q1$est[3]
rec <- if (m_fo <= 0.95 && m_b >= 0.98) "forecast only" else if (m_bb <= 0.95 && m_b < 0.98) "fit and forecast" else
  if (m_fo > 0.95 && m_bb > 0.95) "a correct rate does not improve W1" else "forecast only (fit adds < 2%)"
cat("Q1 recommendation:", rec, "\n")
utils::write.table(q1, file.path(res_dir, "q1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Q2: by true rate ----------------------------------------------------------------------------
q2 <- list()
for (s in sort(unique(L$set))) {
  x <- L[L$set == s, ]; t <- x$true_mult[1]
  arms <- c(FO = sprintf("f1_g%d", t), BOTH = if (t > 1) sprintf("f%d_g%d", t, t) else NA, truth_g1 = "truth_g1",
            truth_gtrue = if (t > 1) "truth_gtrue" else NA, mix_1_4 = "mix_1_4", mix_1_8 = "mix_1_8", mix_1_16 = "mix_1_16")
  for (nm in names(arms)) {
    if (is.na(arms[[nm]]) || !arms[[nm]] %in% x$arm) next
    m <- pair(x, arms[[nm]], "f1_g1")
    up <- x[x$arm == arms[[nm]] & x$obs_halo_change > 0, ]
    q2[[length(q2) + 1L]] <- cbind(set = s, arm = nm, summ(m, paste(arms[[nm]], "vs f1_g1"))[, -1],
                                   direction_right = mean(up$pred_halo_change > 0),
                                   base_direction_right = mean(x$pred_halo_change[x$arm == "f1_g1" & x$obs_halo_change > 0] > 0))
  }
}
q2 <- do.call(rbind, q2)
cat("\n== Q2: each arm against base (f1_g1), by true rate\n"); print(q2, row.names = FALSE, digits = 3)
utils::write.table(q2, file.path(res_dir, "q2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== medians of W1, angle and new-karyotype share by set and arm\n")
print(reshape(stats::aggregate(cbind(w1, angle, pred_new, obs_new) ~ set + arm, L[L$arm %in% c("f1_g1", "f1_g4", "f1_g16", "mix_1_16", "truth_g1", "truth_gtrue"), ], stats::median),
              idvar = "set", timevar = "arm", direction = "wide"), row.names = FALSE, digits = 3)

## Q3: misspecification --------------------------------------------------------------------
grid <- c(1, 2, 4, 8, 16)
q3 <- list()
for (s in sort(unique(L$set))) {
  x <- L[L$set == s, ]; t <- x$true_mult[1]
  for (g in setdiff(grid, t)) {
    m <- pair(x, sprintf("f1_g%d", g), sprintf("f1_g%d", t))
    if (nrow(m)) q3[[length(q3) + 1L]] <- data.frame(set = s, usage = "forecast only", given = g, true = t, factor = g / t,
      inputs = nrow(m), w1_ratio = w1r(m), angle_change = angd(m), new_err_ratio = newr(m))
    if (s != "main") {
      a <- if (g == 1) "f1_g1" else sprintf("f%d_g%d", g, g); b <- if (t == 1) "f1_g1" else sprintf("f%d_g%d", t, t)
      m <- pair(x, a, b)
      if (nrow(m)) q3[[length(q3) + 1L]] <- data.frame(set = s, usage = "fit and forecast", given = g, true = t, factor = g / t,
        inputs = nrow(m), w1_ratio = w1r(m), angle_change = angd(m), new_err_ratio = newr(m))
    }
  }
}
q3 <- do.call(rbind, q3)
cat("\n== Q3: given rate against the true rate (W1 median ratio, angle change)\n"); print(q3, row.names = FALSE, digits = 3)
utils::write.table(q3, file.path(res_dir, "q3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Q4: mixture --------------------------------------------------------------------------------
mixes <- c("mix_1_4", "mix_1_8", "mix_1_16")
q4 <- function(h) do.call(rbind, lapply(mixes, function(mx) {
  y <- L[L$half == h, ]
  gain <- pair(y[y$source == "misseg" & y$true_mult >= 4, ], mx, "f1_g1")
  c_main <- pair(y[y$set == "main", ], mx, "f1_g1"); c_vm <- pair(y[y$set == "vM x1", ], mx, "f1_g1")
  sg <- stat_boot(gain, w1r)
  data.frame(half = h, mixture = mx, gain_w1 = sg[["est"]], gain_low = sg[["low"]], gain_high = sg[["high"]],
             cost_main = w1r(c_main), cost_vm1 = w1r(c_vm), angle_main = angd(c_main), angle_vm1 = angd(c_vm),
             gain_angle = angd(gain), gain_new_err = newr(gain))
}))
dev4 <- q4("dev")
cat("\n== Q4 development half\n"); print(dev4, row.names = FALSE, digits = 3)
okd <- dev4[dev4$cost_main <= 1.02 & dev4$cost_vm1 <= 1.02, ]
if (!nrow(okd)) {
  q4_verdict <- "Q4 FAILED on the development half: no range costs <= 1.02 at x1"
} else {
  chosen <- okd$mixture[which.min(okd$gain_w1)]
  con <- q4("confirm"); cat("\n== Q4 confirmation half\n"); print(con, row.names = FALSE, digits = 3)
  r <- con[con$mixture == chosen, ]
  pass <- r$cost_main <= 1.02 && r$cost_vm1 <= 1.02 && r$gain_w1 <= 0.95 && r$gain_high < 1 && r$angle_main <= 1 && r$angle_vm1 <= 1
  q4_verdict <- sprintf("Q4 %s: chosen %s on development; confirmation gain %.3f (%.3f, %.3f), cost main %.3f, cost vM x1 %.3f, angle change at x1 %+.2f / %+.2f deg",
                        if (pass) "PASSED" else "FAILED", chosen, r$gain_w1, r$gain_low, r$gain_high, r$cost_main, r$cost_vm1, r$angle_main, r$angle_vm1)
  utils::write.table(con, file.path(res_dir, "q4_confirm.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
utils::write.table(dev4, file.path(res_dir, "q4_dev.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n", q4_verdict, "\n")

## KARST against ALFA-K (added with the ALFA-K arms; descriptive) ---------------------------
if (any(L$arm == "alfak_g1")) {
  vs <- function(m, label) {
    sb <- stat_boot(m, w1r)
    data.frame(comparison = label, inputs = nrow(m), sims = length(unique(m$unit)), w1_ratio = sb[["est"]], low = sb[["low"]], high = sb[["high"]],
               karst_w1_smaller = mean(m$w1_a < m$w1_b), angle_change = angd(m), karst_angle_smaller = mean(m$angle_a < m$angle_b, na.rm = TRUE),
               new_err_ratio = newr(m), halo_err_ratio = halr(m))
  }
  ka <- list()
  for (s in sort(unique(L$set))) {
    x <- L[L$set == s, ]; t <- x$true_mult[1]
    cmp <- list(c("f1_g1", "alfak_g1", "KARST default vs ALFA-K default"), c("truth_g1", "alfak_g1", "true landscape x1 vs ALFA-K default"))
    if (t > 1) cmp <- c(cmp, list(c(sprintf("f1_g%d", t), "alfak_gtrue", "both forecast at the true rate (fits at 5e-5)"),
                                  c(sprintf("f%d_g%d", t, t), "alfakT_gtrue", "both fitted and forecast at the true rate"),
                                  c(sprintf("f%d_g%d", t, t), "alfak_g1", "KARST at the true rate vs ALFA-K default")))
    for (cc in cmp) { m <- pair(x, cc[1], cc[2]); if (nrow(m)) ka[[length(ka) + 1L]] <- cbind(set = s, vs(m, cc[3])) }
  }
  ka <- do.call(rbind, ka)
  cat("\n== KARST (a) against ALFA-K (b): W1 = median paired ratio a/b (95% interval by simulation); shares of inputs where KARST is better\n")
  print(ka, row.names = FALSE, digits = 3)
  utils::write.table(ka, file.path(res_dir, "karst_vs_alfak.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n== medians by set: KARST and ALFA-K arms\n")
  print(reshape(stats::aggregate(cbind(w1, angle) ~ set + arm, L[L$arm %in% c("f1_g1", "alfak_g1", "alfak_gtrue", "alfakT_gtrue", "truth_g1"), ], stats::median),
                idvar = "set", timevar = "arm", direction = "wide"), row.names = FALSE, digits = 3)
  up <- L[L$obs_halo_change > 0 & L$arm %in% c("f1_g1", "alfak_g1", "alfak_gtrue", "alfakT_gtrue"), ]
  cat("\n== near out-of-domain share rose: share predicting a rise, by set and arm\n")
  print(reshape(stats::aggregate(cbind(rise = pred_halo_change > 0) ~ set + arm, up, mean), idvar = "set", timevar = "arm", direction = "wide"), row.names = FALSE, digits = 3)
}

## Salehi --------------------------------------------------------------------------------------
S <- list(); rep_w <- c()
for (f in list.files(file.path(out_dir, "salehi"), "[.]rds$", full.names = TRUE)) {
  x <- readRDS(f); if (!identical(x$status, "ok")) { cat("Salehi", basename(f), x$status, "\n"); next }
  rep_w[x$task$window] <- x$reproduce_max_rel
  S[[length(S) + 1L]] <- x$rows
}
S <- do.call(rbind, S)
cat("\n== Salehi: x1 fit reproduces codex's karst_native (<= 1e-6) in", sum(rep_w <= 1e-6), "of", length(rep_w), "windows\n")
if (any(rep_w > 1e-6)) print(rep_w[rep_w > 1e-6])
S <- S[S$window %in% names(rep_w)[rep_w <= 1e-6], ]
S$family <- ifelse(grepl("^SA906", S$window), "SA906", ifelse(grepl("^SA609|^SA000", S$window), "SA609", "other"))
utils::write.table(S, file.path(res_dir, "salehi_rows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("daughters:", length(unique(paste(S$window, S$daughter))), "\n")
med <- stats::aggregate(cbind(angle, w1, overlap, new_err, pred_new, obs_new) ~ arm, S, stats::median)
up <- S[S$obs_halo_change > 0, ]
med$direction_right <- tapply(up$pred_halo_change > 0, up$arm, mean)[med$arm]
cat("\n== Salehi medians over daughters, by arm (direction: share predicting a rise where the near out-of-domain share rose)\n"); print(med, row.names = FALSE, digits = 3)
for (fam in c("SA906", "SA609", "other")) {
  cat("\n--", fam, "\n"); print(stats::aggregate(cbind(angle, w1, overlap, new_err) ~ arm, S[S$family == fam, ], stats::median), row.names = FALSE, digits = 3)
}
al <- S[S$arm == "alfak_g1", ]
wins <- do.call(rbind, lapply(setdiff(unique(S$arm), "alfak_g1"), function(a) {
  x <- S[S$arm == a, ]; mm <- match(paste(x$window, x$daughter), paste(al$window, al$daughter))
  data.frame(arm = a, angle = sum(x$angle < al$angle[mm]), w1 = sum(x$w1 < al$w1[mm], na.rm = TRUE), overlap = sum(x$overlap > al$overlap[mm]))
}))
cat("\n== Salehi wins against ALFA-K (x1), of", nrow(al), "daughters\n"); print(wins, row.names = FALSE)
writeLines(c(paste("Q1:", rec), q4_verdict), file.path(res_dir, "verdict.txt"))
