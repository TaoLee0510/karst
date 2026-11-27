# What the CV score (R2R) of each method measures, on the ABM fits.
#
#   Rscript scripts/cv_score_analysis.R
#
# Writes BENCH_RESULTS/analysis/cv_fits.Rds (one row per fit) and prints the
# tables. Needs karst on the library path (fold membership).
#
# R2R = 1 - sum((p' - o')^2) / sum(o'^2), p', o' centred
#     = 2 r s - s^2,  r = cor(o, p), s = sd(p) / sd(o)
# so R2R < 0 exactly when s > 2 r, and R2R <= r^2.
# Per fit: this decomposition, the share of the error that is a common offset
# per fold, and the same quantities against the true fitness of the held-out
# karyotypes. Per input: true landscape accuracy of both methods.
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
options(width = 220)
neighbours <- karst:::alfak_neighbours

centred_stats <- function(target, x) {
  tc <- target - mean(target)
  xc <- x - mean(x)
  c(r = cor(target, x), s = sqrt(sum(xc^2) / sum(tc^2)), R2R = 1 - sum((xc - tc)^2) / sum(tc^2),
    slope = sum(xc * tc) / sum(xc^2))  # regression of target on x; < 1: x over-dispersed
}

fit_rows <- function(sim_id) {
  sim <- read_abm_sim(sim_id, paths$abm_raw)
  w <- abm_sim_info(sim_id)$w
  dirs <- Sys.glob(file.path(paths$results, "fits", "abm", sim_id, "*", "*", "minobs_*"))
  rows <- lapply(dirs, function(d) {
    if (!file.exists(file.path(d, "score.Rds"))) return(NULL)
    p_ <- rev(strsplit(d, "/", fixed = TRUE)[[1]])
    x <- readRDS(file.path(d, "xval.Rds"))
    o <- unname(x$observations)
    p <- unname(x$predictions)
    k <- x$xval_data$k
    ok <- is.finite(o) & is.finite(p)
    sc <- readRDS(file.path(d, "score.Rds"))$summary
    base <- data.frame(sim_id = sim_id, w = w, design = p_[3], method = p_[2], minobs = sub("minobs_", "", p_[1]),
                       Rxv = sc$Rxv, r_true = sc$r, R_true = sc$R, Rfq_true = sc$Rfq, nfq = sc$nfq,
                       n = sum(ok), stringsAsFactors = FALSE)
    if (sum(ok) < 3) return(base)
    o <- o[ok]; p <- p[ok]; k <- k[ok]
    truth <- grf_fitness(k, w, sim$true_landscape)
    po <- centred_stats(o, p)     # CV score as computed: prediction vs the method's own full-data value
    pt <- centred_stats(truth, p) # held-out prediction vs truth
    ot <- centred_stats(truth, o) # full-data value vs truth
    # folds as in both packages: each frequent karyotype with its neighbours, first fold wins
    fq <- colnames(readRDS(file.path(d, "bootstrap_res.Rds"))$final_fitness)
    fold <- unlist(lapply(seq_along(fq), function(i) {
      m <- c(fq[i], neighbours(fq[i]))
      stats::setNames(rep(i, length(m)), m)
    }))
    fold <- fold[!duplicated(names(fold))]
    fid <- fold[k]
    e <- (p - mean(p)) - (o - mean(o))
    between <- if (all(is.na(fid))) NA else {
      keep <- !is.na(fid)
      sum(tapply(e[keep], fid[keep], function(v) length(v) * mean(v)^2)) / sum(e[keep]^2)
    }
    cbind(base, r = po[["r"]], s = po[["s"]], R2R = po[["R2R"]],
          r_pt = pt[["r"]], s_pt = pt[["s"]], R2R_pt = pt[["R2R"]], slope_pt = pt[["slope"]],
          r_ot = ot[["r"]], s_ot = ot[["s"]], R2R_ot = ot[["R2R"]], slope_ot = ot[["slope"]],
          share_fq = mean(k %in% fq), bias_fq = mean(e[k %in% fq]) / sd(o), n_folds = length(unique(fid[!is.na(fid)])), between = between,
          sd_obs = sd(o), sd_pred = sd(p), sd_true = sd(truth))
  })
  bind(rows)
}
bind <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
}

sims <- sort(list.files(file.path(paths$results, "fits", "abm")))
res <- bind(lapply(sims, function(s) tryCatch(fit_rows(s), error = function(e) {
  message(s, ": ", conditionMessage(e)); NULL
})))
dir.create(file.path(paths$results, "analysis"), showWarnings = FALSE)
saveRDS(res, file.path(paths$results, "analysis", "cv_fits.Rds"))
med <- function(x) round(median(x, na.rm = TRUE), 3)
q <- function(x, p) round(unname(quantile(x, p, na.rm = TRUE)), 2)
cat("fits:", nrow(res), " with < 3 pairs:", sum(res$n < 3), "\n")

cat("\n== A. karst, all wavelengths: CV score and its parts, by design and wavelength\n")
x2 <- res[res$method == "karst" & res$n >= 3, ]
a <- do.call(rbind, lapply(split(x2, list(x2$design, x2$w), drop = TRUE), function(d) data.frame(
  design = d$design[1], w = d$w[1], fits = nrow(d), folds = med(d$n_folds),
  R2R_q10 = q(d$R2R, .1), R2R_q25 = q(d$R2R, .25), R2R_med = med(d$R2R), R2R_q75 = q(d$R2R, .75), min = q(d$R2R, 0),
  lt0 = round(mean(d$R2R < 0), 2), r = med(d$r), s = med(d$s), between = med(d$between), bias_fq = med(d$bias_fq),
  heldout_r_true = med(d$r_pt), heldout_slope = med(d$slope_pt), full_r_true = med(d$r_ot), full_s_true = med(d$s_ot),
  r_true = med(d$r_true))))
print(a[order(a$design, a$w), ], row.names = FALSE)

cat("\n== B. Paired inputs (same simulation, design, minobs; both fits ok), by design and wavelength\n")
key <- c("sim_id", "design", "minobs")
both <- merge(res[res$method == "karst", key], res[res$method == "alfak", key], by = key)
pr <- merge(res, both, by = key)
pr <- pr[pr$n >= 3, ]
b <- do.call(rbind, lapply(split(pr, list(pr$method, pr$design, pr$w), drop = TRUE), function(d) data.frame(
  design = d$design[1], w = d$w[1], method = d$method[1], fits = nrow(d), pairs = med(d$n),
  R2R_med = med(d$R2R), R2R_q10 = q(d$R2R, .1), lt0 = round(mean(d$R2R < 0), 2), r = med(d$r), s = med(d$s),
  between = med(d$between), sd_obs_vs_true = med(d$s_ot),
  heldout_r_true = med(d$r_pt), full_r_true = med(d$r_ot), landscape_r_true = med(d$r_true), landscape_R_true = med(d$R_true),
  pass_xv = round(mean(d$Rxv > 0, na.rm = TRUE), 2),
  rho_xv_rtrue = round(cor(d$Rxv, d$r_true, method = "spearman", use = "complete"), 2),
  rtrue_xv_pos = med(d$r_true[d$Rxv > 0]), rtrue_xv_neg = med(d$r_true[d$Rxv <= 0]))))
print(b[order(b$design, b$w, b$method), ], row.names = FALSE)

cat("\n== C. Paired wins of karst over alfak on the same input\n")
wide <- merge(pr[pr$method == "karst", c(key, "w", "r_true", "R_true", "r_pt", "Rxv")],
              pr[pr$method == "alfak", c(key, "r_true", "R_true", "r_pt", "Rxv")], by = key, suffixes = c(".a2", ".aR"))
cw <- do.call(rbind, lapply(split(wide, list(wide$design, wide$w), drop = TRUE), function(d) data.frame(
  design = d$design[1], w = d$w[1], inputs = nrow(d),
  landscape_r = round(mean(d$r_true.a2 > d$r_true.aR), 2), landscape_R = round(mean(d$R_true.a2 > d$R_true.aR), 2),
  heldout_r_true = round(mean(d$r_pt.a2 > d$r_pt.aR, na.rm = TRUE), 2), cv_score = round(mean(d$Rxv.a2 > d$Rxv.aR, na.rm = TRUE), 2))))
print(cw[order(cw$design, cw$w), ], row.names = FALSE)

cat("\n== D. karst fits with R2R < -3 vs the rest (all wavelengths)\n")
ext <- x2$R2R < -3
cat("extreme:", sum(ext), "of", nrow(x2), "\n")
print(table(design = x2$design[ext], w = x2$w[ext]))
cols <- c("nfq", "n", "n_folds", "r", "s", "between", "bias_fq", "sd_obs", "sd_pred", "r_pt", "slope_pt", "r_ot", "r_true", "R_true")
print(rbind(extreme = sapply(x2[ext, cols], med), rest = sapply(x2[!ext, cols], med)))
cat("\nShare of extremes by number of folds:\n")
x2$fold_bin <- cut(x2$n_folds, c(0, 2, 4, 8, 16, Inf))
print(round(tapply(x2$R2R < -3, x2$fold_bin, mean), 2))
print(table(x2$fold_bin))
