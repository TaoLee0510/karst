# Part 2 of PLAN.md: summary of the own benchmark at missegregation multiples x1-x64.
# Same aggregation as the round report (doc/report/build/own_benchmark_summary.R):
# mean over designs within a landscape, then over landscapes; error ratio = mean RMSE /
# mean null RMSE. Main table distances 0-2. Core designs for KARST against ALFA-K
# (ALFA-K was fitted on core only), all designs for KARST alone. Bootstrap: landscapes
# resampled within wavelength (2000 draws, seed 1).
#   Rscript own_grid_summary.R TASK_FILE ID OUT_DIR   (reads OUT_DIR/{metrics,direction}_heldout_grid.csv)
args <- commandArgs(TRUE)
out <- args[[3]]
sdir <- file.path(out, "summary"); dir.create(sdir, showWarnings = FALSE)
m <- utils::read.csv(file.path(out, "metrics_heldout_grid.csv"), stringsAsFactors = FALSE)
dr <- utils::read.csv(file.path(out, "direction_heldout_grid.csv"), stringsAsFactors = FALSE)
idx <- unique(m[c("dataset_id", "landscape_id", "wavelength_class", "scenario", "role")])
dr <- merge(dr, idx, by = "dataset_id")
mults <- c(1, 2, 4, 8, 16, 32, 64)
lab <- c(karst_sparse_joint_balance = "KARST", alfak_minobs10 = "ALFA-K m10", zero = "zero")
waves <- c("w0.2", "w0.4", "w0.8", "w1.2", "w1.6")
dists <- c("d0", "d1", "d2")
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
f3 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
md <- function(x) c(paste("|", paste(names(x), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
                    apply(x, 1, function(r) paste("|", paste(trimws(r), collapse = " | "), "|")))
L <- c("# Own benchmark at missegregation multiples (KARST 1.1.4 rc, held-out main table)", "",
       sprintf("Generated %s from `metrics_heldout_grid.csv` and `direction_heldout_grid.csv`.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "")

## completeness -----------------------------------------------------------------------------------
fits <- m[!duplicated(m[c("dataset_id", "method", "multiple")]), ]
fits$n <- 1; fits$failed <- as.integer(fits$status != "ok")
comp <- merge(stats::aggregate(cbind(n, failed) ~ method + multiple, fits, sum),
              stats::aggregate(fits["seconds"], fits[c("method", "multiple")], function(v) stats::median(v, na.rm = TRUE)),
              by = c("method", "multiple"))
names(comp)[names(comp) == "seconds"] <- "median_seconds"
comp$method <- lab[comp$method]; comp <- comp[order(comp$method, comp$multiple), ]
utils::write.table(comp, file.path(sdir, "completeness.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
comp$median_seconds <- formatC(comp$median_seconds, format = "f", digits = 1)
L <- c(L, "## Fits scored (failed fits are scored as the zero prediction)", "", md(comp), "")

## landscape-level tables -------------------------------------------------------------------------
m <- m[m$stratum %in% dists, ]
land_of <- function(x, cols) {
  cell <- stats::aggregate(x[cols], x[c("method", "multiple", "wavelength_class", "landscape_id", "scenario", "stratum")], mean_na)
  stats::aggregate(cell[cols], cell[c("method", "multiple", "wavelength_class", "landscape_id", "stratum")], mean_na)
}
cols <- c("rmse", "rmse_zero", "spearman", "coverage95")
land_core <- land_of(m[m$role == "core", ], cols)
land_all <- land_of(m, cols)
dcell <- stats::aggregate(direction ~ method + multiple + wavelength_class + landscape_id + scenario + stratum, dr, mean_na)
dl_core <- stats::aggregate(direction ~ method + multiple + wavelength_class + landscape_id + stratum, dcell[dcell$scenario %in% unique(m$scenario[m$role == "core"]), ], mean_na)
dl_all <- stats::aggregate(direction ~ method + multiple + wavelength_class + landscape_id + stratum, dcell, mean_na)

# value of a statistic per (method, multiple, wave or pooled, stratum)
stat <- function(land, v, ratio = FALSE) {
  key <- c("method", "multiple", "stratum")
  per_w <- if (ratio) {
    a <- stats::aggregate(land[c("rmse", "rmse_zero")], land[c(key, "wavelength_class")], mean_na); a$value <- a$rmse / a$rmse_zero; a
  } else { a <- stats::aggregate(land[v], land[c(key, "wavelength_class")], mean_na); a$value <- a[[v]]; a }
  pooled <- if (ratio) {
    a <- stats::aggregate(land[c("rmse", "rmse_zero")], land[key], mean_na); a$value <- a$rmse / a$rmse_zero; a$wavelength_class <- "all"; a
  } else { a <- stats::aggregate(land[v], land[key], mean_na); a$value <- a[[v]]; a$wavelength_class <- "all"; a }
  rbind(per_w[c(key, "wavelength_class", "value")], pooled[c(key, "wavelength_class", "value")])
}
tab_mult <- function(s, method, wave = "all") {
  x <- s[s$method == method & s$wavelength_class == wave, ]
  t <- data.frame(multiple = paste0("x", mults))
  for (d in dists) t[[sub("d", "distance ", d)]] <- f3(vapply(mults, function(k) { v <- x$value[x$multiple == k & x$stratum == d]; if (length(v)) v else NA_real_ }, 0))
  t
}
sets <- list(
  list(name = "error ratio (RMSE / null RMSE; lower is better)", file = "ratio", land = land_core, v = "rmse", ratio = TRUE, dland = NULL),
  list(name = "direction accuracy of one-copy changes (ties 1/2; higher is better)", file = "direction", land = dl_core, v = "direction", ratio = FALSE),
  list(name = "Spearman correlation with the truth (higher is better)", file = "spearman", land = land_core, v = "spearman", ratio = FALSE))
L <- c(L, "## Core designs, all wavelengths pooled (500 landscapes weighted equally)", "")
for (st in sets) {
  s <- stat(st$land, st$v, st$ratio)
  utils::write.table(s, file.path(sdir, paste0("core_", st$file, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  both <- merge(tab_mult(s, "karst_sparse_joint_balance"), tab_mult(s, "alfak_minobs10"), by = "multiple", suffixes = c(" KARST", " ALFA-K"))
  both <- both[match(paste0("x", mults), both$multiple), ]
  L <- c(L, sprintf("**%s**", st$name), "", md(both), "")
  if (st$file == "ratio") {
    z <- s[s$method == "zero" & s$wavelength_class == "all", ]
    L <- c(L, sprintf("Zero prediction: %s (by definition 1).", paste(f3(z$value), collapse = ", ")), "")
  }
}
cv <- stat(land_core[land_core$method == "karst_sparse_joint_balance", ], "coverage95")
utils::write.table(cv, file.path(sdir, "core_coverage95.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
L <- c(L, "**KARST coverage of the 95% interval (model sd only; ALFA-K reports none)**", "", md(tab_mult(cv, "karst_sparse_joint_balance")), "")

## per wavelength ---------------------------------------------------------------------------------
L <- c(L, "## Core designs by wavelength", "")
for (st in sets[1:2]) {
  s <- stat(st$land, st$v, st$ratio)
  for (meth in c("karst_sparse_joint_balance", "alfak_minobs10")) {
    t <- do.call(rbind, lapply(waves, function(w) cbind(wavelength = sub("w", "", w), tab_mult(s, meth, w))))
    wide <- stats::reshape(t, idvar = "wavelength", timevar = "multiple", direction = "wide")
    for (d in dists) {
      dd <- sub("d", "distance ", d)
      tt <- wide[c("wavelength", paste0(dd, ".x", mults))]; names(tt) <- c("λ", paste0("x", mults))
      L <- c(L, sprintf("**%s, %s, %s**", lab[[meth]], st$name, dd), "", md(tt), "")
    }
  }
}

## paired bootstrap -------------------------------------------------------------------------------
set.seed(1)
B <- 2000
boot_diff <- function(a, b, value, ratio) {
  # a, b: landscape tables (one row per landscape) for the same landscapes; returns mean, 2.5%, 97.5%
  k <- merge(a, b, by = c("wavelength_class", "landscape_id"), suffixes = c(".a", ".b"))
  k <- k[stats::complete.cases(k[grep("\\.(a|b)$", names(k))]), ]
  est <- function(r) if (ratio) (mean(r$rmse.a) - mean(r$rmse.b)) / mean(r$rmse_zero.a) else mean(r[[paste0(value, ".a")]]) - mean(r[[paste0(value, ".b")]])
  groups <- split(seq_len(nrow(k)), k$wavelength_class)
  draws <- replicate(B, est(k[unlist(lapply(groups, function(g) g[sample.int(length(g), length(g), TRUE)])), ]))
  win <- if (ratio) mean(k$rmse.a < k$rmse.b) else mean(k[[paste0(value, ".a")]] > k[[paste0(value, ".b")]])
  c(diff = est(k), lo = unname(stats::quantile(draws, 0.025)), hi = unname(stats::quantile(draws, 0.975)), share_better = win, n = nrow(k))
}
pick <- function(land, meth, k, d, v) { x <- land[land$method == meth & land$multiple == k & land$stratum == d, ]; x[c("wavelength_class", "landscape_id", intersect(c("rmse", "rmse_zero", v), names(x)))] }
pb <- list()
for (d in dists) for (k in mults) {
  pb[[length(pb) + 1]] <- data.frame(comparison = "KARST - ALFA-K m10 (core)", metric = "error ratio", multiple = k, stratum = d,
    t(boot_diff(pick(land_core, "karst_sparse_joint_balance", k, d, "rmse"), pick(land_core, "alfak_minobs10", k, d, "rmse"), "rmse", TRUE)))
  pb[[length(pb) + 1]] <- data.frame(comparison = "KARST - ALFA-K m10 (core)", metric = "direction", multiple = k, stratum = d,
    t(boot_diff(pick(dl_core, "karst_sparse_joint_balance", k, d, "direction"), pick(dl_core, "alfak_minobs10", k, d, "direction"), "direction", FALSE)))
  if (k > 1) {
    pb[[length(pb) + 1]] <- data.frame(comparison = "KARST xk - KARST x1 (all designs)", metric = "error ratio", multiple = k, stratum = d,
      t(boot_diff(pick(land_all, "karst_sparse_joint_balance", k, d, "rmse"), pick(land_all, "karst_sparse_joint_balance", 1, d, "rmse"), "rmse", TRUE)))
    pb[[length(pb) + 1]] <- data.frame(comparison = "KARST xk - KARST x1 (all designs)", metric = "direction", multiple = k, stratum = d,
      t(boot_diff(pick(dl_all, "karst_sparse_joint_balance", k, d, "direction"), pick(dl_all, "karst_sparse_joint_balance", 1, d, "direction"), "direction", FALSE)))
    pb[[length(pb) + 1]] <- data.frame(comparison = "ALFA-K xk - ALFA-K x1 (core)", metric = "error ratio", multiple = k, stratum = d,
      t(boot_diff(pick(land_core, "alfak_minobs10", k, d, "rmse"), pick(land_core, "alfak_minobs10", 1, d, "rmse"), "rmse", TRUE)))
  }
}
pb <- do.call(rbind, pb)
utils::write.table(pb, file.path(sdir, "paired_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
fmt_pb <- function(comp, metric) {
  x <- pb[pb$comparison == comp & pb$metric == metric, ]
  t <- data.frame(multiple = paste0("x", sort(unique(x$multiple))))
  for (d in dists) t[[sub("d", "distance ", d)]] <- vapply(sort(unique(x$multiple)), function(k) {
    r <- x[x$multiple == k & x$stratum == d, ]
    sprintf("%+.3f [%+.3f, %+.3f] (%.0f%%)", r$diff, r$lo, r$hi, 100 * r$share_better)
  }, "")
  t
}
L <- c(L, "## Paired differences, all wavelengths (landscape bootstrap 95% interval; share of landscapes where the first is better)", "",
       "A negative error-ratio difference and a positive direction difference favour the first.", "")
for (cm in list(c("KARST - ALFA-K m10 (core)", "error ratio"), c("KARST - ALFA-K m10 (core)", "direction"),
                c("KARST xk - KARST x1 (all designs)", "error ratio"), c("KARST xk - KARST x1 (all designs)", "direction"),
                c("ALFA-K xk - ALFA-K x1 (core)", "error ratio")))
  L <- c(L, sprintf("**%s, %s**", cm[1], cm[2]), "", md(fmt_pb(cm[1], cm[2])), "")

## KARST on all 14 designs, and stress designs alone --------------------------------------------------
sa <- stat(land_all[land_all$method == "karst_sparse_joint_balance", ], "rmse", TRUE)
ss <- stat(land_of(m[m$role != "core" & m$method == "karst_sparse_joint_balance", ], cols), "rmse", TRUE)
L <- c(L, "## KARST error ratio on all 14 designs and on the 6 stress designs (all wavelengths)", "",
       "**All designs**", "", md(tab_mult(sa, "karst_sparse_joint_balance")), "",
       "**Stress designs**", "", md(tab_mult(ss, "karst_sparse_joint_balance")), "")
writeLines(L, file.path(sdir, "own_grid_summary.md"))
cat(L, sep = "\n")
