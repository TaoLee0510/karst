#!/usr/bin/env Rscript
# Evaluation of dev/realdata/noise/PLAN_v4.md: the four criteria for T1 and T2
# on the confirmation half (simulations) and the Salehi confirm group, the
# choice, and the report-only tables (classes, phi_hat against the injected
# phi, the higher_beta class, the development half, NUMBAT, parallel lineages).
#   Rscript dev/realdata/noise/v4/v4_eval.R [dev]   ("dev": the development half instead, report only)
# Writes realdata/noise/v4/evaluation[_dev].rds.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
half <- if (length(commandArgs(TRUE)) && commandArgs(TRUE)[[1]] == "dev") "dev" else "confirm"
v4 <- file.path(dev_results, "realdata", "noise", "v4")
set.seed(20261120)
read_dir <- function(dir, f) {
  files <- list.files(dir, pattern = "[.][Rr]ds$", recursive = TRUE, full.names = TRUE)
  do.call(rbind, lapply(files, function(p) f(readRDS(p))))
}

## ---------- screen ----------
scr <- read_dir(file.path(v4, "screen"), function(x) {
  if (!x$kind %in% c("eset", "c2")) return(NULL)
  data.frame(kind = x$kind, set = x$set, variant = x$variant, dataset_id = x$dataset_id, landscape_id = x$landscape_id,
             half = x$half, depth = x$depth, bottleneck = x$bottleneck, truth = x$truth, status = x$status,
             class = if (identical(x$status, "ok")) x$class else NA, beta_hat = if (identical(x$status, "ok")) x$beta_hat else NA,
             phi_hat = if (identical(x$status, "ok")) x$phi_hat else NA, stringsAsFactors = FALSE)
})
scr$group <- ifelse(scr$kind == "c2", "C2", scr$set)
cat("screened simulations:", nrow(scr), " failed:", sum(scr$status != "ok"), "\n")
cls <- c("normal", "higher_beta", "noise", "selection")
share_tab <- function(x, by) {
  t <- do.call(rbind, lapply(split(x, x[by], drop = TRUE), function(g)
    cbind(g[1, by, drop = FALSE], n = nrow(g), t(setNames(round(as.numeric(prop.table(table(factor(g$class, levels = cls)))), 3), cls)))))
  t[do.call(order, t[by]), ]
}
cat("\n== classes by set, injection and cells per time point (both halves) ==\n")
print(share_tab(scr[scr$status == "ok", ], c("group", "variant", "depth")), row.names = FALSE)

## ---------- simulation fits ----------
fits <- read_dir(file.path(v4, "fits"), function(x) {
  if (!x$method %in% c("T0", "T1", "T2", "T2r") || !identical(x$library, "karst_v4")) return(NULL)
  m <- if (identical(x$status, "ok")) x$metrics else NULL
  r <- function(s) if (!is.null(m) && any(m$stratum == s)) m$rmse[m$stratum == s] else NA_real_
  data.frame(kind = x$kind, set = x$set, variant = x$variant, dataset_id = x$dataset_id, method = x$method,
             half = x$half, landscape_id = x$landscape_id, class = x$class, depth = x$depth, bottleneck = x$bottleneck,
             truth_beta = x$truth_beta, beta_hat = x$beta_hat, phi_hat = x$phi_hat, status = x$status,
             rmse_d0 = r("d0"), rmse_d1 = r("d1"), rmse_d2 = r("d2"),
             false_selection = if (identical(x$status, "ok")) x$false_selection else NA_real_,
             seconds = if (is.null(x$seconds)) NA_real_ else x$seconds, stringsAsFactors = FALSE)
})
fits$group <- ifelse(fits$kind == "c2", "C2", fits$set)
cat("\nfits:", nrow(fits), "\n"); print(table(fits$method, fits$status))
fits$key <- paste(fits$kind, fits$set, fits$variant, fits$dataset_id)
base <- fits[fits$method == "T0", ]
paired <- function(method) {
  x <- fits[fits$method == method & fits$status == "ok", ]
  b <- base[match(x$key, base$key), ]
  ok <- !is.na(b$status) & b$status == "ok"
  x <- x[ok, ]; b <- b[ok, ]
  rr <- cbind(x$rmse_d0 / b$rmse_d0, x$rmse_d1 / b$rmse_d1, x$rmse_d2 / b$rmse_d2)
  rr[!is.finite(rr)] <- NA
  x$ratio <- rowMeans(rr, na.rm = TRUE)
  x$fs_diff <- x$false_selection - b$false_selection
  x$fs_base <- b$false_selection
  x
}
P <- list(T1 = paired("T1"), T2 = paired("T2"), T2r = paired("T2r"))
med_ci <- function(v, cl, reps = 2000L) {
  ok <- is.finite(v); v <- v[ok]; cl <- cl[ok]
  u <- unique(cl); by <- split(v, cl)
  b <- replicate(reps, stats::median(unlist(by[sample(u, length(u), replace = TRUE)], use.names = FALSE)))
  c(median = stats::median(v), lo = unname(stats::quantile(b, 0.025)), hi = unname(stats::quantile(b, 0.975)), n = length(v))
}
crit <- do.call(rbind, lapply(names(P), function(m) {
  x <- P[[m]][P[[m]]$half == half & P[[m]]$class == "noise", ]
  s1 <- x[x$variant == "none" & x$group %in% c("e1", "e2", "e4", "C2"), ]
  s2 <- x[x$variant == "none" & x$group == "e3", ]
  bn <- x[x$variant %in% c("dm20", "reg10"), ]
  b <- med_ci(bn$ratio, bn$landscape_id)
  data.frame(method = m,
             safe1_median = stats::median(s1$ratio, na.rm = TRUE), safe1_n = sum(is.finite(s1$ratio)),
             safe1 = stats::median(s1$ratio, na.rm = TRUE) <= 1.02,
             safe2_diff = mean(s2$fs_diff, na.rm = TRUE), safe2_n = sum(is.finite(s2$fs_diff)),
             safe2 = mean(s2$fs_diff, na.rm = TRUE) <= 0.02,
             useful_median = b[["median"]], useful_lo = b[["lo"]], useful_hi = b[["hi"]], useful_n = b[["n"]],
             useful = b[["median"]] <= 0.95 && b[["hi"]] < 1)
}))
cat(sprintf("\n== criteria on the %s half (noise class) ==\n", half)); print(crit, digits = 3, row.names = FALSE)

## ---------- Salehi ----------
salehi_a <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
sal_scr <- read_dir(file.path(v4, "screen", "salehi"), function(x)
  data.frame(id = x$dataset_id, line = x$line, group = x$group, n0 = x$n0, n1 = x$n1, status = x$status,
             class = if (identical(x$status, "ok")) x$class else NA, beta_hat = if (identical(x$status, "ok")) x$beta_hat else NA,
             phi_hat = if (identical(x$status, "ok")) x$phi_hat else NA, stringsAsFactors = FALSE))
real <- NULL
if (!is.null(sal_scr)) {
  cat("\n== Salehi classes by line ==\n"); print(table(sal_scr$line, factor(sal_scr$class, levels = cls)))
  hz <- read_dir(file.path(v4, "salehi", "horizon"), function(x) {
    if (!identical(x$status, "ok")) return(NULL)
    r <- x$rows; r$method <- x$task$method; r$id <- x$task$id
    r[c("id", "uid", "days", "method", "ll_fc_s_e2")]
  })
  if (!is.null(hz)) {
    hz <- hz[hz$days <= 60 & is.finite(hz$ll_fc_s_e2), ]
    hz$group <- sal_scr$group[match(hz$id, sal_scr$id)]
    hz$k <- paste(hz$id, hz$uid)
    b0 <- hz[hz$method == "T0", ]
    boot_mean <- function(g, ids, reps = 2000L) {
      u <- unique(ids); by <- split(g, ids)
      b <- replicate(reps, mean(unlist(by[sample(u, length(u), replace = TRUE)], use.names = FALSE)))
      stats::quantile(b, c(0.025, 0.975), names = FALSE)
    }
    real <- do.call(rbind, lapply(c("T1", "T2", "T2r"), function(m) do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
      x <- hz[hz$method == m & hz$group == gr, ]
      x$gain <- x$ll_fc_s_e2 - b0$ll_fc_s_e2[match(x$k, b0$k)]
      x <- x[is.finite(x$gain), ]
      if (!nrow(x)) return(data.frame(method = m, group = gr, passages = 0, trajectories = 0, gain = NA, lo = NA, hi = NA))
      ci <- boot_mean(x$gain, x$id)
      data.frame(method = m, group = gr, passages = nrow(x), trajectories = length(unique(x$id)), gain = mean(x$gain),
                 lo = ci[1], hi = ci[2])
    }))))
    cat("\n== Salehi noise-class trajectories: log score gain over T0 per later passage (<= 60 days) ==\n")
    print(real, digits = 3, row.names = FALSE)
  }
  # parallel lineages: same line and treatment at both ends, disjoint samples; same class?
  sa <- salehi_a[salehi_a$id %in% sal_scr$id, ]
  sa$class <- sal_scr$class[match(sa$id, sal_scr$id)]
  keys <- split(sa, list(sa$line, sa$treat_first, sa$treat_last), drop = TRUE)
  agree <- unlist(lapply(keys, function(g) {
    if (nrow(g) < 2) return(NULL)
    cmb <- utils::combn(nrow(g), 2)
    vapply(seq_len(ncol(cmb)), function(j) {
      a <- g[cmb[1, j], ]; b <- g[cmb[2, j], ]
      if (length(intersect(c(a$first_uid, a$last_uid), c(b$first_uid, b$last_uid)))) return(NA)
      (a$class == "noise") == (b$class == "noise")
    }, logical(1))
  }))
  cat(sprintf("\nparallel lineages with the same noise call: %.3f (%d pairs)\n", mean(agree, na.rm = TRUE), sum(!is.na(agree))))
}

## ---------- verdict ----------
verdict <- crit[c("method", "safe1", "safe2", "useful")]
verdict$real <- vapply(verdict$method, function(m) {
  if (is.null(real)) return(NA)
  x <- real[real$method == m & real$group == "confirm", ]
  if (!nrow(x) || !is.finite(x$gain)) return(NA)
  x$gain >= 0 && x$lo > -0.1
}, logical(1))
verdict$pass <- verdict$safe1 & verdict$safe2 & verdict$useful & verdict$real
cat(sprintf("\n== VERDICT (%s half; Salehi confirm group) ==\n", half)); print(verdict, row.names = FALSE)
passing <- verdict$method[which(verdict$pass)]
choice <- if (!length(passing)) "none" else passing[which.min(crit$useful_median[match(passing, crit$method)])]
cat("CHOICE:", choice, "\n")

## ---------- report only ----------
cat("\n== phi_hat against the injection (noise class, both halves) ==\n")
ns <- scr[scr$status == "ok" & scr$class == "noise", ]
print(stats::aggregate(phi_hat ~ variant + depth, ns, function(v) round(stats::quantile(v, c(0.25, 0.5, 0.75)), 2)))
cat("\n== by set, injection and cells (", half, " half, noise class): median error ratio to T0 ==\n", sep = "")
for (m in names(P)) {
  x <- P[[m]][P[[m]]$half == half & P[[m]]$class == "noise", ]
  if (!nrow(x)) next
  cat(m, "\n"); print(stats::aggregate(ratio ~ group + variant + depth, x, function(v) round(stats::median(v), 3)), row.names = FALSE)
}
hb <- P$T1[P$T1$half == half & P$T1$class == "higher_beta", ]
if (nrow(hb)) {
  cat("\n== higher_beta class, T1 (report only): median error ratio by set and true rate ==\n")
  print(stats::aggregate(ratio ~ group + truth_beta, hb, function(v) c(median = round(stats::median(v), 3), n = length(v))))
}
num <- read_dir(file.path(v4, "screen", "numbat"), function(x)
  data.frame(id = x$dataset_id, cancer = x$cancer, n0 = x$n0, n1 = x$n1, status = x$status,
             class = if (identical(x$status, "ok")) x$class else NA, stringsAsFactors = FALSE))
if (!is.null(num)) { cat("\n== NUMBAT classes ==\n"); print(table(factor(num$class, levels = cls), useNA = "ifany")) }
dev_save_rds(list(half = half, screen = scr, fits = fits, paired = P, criteria = crit, real = real, verdict = verdict,
                  choice = choice, salehi_screen = sal_scr, numbat = num),
             file.path(v4, if (half == "dev") "evaluation_dev.rds" else "evaluation.rds"))
