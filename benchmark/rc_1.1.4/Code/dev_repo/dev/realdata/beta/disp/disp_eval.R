#!/usr/bin/env Rscript
# disp/PLAN.md: the base check against the decomposition study, diagnostics
# D1-D4, the interventions' gains over base F x1 with trajectory-cluster
# bootstrap intervals, and the pre-registered verdicts. Writes
# realdata/beta/disp/evaluation.rds.
#   Rscript dev/realdata/beta/disp/disp_eval.R
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
options(width = 220)
root <- file.path(dev_results, "realdata", "beta")
salehi_a <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
group_of <- function(line) ifelse(line %in% c("SA609", "SA039"), "dev", "confirm")
f_only_loo <- c(dev = 1.0263, confirm = 0.4531)   # decomp/report.md, estimator F
days_bin <- function(d) cut(d, c(0, 15, 30, 45, 60), labels = c("15", "30", "45", "60"))

## ---------- read ----------
files <- list.files(file.path(root, "disp", "forecasts"), pattern = "[.]Rds$", recursive = TRUE, full.names = TRUE)
recs <- lapply(files, readRDS)
status <- vapply(recs, function(x) x$status, character(1))
combo <- vapply(recs, function(x) paste0(x$task$variant, "_F", x$task$F), character(1))
cat("forecast tasks:", length(recs), "\n"); print(table(combo, status))
ok <- recs[status == "ok"]
tag <- function(x, d) { d$variant <- x$task$variant; d$F <- x$task$F; d }
rows <- do.call(rbind, lapply(ok, function(x) tag(x, x$rows)))
kary <- do.call(rbind, lapply(ok, function(x) tag(x, cbind(id = x$task$id, x$kary))))
dmass <- do.call(rbind, lapply(ok, function(x) tag(x, cbind(id = x$task$id, x$dist_mass))))
for (nm in c("rows", "kary", "dmass")) {
  d <- get(nm)
  d$line <- salehi_a$line[match(d$id, salehi_a$id)]
  d$group <- group_of(d$line)
  d <- d[d$days <= 60, ]
  assign(nm, d)
}
rows$key <- paste(rows$id, rows$uid)
rows$combo <- paste0(rows$variant, "_F", rows$F)
n_combo <- tapply(rows$combo, rows$key, function(x) length(unique(x)))
full <- names(n_combo)[n_combo == 9L]
cat("later passages (<= 60 days) with all 9 forecasts:", length(full), "of", length(n_combo), "\n")
rows <- rows[rows$key %in% full, ]

## ---------- base check against decomp ----------
dec <- readRDS(file.path(root, "decomp", "evaluation.rds"))$rows
chk <- do.call(rbind, lapply(c(1, 16), function(f) {
  a <- rows[rows$variant == "base" & rows$F == f, c("key", "ll_fc_s_e2")]
  b <- dec[dec$L == 1 & dec$F == f, c("key", "ll_fc_s_e2")]
  m <- merge(a, b, by = "key", suffixes = c("_disp", "_decomp"))
  data.frame(F = f, n = nrow(m), max_abs_diff = max(abs(m$ll_fc_s_e2_disp - m$ll_fc_s_e2_decomp)))
}))
cat("\n== base check (must be 0) ==\n"); print(chk)
if (any(chk$max_abs_diff > 1e-8) || any(chk$n < length(full))) stop("base does not reproduce the decomposition study (PLAN stop rule)")

base1 <- rows[rows$variant == "base" & rows$F == 1, c("key", "ll_fc_s_e2")]
rows$gain <- rows$ll_fc_s_e2 - base1$ll_fc_s_e2[match(rows$key, base1$key)]

## ---------- D1: where the F x16 gain comes from ----------
kary$key <- paste(kary$id, kary$uid)
k1 <- kary[kary$key %in% full, ]
d1 <- do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
  x <- k1[k1$group == gr & k1$variant == "base" & k1$in_support, ]
  cells <- stats::aggregate(y ~ id + uid, x[x$F == 1, ], sum)
  tot <- sum(cells$y)
  do.call(rbind, lapply(c("A", "B", "C", "D"), function(cl) {
    a <- x[x$class == cl & x$F == 1, ]; b <- x[x$class == cl & x$F == 16, ]
    data.frame(group = gr, class = cl, cell_share = sum(a$y) / tot,
               ll_per_cell_F1 = sum(a$contrib) / sum(a$y), ll_per_cell_F16 = sum(b$contrib) / sum(b$y),
               contribution_to_gain = (sum(b$contrib) - sum(a$contrib)) / tot)
  }))
}))
cat("\n== D1: base F x1 -> x16, by class of the later passage's cells (in the landscape) ==\n"); print(d1, digits = 3, row.names = FALSE)
# per passage, the gain is the per-passage mean; D1's contribution_to_gain is pooled over cells (weights differ)

## ---------- D2: novelty vs Good-Turing unseen share ----------
kb <- k1[k1$variant == "base" & k1$F == 1, ]
nov <- stats::aggregate(cbind(cells = y, novel = y * (class != "A")) ~ id + uid + days + line + group, kb, sum)
gt <- unique(rows[rows$variant == "base" & rows$F == 1, c("id", "uid", "gt_unseen", "n_parent", "coverage_s")])
nov <- merge(nov, gt, by = c("id", "uid"))
nov$novel_share <- nov$novel / nov$cells
nov$bin <- days_bin(nov$days)
d2 <- stats::aggregate(cbind(novel_share, gt_unseen, outside_landscape = 1 - coverage_s) ~ group + bin, nov, mean)
d2$n <- stats::aggregate(novel_share ~ group + bin, nov, length)$novel_share
cat("\n== D2: share of later-passage cells in karyotypes the parent sample lacks, vs the parent's Good-Turing unseen share ==\n")
print(d2, digits = 3, row.names = FALSE)
d2_class <- stats::aggregate(y ~ group + class, kb, sum)
d2_class$share <- d2_class$y / ave(d2_class$y, d2_class$group, FUN = sum)
cat("\nall cells (inside and outside the landscape) by class:\n"); print(d2_class, digits = 3, row.names = FALSE)

## ---------- D3: distance spectrum ----------
kb$dclass <- cut(kb$dist, c(-Inf, 0, 1, 2, Inf), labels = c("0", "1", "2", "3+"))
obs <- stats::aggregate(y ~ group + dclass, kb, sum)
obs$share <- obs$y / ave(obs$y, obs$group, FUN = sum)
obs_in <- stats::aggregate(y ~ group + dclass, kb[kb$in_support, ], sum)
obs_in$share_in_landscape <- obs_in$y / ave(obs_in$y, obs_in$group, FUN = sum)
dm <- dmass[dmass$variant == "base" & paste(dmass$id, dmass$uid) %in% full, ]
dm$mass[is.na(dm$mass)] <- 0   # no landscape karyotype at that distance
fm <- stats::aggregate(mass ~ group + dist + F, dm, mean)
fm <- stats::reshape(fm, idvar = c("group", "dist"), timevar = "F", direction = "wide")
names(fm)[names(fm) == "dist"] <- "dclass"
d3 <- merge(merge(obs[c("group", "dclass", "share")], obs_in[c("group", "dclass", "share_in_landscape")], all = TRUE),
            fm, by = c("group", "dclass"), all = TRUE)
cat("\n== D3: distance to the nearest parent-sample karyotype: observed cell share (all; in landscape) vs forecast mass (base F x1, x16) ==\n")
print(d3, digits = 3, row.names = FALSE)

## ---------- D4: drift dispersion on class A ----------
ka <- k1[k1$variant == "base" & k1$F == 1 & k1$class == "A" & k1$in_support, ]
ncell <- stats::aggregate(y ~ id + uid, k1[k1$variant == "base" & k1$F == 1 & k1$in_support, ], sum)
ka$n <- ncell$y[match(paste(ka$id, ka$uid), paste(ncell$id, ncell$uid))]
ka$e <- ka$n * ka$q
phi <- do.call(rbind, lapply(split(ka, paste(ka$id, ka$uid)), function(x) {
  x <- x[x$e >= 5, ]
  if (nrow(x) < 2) return(NULL)
  data.frame(id = x$id[1], uid = x$uid[1], days = x$days[1], group = x$group[1], n = x$n[1],
             phi = sum((x$y - x$e)^2 / x$e) / (nrow(x) - 1), terms = nrow(x))
}))
phi$passages <- phi$days / ifelse(salehi_a$line[match(phi$id, salehi_a$id)] %in% c("SA906", "SA039"), 5, 15)
phi$ne <- ifelse(phi$phi > 1, phi$n * phi$passages / (phi$phi - 1), Inf)
phi$bin <- days_bin(phi$days)
d4 <- stats::aggregate(cbind(phi, ne) ~ group + bin, phi, stats::median)
d4$n <- stats::aggregate(phi ~ group + bin, phi, length)$phi
cat("\n== D4: Pearson dispersion of class-A karyotypes (expected >= 5) under base F x1; implied N_e per passage (medians) ==\n")
print(d4, digits = 3, row.names = FALSE)

## ---------- interventions ----------
set.seed(20261120)
boot_ci <- function(g, ids, reps = 2000L) {
  u <- unique(ids); by_id <- split(g, ids)
  b <- replicate(reps, { s <- sample(u, length(u), replace = TRUE); mean(unlist(by_id[s], use.names = FALSE)) })
  stats::quantile(b, c(0.025, 0.975), names = FALSE)
}
rows$line <- salehi_a$line[match(rows$id, salehi_a$id)]
rows$group <- group_of(rows$line)
iv <- do.call(rbind, lapply(split(rows, rows$combo), function(x) do.call(rbind, lapply(c("dev", "confirm"), function(gr) {
  y <- x[x$group == gr, ]
  ci <- boot_ci(y$gain, y$id)
  data.frame(combo = x$combo[1], variant = x$variant[1], F = x$F[1], group = gr, passages = nrow(y),
             gain = mean(y$gain), lo = ci[1], hi = ci[2])
}))))
iv <- iv[order(iv$group, iv$variant, iv$F), ]
cat("\n== gain over base F x1 (per later passage), 95% trajectory bootstrap ==\n"); print(iv, digits = 3, row.names = FALSE)
by_line <- stats::aggregate(gain ~ combo + line, rows, mean)
by_line <- stats::reshape(by_line, idvar = "combo", timevar = "line", direction = "wide")
cat("\n== gain by line ==\n"); print(by_line, digits = 2, row.names = FALSE)

verdict <- do.call(rbind, lapply(c("smooth", "drift1e3", "drift1e4"), function(v) {
  c1 <- iv[iv$variant == v & iv$F == 1 & iv$group == "confirm", ]
  d1r <- iv[iv$variant == v & iv$F == 1 & iv$group == "dev", ]
  pass <- c1$gain >= 0.5 * f_only_loo[["confirm"]] && c1$lo > 0 && d1r$gain >= 0
  data.frame(intervention = v, confirm_gain = c1$gain, confirm_lo = c1$lo, threshold = 0.5 * f_only_loo[["confirm"]],
             dev_gain = d1r$gain, explains = pass)
}))
cat("\n== verdicts (PLAN) ==\n"); print(verdict, digits = 3, row.names = FALSE)

dev_save_rds(list(check = chk, d1 = d1, d2 = d2, d2_class = d2_class, d3 = d3, d4 = d4, phi = phi, interventions = iv,
                  by_line = by_line, verdict = verdict, rows = rows),
             file.path(root, "disp", "evaluation.rds"))
