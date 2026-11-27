#!/usr/bin/env Rscript
# Steps 1-2 of dev/realdata/noise/PLAN_v5.md, local, from the records of
# step12_task.R (DEV_RESULTS/realdata/noise/v5/anchor/<tag>.rds).
#   Rscript dev/realdata/noise/v5/step12_eval.R
# Writes DEV_RESULTS/realdata/noise/v5/step12_evaluation.rds and prints the tables.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
source(file.path(dev_repo, "dev", "realdata", "noise", "R", "envelope_v4.R"))
root <- file.path(dev_results, "realdata", "noise", "v5", "anchor")
files <- list.files(root, pattern = "[.]rds$", full.names = TRUE)
recs <- unlist(lapply(files, function(f) lapply(readRDS(f), function(r) c(r, tag = sub("[.]rds$", "", basename(f))))),
               recursive = FALSE)
sim <- Filter(function(r) r$kind %in% c("eset", "c2"), recs)
real <- Filter(function(r) r$kind %in% c("salehi", "numbat"), recs)
cat("records:", length(recs), " simulations:", length(sim), " real:", length(real), "\n")
cat("not ok:", sum(!vapply(recs, function(r) identical(r$status, "ok"), logical(1))), "\n")
set_of <- function(r) {
  if (r$kind == "c2") return(sub("^c2_r1_(.*)_[a-z0-9]+$", "C2 \\1", r$tag))
  if (r$set == "e4") return(if (abs(r$truth - 2e-4) < 1e-12) "e4 x4" else "e4 x16")
  r$set
}
ok_sim <- Filter(function(r) identical(r$status, "ok"), sim)

## ---------- step 1 ----------
high <- do.call(rbind, lapply(ok_sim, function(r) {
  if (!nrow(r$high)) return(NULL)
  data.frame(set = set_of(r), variant = r$variant, half = r$half, depth = r$depth, bottleneck = r$bottleneck,
             dt = r$dt, dataset_id = r$dataset_id, r$high[, c("karyotype", "x0", "x1", "z", "abnormal", "testable", "z2",
                                                              "falling", "f_true")],
             f_mean = r$f_mean, stringsAsFactors = FALSE)
}))
h0 <- high[high$variant == "none", ]
cat("\n== Step 1.1: neutral e3, share of high-frequency karyotypes with a significant change ==\n")
e3 <- h0[h0$set == "e3", ]
s11 <- do.call(rbind, lapply(split(e3, list(e3$bottleneck, e3$depth), drop = TRUE), function(x) {
  per <- tapply(x$abnormal, x$dataset_id, any)
  data.frame(bottleneck = x$bottleneck[1], depth = x$depth[1], high = nrow(x), flagged = sum(x$abnormal),
             share = mean(x$abnormal), share_testable = mean(x$abnormal[x$testable]), datasets = length(per),
             datasets_with_one = mean(per))
}))
all_e3 <- unique(vapply(Filter(function(r) r$variant == "none" && r$set == "e3", ok_sim), `[[`, "", "dataset_id"))
s11 <- s11[order(s11$bottleneck, s11$depth), ]
print(s11, digits = 3, row.names = FALSE)
cat(sprintf("e3 overall: %.3f of %d high-frequency karyotypes (%d datasets, %d with none high)\n", mean(e3$abnormal),
            nrow(e3), length(all_e3), length(all_e3) - length(unique(e3$dataset_id))))
cat("by interval (days):\n"); print(tapply(e3$abnormal, e3$dt, mean), digits = 3)

cat("\n== Step 1.2: selection sets, direction of significant changes vs the truth ==\n")
sel <- h0[h0$set != "e3" & h0$abnormal, ]
sel$agree <- ifelse(sel$z > 0, sel$f_true > sel$f_mean, sel$f_true < sel$f_mean)
s12 <- do.call(rbind, lapply(split(sel, sel$set), function(x)
  data.frame(set = x$set[1], flagged = nrow(x), rises = sum(x$z > 0), agree = mean(x$agree),
             agree_rise = mean(x$agree[x$z > 0]), agree_fall = if (any(x$z < 0)) mean(x$agree[x$z < 0]) else NA)))
print(s12, digits = 3, row.names = FALSE)
cat(sprintf("all selection sets: %.3f of %d\n", mean(sel$agree), nrow(sel)))
cat("by bottleneck (e sets):\n"); es <- sel[!grepl("^C2", sel$set), ]; print(tapply(es$agree, es$bottleneck, mean), digits = 3)
cat("by depth:\n"); print(tapply(sel$agree, sel$depth, mean), digits = 3)
cat("share of high-frequency karyotypes flagged, selection sets:\n")
hs <- h0[h0$set != "e3", ]; print(tapply(hs$abnormal, hs$set, mean), digits = 3)

cat("\n== Step 1.3: Salehi parallel lineages ==\n")
sal <- Filter(function(r) r$kind == "salehi" && identical(r$status, "ok"), real)
sal_high <- lapply(sal, `[[`, "high"); names(sal_high) <- vapply(sal, `[[`, "", "dataset_id")
sa <- readRDS(file.path(dev_results, "realdata", "A", "salehi.rds"))
sa <- sa[sa$id %in% names(sal_high), ]
keys <- split(sa, list(sa$line, sa$treat_first, sa$treat_last), drop = TRUE)
pairs <- do.call(rbind, lapply(keys, function(g) {
  if (nrow(g) < 2) return(NULL)
  cmb <- utils::combn(nrow(g), 2)
  do.call(rbind, lapply(seq_len(ncol(cmb)), function(j) {
    a <- g[cmb[1, j], ]; b <- g[cmb[2, j], ]
    if (length(intersect(c(a$first_uid, a$last_uid), c(b$first_uid, b$last_uid)))) return(NULL)
    ha <- sal_high[[a$id]]; hb <- sal_high[[b$id]]
    common <- intersect(ha$karyotype, hb$karyotype)
    if (!length(common)) return(data.frame(line = a$line, pair = paste(a$id, b$id), karyotype = NA, za = NA, zb = NA,
                                           fa = NA, fb = NA))
    ia <- match(common, ha$karyotype); ib <- match(common, hb$karyotype)
    data.frame(line = a$line, pair = paste(a$id, b$id), karyotype = common, za = ha$z[ia], zb = hb$z[ib],
               fa = ha$abnormal[ia], fb = hb$abnormal[ib])
  }))
}))
n_pairs <- length(unique(pairs$pair))
both <- pairs[!is.na(pairs$fa) & pairs$fa & pairs$fb, ]
one <- pairs[!is.na(pairs$fa) & xor(pairs$fa, pairs$fb), ]
anyf <- pairs[!is.na(pairs$fa) & (pairs$fa | pairs$fb), ]
cat(sprintf("pairs: %d; pairs sharing a high-frequency karyotype: %d\n", n_pairs,
            length(unique(pairs$pair[!is.na(pairs$karyotype)]))))
cat(sprintf("shared high-frequency karyotypes flagged in both: %d (%d pairs), same direction: %.3f\n", nrow(both),
            length(unique(both$pair)), mean(sign(both$za) == sign(both$zb))))
cat(sprintf("flagged in one, high but not flagged in the other: %d of %d flagged in at least one (%.3f)\n",
            nrow(one), nrow(anyf), nrow(one) / max(1, nrow(anyf))))
if (nrow(both)) { cat("by line, same direction:\n"); print(tapply(sign(both$za) == sign(both$zb), both$line, function(v) sprintf("%.2f (%d)", mean(v), length(v)))) }

cat("\n== Step 1.4: references with a true net rate above rho = 1 (simulations, none) ==\n")
refs <- h0[!h0$falling, ]
print(tapply(refs$f_true > 1, refs$set, mean), digits = 3)
cat("quantiles of the references' true net rate:\n"); print(stats::quantile(refs$f_true, c(0, 0.5, 0.9, 0.99, 1)), digits = 3)

## ---------- step 2 ----------
B <- do.call(rbind, lapply(ok_sim, function(r) data.frame(set = set_of(r), variant = r$variant, half = r$half,
                                                         landscape_id = r$landscape_id, depth = r$depth,
                                                         bottleneck = r$bottleneck, dt = r$dt, truth = r$truth,
                                                         dataset_id = r$dataset_id, r$bounds, stringsAsFactors = FALSE)))
B$exceeds <- !is.na(B$L) & B$L > B$truth * (1 + 1e-9)
crit <- function(x) {
  per_set <- tapply(x$exceeds, x$set, mean)
  hi <- x[x$truth >= 4 * 5e-5 - 1e-12, ]
  data.frame(datasets = nrow(x), valid_all = mean(x$exceeds), valid_worst_set = max(per_set),
             worst_set = names(per_set)[which.max(per_set)], useful = mean(hi$triggered), useful_n = nrow(hi),
             triggered_all = mean(x$triggered), pass_valid = mean(x$exceeds) <= 0.05 && max(per_set) <= 0.10,
             pass_useful = mean(hi$triggered) >= 0.20, stringsAsFactors = FALSE)
}
b0 <- B[B$variant == "none", ]
dev_c <- do.call(rbind, lapply(c("a", "b"), function(v) cbind(version = v, crit(b0[b0$half == "dev" & b0$version == v, ]))))
cat("\n== Step 2, development half ==\n"); print(dev_c, digits = 3, row.names = FALSE)
cand <- dev_c[dev_c$pass_valid & dev_c$pass_useful, ]
chosen <- if (!nrow(cand)) NA_character_ else {
  best <- cand[cand$useful == max(cand$useful), ]
  if ("b" %in% best$version) "b" else best$version[1]
}
conf_c <- NULL; verdict <- "step 2 failed on the development half"
if (!is.na(chosen)) {
  conf_c <- cbind(version = chosen, crit(b0[b0$half == "confirm" & b0$version == chosen, ]))
  cat("\n== Step 2, confirmation half, version", chosen, "==\n"); print(conf_c, digits = 3, row.names = FALSE)
  verdict <- if (conf_c$pass_valid && conf_c$pass_useful) "step 2 passed" else "step 2 failed on the confirmation half"
}
cat("\nVERDICT:", verdict, " (chosen version:", chosen, ")\n")

cat("\n== Report only: per set (none, both halves) ==\n")
per <- do.call(rbind, lapply(split(b0, list(b0$set, b0$version), drop = TRUE), function(x)
  data.frame(set = x$set[1], version = x$version[1], n = nrow(x), truth = x$truth[1], no_ref = mean(x$n_ref == 0),
             no_bound = mean(is.na(x$L)), point_pos = mean(x$point > 0, na.rm = TRUE), exceeds = mean(x$exceeds),
             triggered = mean(x$triggered), median_point_x = stats::median(x$point / 5e-5, na.rm = TRUE))))
print(per[order(per$version, per$set), ], digits = 3, row.names = FALSE)
cat("\ntriggered by interval (none, version a):\n"); print(tapply(b0$triggered[b0$version == "a"], b0$dt[b0$version == "a"], mean), digits = 3)
cat("triggered by depth (none, version a):\n"); print(tapply(b0$triggered[b0$version == "a"], b0$depth[b0$version == "a"], mean), digits = 3)
cat("\n== Report only: injected variants ==\n")
inj <- do.call(rbind, lapply(split(B[B$variant != "none", ], list(B$variant[B$variant != "none"], B$version[B$variant != "none"]), drop = TRUE),
                             function(x) cbind(variant = x$variant[1], version = x$version[1], crit(x))))
print(inj, digits = 3, row.names = FALSE)

cat("\n== Report only: real data ==\n")
R <- do.call(rbind, lapply(Filter(function(r) identical(r$status, "ok"), real), function(r)
  data.frame(kind = r$kind, line = if (r$kind == "salehi") r$line else r$cancer,
             group = if (r$kind == "salehi") r$group else NA, dataset_id = r$dataset_id, dt = r$dt, r$bounds,
             stringsAsFactors = FALSE)))
rs <- do.call(rbind, lapply(split(R, list(R$kind, R$line, R$version), drop = TRUE), function(x)
  data.frame(kind = x$kind[1], line = x$line[1], version = x$version[1], n = nrow(x), no_ref = sum(x$n_ref == 0),
             no_bound = sum(is.na(x$L)), point_pos = sum(x$point > 0, na.rm = TRUE), triggered = sum(x$triggered),
             L_median_x = stats::median(x$L / 5e-5, na.rm = TRUE), L_max_x = suppressWarnings(max(x$L / 5e-5, na.rm = TRUE)))))
print(rs[order(rs$version, rs$kind, rs$line), ], digits = 3, row.names = FALSE)

dev_save_rds(list(high = high, pairs = pairs, bounds = B, real = R, step1 = list(e3 = s11, selection = s12),
                  dev = dev_c, confirm = conf_c, chosen = chosen, verdict = verdict, per_set = per, injected = inj,
                  real_summary = rs),
             file.path(dev_results, "realdata", "noise", "v5", "step12_evaluation.rds"))
