source("dev/R/common.R")
ev <- readRDS(file.path(dev_results, "realdata", "noise", "v5", "step12_evaluation.rds"))
B <- ev$bounds; b0 <- B[B$variant == "none" & B$version == "a", ]
H <- ev$high; h0 <- H[H$variant == "none", ]
# share of cells in high-frequency karyotypes, per dataset
hs <- aggregate(cbind(x0, x1) ~ dataset_id + set, h0, sum)
b0$high_share <- with(hs, (x0 + x1))[match(b0$dataset_id, hs$dataset_id)]
recs <- unlist(lapply(list.files(file.path(dev_results, "realdata", "noise", "v5", "anchor"), full.names = TRUE), readRDS), recursive = FALSE)
n_of <- vapply(recs, function(r) if (identical(r$status, "ok")) r$n0 + r$n1 else NA_real_, 1); names(n_of) <- vapply(recs, `[[`, "", "dataset_id")
b0$high_share <- ifelse(is.na(b0$high_share), 0, b0$high_share) / n_of[b0$dataset_id]
q <- function(v) sprintf("%.2f / %.2f / %.2f", stats::quantile(v, 0.25, na.rm = TRUE), stats::median(v, na.rm = TRUE), stats::quantile(v, 0.75, na.rm = TRUE))
cat("== per set (none, version a): share of cells in high-frequency karyotypes; r0; r1 (25/50/75%) ==\n")
for (s in sort(unique(b0$set))) { x <- b0[b0$set == s, ]; cat(sprintf("%-18s truth %.1e  high share %s | r0 %s | r1 %s | no ref %.2f\n", s, x$truth[1], q(x$high_share), q(x$r0[is.finite(x$r0)]), q(x$r1[is.finite(x$r1)]), mean(x$n_ref == 0))) }
R <- ev$real; ra <- R[R$version == "a", ]
nr <- vapply(recs, function(r) if (r$kind %in% c("salehi","numbat") && identical(r$status, "ok")) sum(r$high$x0 + r$high$x1) / (r$n0 + r$n1) else NA_real_, 1)
names(nr) <- names(n_of); ra$high_share <- nr[ra$dataset_id]
cat("\n== real data (version a) ==\n")
for (l in unique(ra$line)) { x <- ra[ra$line == l, ]; cat(sprintf("%-9s n %3d  high share %s | r0 %s | r1 %s\n", l, nrow(x), q(x$high_share), q(x$r0[is.finite(x$r0)]), q(x$r1[is.finite(x$r1)]))) }
# Salehi pairs, two-proportion z
sal <- Filter(function(r) r$kind == "salehi", recs); sh <- lapply(sal, `[[`, "high"); names(sh) <- vapply(sal, `[[`, "", "dataset_id")
P <- ev$pairs[!is.na(ev$pairs$karyotype), ]
ids <- do.call(rbind, strsplit(P$pair, " ", fixed = TRUE))
P$z2a <- mapply(function(i, k) sh[[i]]$z2[match(k, sh[[i]]$karyotype)], ids[, 1], P$karyotype)
P$z2b <- mapply(function(i, k) sh[[i]]$z2[match(k, sh[[i]]$karyotype)], ids[, 2], P$karyotype)
both <- P[abs(P$z2a) > 3 & abs(P$z2b) > 3, ]
cat(sprintf("\nSalehi pairs, |two-proportion z| > 3 in both: %d karyotypes (%d pairs); same direction %.3f; rises among them %.2f\n",
            nrow(both), length(unique(both$pair)), mean(sign(both$z2a) == sign(both$z2b)), mean(both$z2a > 0)))
print(tapply(sign(both$z2a) == sign(both$z2b), P$line[match(both$pair, P$pair)], function(v) sprintf("%.2f (%d)", mean(v), length(v))))
one <- P[xor(abs(P$z2a) > 3, abs(P$z2b) > 3), ]
cat(sprintf("|z| > 3 in one only: %d of %d (%.2f)\n", nrow(one), nrow(one) + nrow(both), nrow(one) / (nrow(one) + nrow(both))))
cat("\nenvelope falls flagged among high-frequency (all sims, none):", sum(h0$abnormal & h0$z < 0), " testable falls:", sum(h0$x0 > 0 & h0$testable & h0$x1 < h0$x0), "\n")
# selection sets, two-proportion direction
sel <- h0[h0$set != "e3" & abs(h0$z2) > 3, ]
sel$agree <- ifelse(sel$z2 > 0, sel$f_true > sel$f_mean, sel$f_true < sel$f_mean)
cat(sprintf("selection sets, |two-proportion z| > 3: %d, agree %.3f (rises %.3f, falls %.3f)\n", nrow(sel), mean(sel$agree), mean(sel$agree[sel$z2 > 0]), mean(sel$agree[sel$z2 < 0])))
print(tapply(sel$agree, sel$bottleneck, mean), digits = 3)
e3 <- h0[h0$set == "e3", ]; cat("e3 share with |two-proportion z| > 3 by bottleneck:\n"); print(tapply(abs(e3$z2) > 3, e3$bottleneck, mean), digits = 3)
