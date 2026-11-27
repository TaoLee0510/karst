source("dev/R/common.R")
root <- file.path(dev_results, "realdata", "noise", "v5", "anchor")
rows <- do.call(rbind, lapply(list.files(root, pattern = "_none[.]rds$", full.names = TRUE), function(f) {
  do.call(rbind, lapply(readRDS(f), function(r) {
    if (!identical(r$status, "ok")) return(NULL)
    set <- if (r$kind == "c2") sub("^c2_r1_(.*)_none[.]rds$", "C2 \\1", basename(f)) else if (r$set == "e4") (if (abs(r$truth - 2e-4) < 1e-12) "e4 x4" else "e4 x16") else r$set
    b <- r$bounds[r$bounds$version == "a", ]
    data.frame(set = set, truth = r$truth, high_share = sum(r$high$x0 + r$high$x1) / (r$n0 + r$n1), n_high = nrow(r$high),
               r0 = b$r0, r1 = b$r1, no_ref = b$n_ref == 0)
  }))
}))
q <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return("—"); sprintf("%.2f / %.2f / %.2f", stats::quantile(v, 0.25), stats::median(v), stats::quantile(v, 0.75)) }
for (s in unique(rows$set[order(rows$truth, rows$set)])) { x <- rows[rows$set == s, ]
  cat(sprintf("%-18s truth %.1e n %3d | high share %s | n_high median %g | r0 %s | r1 %s | no ref %.2f\n", s, x$truth[1], nrow(x), q(x$high_share), stats::median(x$n_high), q(x$r0), q(x$r1), mean(x$no_ref))) }
