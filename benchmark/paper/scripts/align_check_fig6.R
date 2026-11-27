# Where the dominant karyotypes of version B's Fig. 6 come from
# (docs/align_karst_6ae2a24.md, section 7):
#
#   ALIGN=karst Rscript scripts/align_check_fig6.R
#   Rscript scripts/align_check_fig6.R       (version A, for comparison)
#
# For the Fig. 6 results of alfak's and KARST's anchored arms (both designs)
# on the extended lists: the dominant karyotype at each missegregation rate,
# its distance class to the observed karyotypes (0, 1, 2), for alfak whether
# its fitness is alfak's landscape mean ("landscape") or its Kriging
# ("kriging"), its fitness rank on the list and its share. Prints the counts
# at the lowest and highest rate and saves fig6_dominant_source.Rds into the
# results folder. In version A (no `dist`) the class is the landscape's own
# sets: "fq" (frequent), "nn" (their neighbours), "far" (the rest).
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
rows <- list()
for (key in c("alfak_2", "alfak_all", "karst_2", "karst_all")) {
  arm <- fig456_arms[fig456_arms$key == key, ]
  res <- readRDS(file.path(analysis_root(paths), "fig456", "fig6", paste0(key, ".Rds")))
  for (fi in names(res)) {
    m <- res[[fi]]
    l <- readRDS(file.path(fig456_fit_dir(arm, fi, attr(m, "min_obs"), paths), "landscape.Rds"))
    top <- rownames(m)[apply(m, 2, which.max)]
    i <- match(top, l$k)
    rk <- rank(-l$mean, ties.method = "min")
    cls <- if (!is.null(l$dist)) l$dist[i] else ifelse(l$fq[i], "fq", ifelse(l$nn[i], "nn", "far"))
    rows[[length(rows) + 1]] <- data.frame(arm = key, fi = fi, rate = seq_along(top), dist = cls,
      source = if (!is.null(l$source)) l$source[i] else if (grepl("alfak", key)) "landscape" else "karst", fitness_rank = rk[i], n = nrow(l),
      share = apply(m, 2, max), stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, rows)
options(width = 200)
for (r in c(1, 8)) {
  cat("rate", r, "\n")
  x <- d[d$rate == r, ]
  print(table(x$arm, paste(x$source, x$dist)))
  print(aggregate(cbind(fitness_rank, share) ~ arm, x, median))
}
saveRDS(d, file.path(paths$results, "fig6_dominant_source.Rds"))
