# Cluster bootstrap of a rule's gain over current (PLAN.md, K1).
# Gain of a rule over current: the difference is taken per input first, then
# averaged within each cluster (landscape_id; the underlying simulation), and
# clusters are weighted equally. The bootstrap draws whole clusters with
# replacement (2,000 times): every row of a drawn cluster (all its beta arms
# and designs) enters as many times as the cluster is drawn. `draws = TRUE`
# also returns the row multiplicities of the first draw (for the check).
boot_gain <- function(x, col, cluster, B = 2000L, seed = 20261125L, draws = FALSE) {
  x$diff <- x[[col]] - x$current
  x <- x[is.finite(x$diff), , drop = FALSE]
  ids <- unique(x[[cluster]])
  rows_of <- split(seq_len(nrow(x)), factor(x[[cluster]], levels = ids))
  cluster_mean <- vapply(rows_of, function(i) mean(x$diff[i]), numeric(1))
  set.seed(seed)
  drawn <- replicate(B, sample.int(length(ids), replace = TRUE), simplify = FALSE)
  b <- vapply(drawn, function(d) mean(cluster_mean[d]), numeric(1))
  out <- c(gain = mean(cluster_mean), low = unname(stats::quantile(b, 0.025)), high = unname(stats::quantile(b, 0.975)), units = length(ids))
  if (draws) attr(out, "first_draw") <- data.frame(x[unlist(rows_of[drawn[[1]]]), c(cluster, intersect(c("arm_set", "design"), names(x))), drop = FALSE])
  out
}
