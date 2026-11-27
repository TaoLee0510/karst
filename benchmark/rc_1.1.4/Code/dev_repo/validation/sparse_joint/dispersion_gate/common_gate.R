# Shared code of the dispersion gate (PLAN.md): the two dispersion statistics
# of a two-time-point dataset.

# Subsample a vector of counts named by karyotype to n cells, without
# replacement; returns the counts of the subsample (named, zeros dropped).
gate_subsample <- function(v, n) {
  v <- v[v > 0]
  cells <- rep.int(seq_along(v), v)
  s <- tabulate(sample(cells, n), nbins = length(v))
  names(s) <- names(v)
  s[s > 0]
}

# once: share of the 2n subsampled cells whose karyotype occurs exactly once
# in its own subsample (t0 and t1 pooled). new: share of the n t1 cells whose
# karyotype does not occur in the t0 subsample. Mean over `draws` draws
# (seeds 1..draws); NA when a time point has fewer than n cells.
gate_stats <- function(c0, c1, n, draws = 10L) {
  if (sum(c0) < n || sum(c1) < n) return(c(once = NA_real_, new = NA_real_, simpson = NA_real_, top = NA_real_))
  r <- vapply(seq_len(draws), function(d) {
    set.seed(d)
    a <- gate_subsample(c0, n)
    b <- gate_subsample(c1, n)
    pb <- b / n
    c(once = (sum(a == 1L) + sum(b == 1L)) / (2 * n),
      new = sum(b[!names(b) %in% names(a)]) / n,
      simpson = 1 - sum(pb^2), top = max(pb))
  }, numeric(4))
  rowMeans(r)
}

# Counts of the two input time points as named vectors.
gate_columns <- function(counts) {
  counts <- as.matrix(counts)
  list(c0 = stats::setNames(as.numeric(counts[, 1L]), rownames(counts)),
       c1 = stats::setNames(as.numeric(counts[, ncol(counts)]), rownames(counts)))
}
