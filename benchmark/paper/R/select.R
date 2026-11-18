# Selection by cross-validation score, the paper's rule (each input's highest
# score), with ties resolved.
#
# The paper keeps every option that ties for the highest score (S05:
# `Rxv == max(Rxv)`). alfak fits each minobs on its own and its scores never
# tie. KARST derives every minobs from one fit; the cross-validation
# predictions are then often the same for two or three minobs, and their
# scores differ only by floating-point rounding, so that a strict comparison
# keeps one, two or three of them depending on the rounding, and the input is
# counted that many times. Here scores within `tol` of a group's highest count
# as tied, and of the tied options the first in `rank` order is kept (the
# lowest minobs; minobs of a tie can list different karyotypes, so the choice
# is fixed). Every row of the chosen option is kept (e.g. every passage and
# metric of a forecast). For alfak this gives the paper's selection.
#
#   d       rows, with the option's score in column `score`
#   by      columns of the groups (one choice per group)
#   choice  column naming the option (e.g. "minobs", "fi")
#   rank    column ordering the options; numeric where it parses as a number
# An option without a score (NA) is never chosen; a group with no scored
# option gives no row. A missing value in a `by` column is an error: such a
# row belongs to no group (split() would drop it without a word).
select_by_cv <- function(d, by, choice, score, tol = 1e-8, rank = choice) {
  if (is.null(d) || !nrow(d)) return(d)
  na_by <- !stats::complete.cases(d[by])
  if (any(na_by)) stop("select_by_cv: ", sum(na_by), " rows with a missing value in ", paste(by, collapse = ", "))
  groups <- split(seq_len(nrow(d)), d[by], drop = TRUE)
  keep <- unlist(lapply(groups, function(i) {
    i <- i[!is.na(d[[score]][i])]
    if (!length(i)) return(integer())
    s <- d[[score]][i]
    tied <- i[s >= max(s) - tol]
    r <- d[[rank]][tied]
    rn <- suppressWarnings(as.numeric(r))
    first <- d[[choice]][tied][order(rn, r)][1]
    i[d[[choice]][i] == first]
  }), use.names = FALSE)
  d[sort(keep), , drop = FALSE]
}
