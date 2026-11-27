# Verification of batch 1: the own-benchmark scores of the 20-dataset subset
# (rerun from Code/ with STT_DATA_DIR = the subset) against the frozen tables.
#   Rscript own_subset_check.R RERUN_DIR FROZEN_DIR OUT.tsv
args <- commandArgs(TRUE)
chk <- function(file, keys) {
  a <- utils::read.csv(file.path(args[[1]], file), stringsAsFactors = FALSE)
  b <- utils::read.csv(file.path(args[[2]], file), stringsAsFactors = FALSE)
  b <- b[b$dataset_id %in% a$dataset_id, ]
  m <- merge(a, b, by = keys, suffixes = c(".rerun", ".frozen"))
  cols <- setdiff(intersect(names(a), names(b)), c(keys, "seconds"))
  d <- vapply(cols, function(cn) {
    x <- m[[paste0(cn, ".rerun")]]; y <- m[[paste0(cn, ".frozen")]]
    if (is.numeric(x)) max(abs(x - y), 0, na.rm = TRUE) + sum(is.na(x) != is.na(y)) else sum(x != y, na.rm = TRUE) + sum(is.na(x) != is.na(y))
  }, numeric(1))
  data.frame(file = file, rows_rerun = nrow(a), rows_matched = nrow(m), columns = length(cols),
             max_difference = max(d), worst_column = names(d)[which.max(d)], stringsAsFactors = FALSE)
}
res <- rbind(chk("metrics_heldout_grid.csv", c("dataset_id", "multiple", "method", "stratum")),
             chk("direction_heldout_grid.csv", c("dataset_id", "multiple", "method", "stratum")))
utils::write.table(res, args[[3]], sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
