# Verification of batch 1: the reproduction check rerun from Code/ (60 refits
# with karst 50e5889 compared with the frozen fits).
#   Rscript rc_check_summary.R CHECK_DIR OUT.tsv
args <- commandArgs(TRUE)
files <- list.files(args[[1]], pattern = "^check_[0-9]+\\.rds$", full.names = TRUE)
res <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  data.frame(check = basename(f), type = r$task$type, path = r$task$path, karst = r$karst, identical = r$identical,
             max_difference = max(r$diff), stringsAsFactors = FALSE)
}))
utils::write.table(res, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(res), "checks;", sum(res$identical), "identical; largest difference", max(res$max_difference), "\n")
print(table(res$type, res$identical))
