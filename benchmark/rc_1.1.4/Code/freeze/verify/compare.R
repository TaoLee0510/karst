# Compare a rerun (from Code/ and data/) with the frozen results.
#   Rscript compare.R RERUN FROZEN OUT.tsv
# RERUN and FROZEN are two files or two folders (every file under RERUN is
# compared with the same relative path under FROZEN). R objects (.rds/.Rds)
# must be identical after the run-time fields (seconds, times, host, peak
# memory) are dropped; text tables (.tsv/.csv/.md/.txt/.json) must be equal
# line by line except lines naming when they were generated; anything else
# must have the same sha256. Writes one row per file and prints a summary.
#   Rscript compare.R --pairs PAIRS.tsv OUT_DIR
# compares every row of PAIRS.tsv (columns label, rerun, frozen) and writes
# OUT_DIR/<label>.tsv plus OUT_DIR/summary.tsv.
args <- commandArgs(TRUE)
# only fields that record the run itself; a real difference elsewhere shows up
timing <- "^(seconds|.*_seconds|elapsed|started|finished|host|hostname|timestamp|.*peak.*|runtime|run_time|pid|tempdir|work_dir)$"
strip <- function(x, depth = 0L) {
  if (depth > 12L) return(x)
  if (is.environment(x) || is.function(x)) return(NULL)
  if (is.list(x)) {
    if (!is.null(names(x))) x <- x[!grepl(timing, names(x), ignore.case = TRUE)]
    if (is.data.frame(x)) return(x)
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[i] <- list(strip(x[[i]], depth + 1L))
  }
  if (inherits(x, "proc_time") || inherits(x, "POSIXt") || inherits(x, "difftime")) return(NULL)
  x
}
text_lines <- function(f) { l <- readLines(f, warn = FALSE); l[!grepl("^Generated |generated at|^Finished |^Host ", l)] }
one <- function(a, b) {
  if (!file.exists(b)) return("missing_in_frozen")
  ext <- tolower(tools::file_ext(a))
  if (ext == "rds") {
    x <- tryCatch(strip(readRDS(a)), error = function(e) e); y <- tryCatch(strip(readRDS(b)), error = function(e) e)
    if (inherits(x, "error") || inherits(y, "error")) return("unreadable")
    if (identical(x, y)) return("identical")
    d <- all.equal(x, y, tolerance = 0, check.attributes = TRUE)
    return(if (isTRUE(d)) "equal" else paste("DIFFERENT:", substr(paste(d, collapse = "; "), 1, 300)))
  }
  if (ext %in% c("tsv", "csv", "md", "txt", "json")) return(if (identical(text_lines(a), text_lines(b))) "identical" else "DIFFERENT: text")
  if (identical(unname(tools::md5sum(a)), unname(tools::md5sum(b)))) "identical" else "DIFFERENT: bytes"
}
compare_one <- function(rerun, frozen, out) {
  rerun <- normalizePath(rerun); frozen <- normalizePath(frozen)
  files <- if (dir.exists(rerun)) list.files(rerun, recursive = TRUE, all.files = FALSE) else ""
  # only files the rerun wrote: not the mirror's links to results_rc (inputs), logs or job registries
  # (also files inside a linked folder, e.g. the dev inputs A/ and CE/: they resolve outside RERUN)
  if (dir.exists(rerun)) files <- files[!grepl("(^|/)logs/", files) & basename(files) != "jobs.tsv" &
                                          !nzchar(Sys.readlink(file.path(rerun, files))) &
                                          startsWith(normalizePath(file.path(rerun, files), mustWork = FALSE), paste0(rerun, "/"))]
  res <- data.frame(file = files, result = vapply(files, function(f) {
    a <- if (nzchar(f)) file.path(rerun, f) else rerun; b <- if (nzchar(f)) file.path(frozen, f) else frozen; one(a, b)
  }, ""), stringsAsFactors = FALSE)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(rerun, "vs", frozen, ":", nrow(res), "files\n"); print(table(sub(":.*", "", res$result)))
  bad <- res[!res$result %in% c("identical", "equal"), ]
  if (nrow(bad)) print(utils::head(bad, 20), row.names = FALSE)
  res
}
if (identical(args[1], "--pairs")) {
  pairs <- utils::read.delim(args[2], stringsAsFactors = FALSE); dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
  summ <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    r <- compare_one(pairs$rerun[i], pairs$frozen[i], file.path(args[3], paste0(pairs$label[i], ".tsv")))
    data.frame(label = pairs$label[i], files = nrow(r), same = sum(r$result %in% c("identical", "equal")),
               different = sum(!r$result %in% c("identical", "equal")), stringsAsFactors = FALSE)
  }))
  utils::write.table(summ, file.path(args[3], "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(summ, row.names = FALSE)
} else invisible(compare_one(args[[1]], args[[2]], args[[3]]))
