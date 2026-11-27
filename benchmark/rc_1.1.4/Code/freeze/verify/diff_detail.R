# Detail of the differences compare.R found: for every file of COMPARE_DIR's
# pair tables that is neither identical nor equal, the leaves that differ.
#   Rscript diff_detail.R PAIRS.tsv COMPARE_DIR OUT.tsv
# One row per differing leaf: file, leaf (path inside the object or table
# column), kind, detail. Kinds:
#   path       character values that are file paths in both (as-run root vs freeze root)
#   float      numbers, largest relative difference below 1e-6
#   numeric    numbers, larger differences
#   text       other character values
#   structure  components or columns present on one side only, lengths or types that differ
#   rows       tables with a different number of rows
# Files outside the rerun folder (mirror links, linked input folders) are
# counted as "input_link" and not compared.
args <- commandArgs(TRUE)
pairs <- utils::read.delim(args[[1]], stringsAsFactors = FALSE); cdir <- args[[2]]; out <- args[[3]]
timing <- "^(seconds|.*_seconds|elapsed|started|finished|host|hostname|timestamp|.*peak.*|runtime|run_time|pid|tempdir|work_dir)$"
is_path <- function(v) is.character(v) && length(v) > 0L && all(is.na(v) | grepl("^/|^[A-Za-z0-9_.-]+/", v))
rows <- list()
add <- function(file, leaf, kind, detail) rows[[length(rows) + 1L]] <<- data.frame(file = file, leaf = leaf, kind = kind, detail = substr(detail, 1, 300), stringsAsFactors = FALSE)
num_detail <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  na <- xor(is.na(a), is.na(b)); ok <- !is.na(a) & !is.na(b) & is.finite(a) & is.finite(b)
  d <- abs(a[ok] - b[ok]); r <- d / pmax(abs(a[ok]), abs(b[ok]), 1e-300)
  both <- !is.na(a) & !is.na(b)
  list(n = sum(na) + sum(a[both] != b[both]), na = sum(na),
       abs = if (length(d)) max(d) else 0, rel = if (length(r)) max(r) else 0)
}
walk <- function(file, x, y, leaf, depth = 0L) {
  if (identical(x, y)) return(invisible())
  if (depth > 15L) return(add(file, leaf, "structure", "too deep"))
  if (is.environment(x) || is.function(x)) return(invisible())
  if (is.list(x) && is.list(y)) {
    nx <- names(x); ny <- names(y)
    if (!is.null(nx) && !is.null(ny)) {
      tx <- nx[grepl(timing, nx, ignore.case = TRUE)]
      only_x <- setdiff(nx, ny); only_y <- setdiff(ny, nx)
      if (length(only_x) || length(only_y))
        add(file, leaf, "structure", paste0("rerun only: ", paste(head(only_x, 12), collapse = ","), " | frozen only: ", paste(head(only_y, 12), collapse = ",")))
      if (is.data.frame(x) && is.data.frame(y) && nrow(x) != nrow(y)) return(add(file, leaf, "rows", paste(nrow(x), "rows vs", nrow(y))))
      for (n in setdiff(intersect(nx, ny), tx)) if (nzchar(n)) walk(file, x[[n]], y[[n]], paste0(leaf, "$", n), depth + 1L)
      return(invisible())
    }
    if (length(x) != length(y)) return(add(file, leaf, "structure", paste("list length", length(x), "vs", length(y))))
    for (i in seq_along(x)) walk(file, x[[i]], y[[i]], paste0(leaf, "[[", i, "]]"), depth + 1L)
    return(invisible())
  }
  if (inherits(x, "proc_time") || inherits(x, "POSIXt") || inherits(x, "difftime")) return(invisible())
  if (is.atomic(x) && is.atomic(y)) {
    if (length(x) != length(y)) return(add(file, leaf, "structure", paste("length", length(x), "vs", length(y))))
    if (is.numeric(x) && is.numeric(y)) {
      if (identical(as.numeric(x), as.numeric(y))) return(if (!identical(attributes(x), attributes(y))) add(file, leaf, "structure", "attributes") else invisible())
      d <- num_detail(x, y)
      return(add(file, leaf, if (d$na == 0 && d$rel < 1e-6) "float" else "numeric",
                 sprintf("%d of %d differ, %d NA pattern, max abs %.3g, max rel %.3g", d$n, length(x), d$na, d$abs, d$rel)))
    }
    if (is.character(x) && is.character(y)) {
      m <- which(!(x == y) | xor(is.na(x), is.na(y))); if (!length(m)) return(if (!identical(attributes(x), attributes(y))) add(file, leaf, "structure", "attributes") else invisible())
      k <- if (is_path(x[m]) && is_path(y[m]) && identical(basename(x[m]), basename(y[m]))) "path" else "text"
      return(add(file, leaf, k, sprintf("%d of %d differ: '%s' vs '%s'", length(m), length(x), x[m[1]], y[m[1]])))
    }
    if (identical(unclass(x), unclass(y))) return(add(file, leaf, "structure", "attributes"))
    return(add(file, leaf, "structure", paste("types", typeof(x), "vs", typeof(y), "or values differ")))
  }
  add(file, leaf, "structure", paste("classes", paste(class(x), collapse = "/"), "vs", paste(class(y), collapse = "/")))
}
read_tab <- function(f) {
  sep <- if (grepl("[.]csv$", f)) "," else "\t"
  tryCatch(utils::read.table(f, header = TRUE, sep = sep, quote = "", comment.char = "", stringsAsFactors = FALSE,
                             check.names = FALSE, fill = TRUE), error = function(e) NULL)
}
summ <- list()
for (i in seq_len(nrow(pairs))) {
  lab <- pairs$label[i]; rr <- normalizePath(pairs$rerun[i]); fr <- normalizePath(pairs$frozen[i])
  tab <- utils::read.delim(file.path(cdir, paste0(lab, ".tsv")), stringsAsFactors = FALSE, quote = "")
  inside <- startsWith(normalizePath(file.path(rr, tab$file), mustWork = FALSE), paste0(rr, "/"))
  res <- sub(":.*", "", tab$result)
  summ[[i]] <- data.frame(label = lab, compared = sum(inside), identical_or_equal = sum(inside & res %in% c("identical", "equal")),
                          different = sum(inside & !res %in% c("identical", "equal")), input_link = sum(!inside), stringsAsFactors = FALSE)
  for (j in which(inside & !res %in% c("identical", "equal"))) {
    f <- tab$file[j]; a <- file.path(rr, f); b <- file.path(fr, f); file <- paste0(lab, ": ", f)
    if (!file.exists(b)) { add(file, "", "structure", "missing in frozen"); next }
    ext <- tolower(tools::file_ext(f))
    if (ext == "rds") {
      x <- tryCatch(readRDS(a), error = function(e) e); y <- tryCatch(readRDS(b), error = function(e) e)
      if (inherits(x, "error") || inherits(y, "error")) { add(file, "", "structure", "unreadable"); next }
      walk(file, x, y, "")
    } else if (ext %in% c("tsv", "csv")) {
      x <- read_tab(a); y <- read_tab(b)
      if (is.null(x) || is.null(y)) { add(file, "", "text", "not a table"); next }
      walk(file, x, y, "")
    } else if (ext == "pdf") {
      strip_pdf <- function(p) { r <- readBin(p, "raw", file.info(p)$size); s <- rawToChar(r[r != as.raw(0)], multiple = FALSE)
        gsub("/(CreationDate|ModDate) \\(D:[^)]*\\)", "", s, useBytes = TRUE) }
      add(file, "", if (identical(strip_pdf(a), strip_pdf(b))) "pdf_dates" else "pdf_content", "")
    } else {
      la <- readLines(a, warn = FALSE); lb <- readLines(b, warn = FALSE)
      m <- which(la[seq_len(min(length(la), length(lb)))] != lb[seq_len(min(length(la), length(lb)))])
      add(file, "", "text", sprintf("%d vs %d lines, %d differ; first: '%s' vs '%s'", length(la), length(lb), length(m),
                                    if (length(m)) la[m[1]] else "", if (length(m)) lb[m[1]] else ""))
    }
  }
}
d <- do.call(rbind, rows); s <- do.call(rbind, summ)
utils::write.table(d, out, sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(s, sub("[.]tsv$", "_summary.tsv", out), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 250)
print(s, row.names = FALSE)
cat("\nleaves by kind and pair:\n"); print(table(sub(":.*", "", d$file), d$kind))
