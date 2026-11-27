# Check that alfak and karst score the same karyotypes.
#
#   Rscript scripts/check_ksets.R [SOURCE ...]
#
# For every input fitted by both methods (each design) and every
# minobs where both fits succeeded, compares the landscapes' karyotypes and
# their frequent (fq) and neighbour (nn) flags. Writes
# BENCH_RESULTS/tables/ksets.Rds and prints the inputs where they differ.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
sources <- commandArgs(TRUE)
if (!length(sources)) sources <- c("abm", "misseg", "salehi")

rows <- list()
for (src in sources) {
  for (d2 in Sys.glob(file.path(paths$results, "fits", src, "*", "*", "karst", "minobs_*"))) {
    # alfak's folder of the same input and minobs (a path substitution would
    # also hit the project folder, .../Project/karst/)
    dR <- file.path(dirname(dirname(d2)), "alfak", basename(d2))
    ok <- function(d) {
      p <- file.path(d, "run.Rds")
      file.exists(p) && identical(readRDS(p)$status, "ok")
    }
    if (!ok(d2) || !ok(dR)) next
    a <- readRDS(file.path(dR, "landscape.Rds"))
    b <- readRDS(file.path(d2, "landscape.Rds"))
    p <- rev(strsplit(d2, "/", fixed = TRUE)[[1]])
    rows[[length(rows) + 1L]] <- data.frame(
      source = src, dataset = p[4], design = p[3], minobs = sub("minobs_", "", p[1]),
      n_alfak = nrow(a), n_karst = nrow(b), same_k = setequal(a$k, b$k),
      same_fq = setequal(a$k[a$fq], b$k[b$fq]), same_nn = setequal(a$k[a$nn], b$k[b$nn]),
      stringsAsFactors = FALSE)
  }
}
rows <- do.call(rbind, rows)
if (is.null(rows)) stop("No input has both an alfak and a karst fit to compare.")
dir.create(file.path(paths$results, "tables"), recursive = TRUE, showWarnings = FALSE)
saveRDS(rows, file.path(paths$results, "tables", "ksets.Rds"))
cat("Inputs compared:", nrow(rows), "\n")
print(table(source = rows$source, design = rows$design, identical = rows$same_k & rows$same_fq & rows$same_nn))
bad <- rows[!(rows$same_k & rows$same_fq & rows$same_nn), ]
if (nrow(bad)) print(bad)
