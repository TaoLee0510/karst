# Check that alfak and karst score the same karyotypes.
#
#   Rscript scripts/check_ksets.R [SOURCE ...]
#
# For every input fitted by alfak and by a KARST method (karst,
# karst_pairwise, karst_state_space; each design) and every minobs where
# both fits succeeded, compares the landscapes' karyotypes and their
# frequent (fq) and neighbour (nn) flags. Writes
# BENCH_RESULTS/tables/ksets.Rds (the rows of the sources not checked are
# kept from the table already there), prints the inputs where they differ and
# then stops with an error, so that the analyses after it in
# hpc/collect_analyse.sh do not run (they also check the table themselves:
# require_same_ksets(), R/xv_source.R).

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
sources <- commandArgs(TRUE)
if (!length(sources)) sources <- c("abm", "misseg", "salehi")

# folders in parallel (forked workers, as many as the job's CPUs: version
# B's landscapes hold up to 2 million karyotypes)
cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), na.rm = TRUE)
check_one <- function(src, d2) {
  # alfak's folder of the same input and minobs (a path substitution would
  # also hit the project folder, .../Project/karst/)
  dR <- file.path(dirname(dirname(d2)), "alfak", basename(d2))
  ok <- function(d) {
    p <- file.path(d, "run.Rds")
    file.exists(p) && identical(readRDS(p)$status, "ok")
  }
  if (!ok(d2) || !ok(dR)) return(NULL)
  a <- readRDS(file.path(dR, "landscape.Rds"))
  b <- readRDS(file.path(d2, "landscape.Rds"))
  p <- rev(strsplit(d2, "/", fixed = TRUE)[[1]])
  data.frame(
    source = src, dataset = p[4], design = p[3], method = p[2], minobs = sub("minobs_", "", p[1]),
    n_alfak = nrow(a), n_karst = nrow(b), same_k = setequal(a$k, b$k),
    same_fq = setequal(a$k[a$fq], b$k[b$fq]), same_nn = setequal(a$k[a$nn], b$k[b$nn]),
    stringsAsFactors = FALSE)
}
rows <- list()
for (src in sources) {
  dirs <- Sys.glob(file.path(paths$results, "fits", src, "*", "*", karst_xv_methods, "minobs_*"))
  part <- parallel::mclapply(dirs, function(d2) check_one(src, d2), mc.cores = cores)
  failed <- vapply(part, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sum(failed), " folders could not be checked, e.g. ", dirs[which(failed)[1]], ": ",
                        part[[which(failed)[1]]])
  rows <- c(rows, part)
  cat(src, ":", length(dirs), "folders checked\n")
}
rows <- do.call(rbind, rows)
if (is.null(rows)) stop("No input has both an alfak and a karst fit to compare.")
dir.create(file.path(paths$results, "tables"), recursive = TRUE, showWarnings = FALSE)
path <- file.path(paths$results, "tables", "ksets.Rds")
if (file.exists(path)) {
  old <- readRDS(path)
  rows <- rbind(old[!old$source %in% sources, names(rows), drop = FALSE], rows)
}
saveRDS(rows, path)
cat("Inputs compared:", nrow(rows), "\n")
print(table(source = rows$source, design = rows$design, identical = rows$same_k & rows$same_fq & rows$same_nn))
print(table(source = rows$source, method = rows$method))
bad <- rows[!(rows$same_k & rows$same_fq & rows$same_nn), ]
if (nrow(bad)) {
  print(bad)
  stop(nrow(bad), " fits have landscapes whose karyotypes differ from alfak's")
}
