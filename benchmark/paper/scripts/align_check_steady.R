# Check of align_find_steady_state() (R/align.R), which version B's Fig. 6
# runs on the extended lists in place of alfakR::find_steady_state(Nmax = 1):
#
#   ALIGN=karst Rscript scripts/align_check_steady.R [N] [M]
#
# 1. On N Salehi landscapes of the round (alfak's and KARST's, on alfak's
#    lists, short enough for alfakR's builder), at Fig. 6's eight rates: W
#    and the steady state against alfakR's (identical(), largest difference).
# 2. On version B's M longest Salehi landscapes: the time and size of W and
#    of the eight steady states (the cost of Fig. 6 in B).
# Saves align_checks_steady.Rds into ALIGN_RESULTS or <round>/align_karst.

for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(TRUE)
n_check <- if (length(args) >= 1) as.integer(args[1]) else 6L
n_long <- if (length(args) >= 2) as.integer(args[2]) else 2L
paths <- bench_paths()
out <- if (nzchar(Sys.getenv("ALIGN_RESULTS"))) Sys.getenv("ALIGN_RESULTS") else paths$results
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rates <- 0.00002 * 2^(1:8)

## 1. against alfakR on the round's lists
base <- file.path(paths$results_base, "fits", "salehi")
cand <- Sys.glob(file.path(base, "*", "*", "*", "minobs_*", "landscape.Rds"))
set.seed(1)
cand <- sample(cand, min(300L, length(cand)))
size <- vapply(cand, function(f) nrow(readRDS(f)), integer(1))
cand <- cand[size <= 20000]
size <- size[size <= 20000]
set.seed(1)
pick <- cand[unique(c(which.max(size), sample(seq_along(cand), n_check - 1L)))]
same <- do.call(rbind, lapply(pick, function(f) {
  ls <- readRDS(f)
  do.call(rbind, lapply(rates, function(p) {
    W0 <- alfakR:::build_W_rcpp(ls$k, p, 1)
    W1 <- align_W(ls$k, p)
    s0 <- alfakR::find_steady_state(ls, p, Nmax = 1)
    s1 <- align_find_steady_state(ls, p)
    data.frame(landscape = sub(paste0(base, "/"), "", dirname(f), fixed = TRUE), karyotypes = nrow(ls), p = p,
               identical_W = identical(W0, W1), max_diff_W = max(abs(W0 - W1)), identical_steady = identical(s0, s1),
               max_diff_steady = max(abs(s0 - s1)), same_dominant = which.max(s0) == which.max(s1))
  }))
}))
options(width = 200)
cat("against alfakR::find_steady_state on", length(pick), "landscapes x", length(rates), "rates\n")
print(same, row.names = FALSE)

## 2. cost on version B's longest lists
bdir <- file.path(paths$results, "fits", "salehi")
recs <- Sys.glob(file.path(bdir, "*", "*", "*", "task_*.Rds"))
E <- vapply(recs, function(f) readRDS(f)$E, numeric(1))
long <- recs[order(E, decreasing = TRUE)]
long <- long[!duplicated(dirname(dirname(long)))][seq_len(min(n_long, length(long)))]
cost <- do.call(rbind, lapply(long, function(r) {
  f <- Sys.glob(file.path(dirname(r), "minobs_*", "landscape.Rds"))[1]
  ls <- readRDS(f)
  do.call(rbind, lapply(rates, function(p) {
    gc(reset = TRUE)
    t0 <- proc.time()[["elapsed"]]
    W <- align_W(ls$k, p)
    tw <- proc.time()[["elapsed"]] - t0
    s <- align_find_steady_state(ls, p)
    g <- gc()
    data.frame(landscape = sub(paste0(bdir, "/"), "", dirname(f), fixed = TRUE), karyotypes = nrow(ls), p = p,
               nonzero_W = length(W@x), seconds_W = round(tw, 1),
               seconds = round(proc.time()[["elapsed"]] - t0, 1), peak_mb = round(sum(g[, ncol(g)])),
               above_zero = sum(s > 0))
  }))
}))
cat("\ncost on version B's longest lists\n")
print(cost, row.names = FALSE)
saveRDS(list(same = same, cost = cost), file.path(out, "align_checks_steady.Rds"))
