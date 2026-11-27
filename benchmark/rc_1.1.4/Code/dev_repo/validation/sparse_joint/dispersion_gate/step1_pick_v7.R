#!/usr/bin/env Rscript
# Step 1 of PLAN.md, pick step: the cell-half hold-out of one v7 dataset at
# one multiple of 5e-5. Each time point's cells are split at random into
# halves A and B (seed from the dataset); KARST (library 1.1.3) is fitted to
# A as in step1_fit_v7.R, and B is scored by the multinomial negative
# log-likelihood of its counts over the fit's domain states, with the fit's
# expected counts normalised within the domain (each time point; B cells
# outside the domain are left out and counted).
#   Rscript validation/sparse_joint/dispersion_gate/step1_pick_v7.R TASK_FILE ID STAGE V7_RESULTS OUT
# Output: OUT/x<m>/<dataset_id>.rds.
args <- commandArgs(TRUE)
task_file <- args[[1]]; index <- args[[2]]; stage <- args[[3]]; v7_results <- args[[4]]; out <- args[[5]]
Sys.setenv(DEV_REPO = stage, DEV_RESULTS = v7_results, OMP_NUM_THREADS = "1")
source(file.path(stage, "dev", "R", "common.R"))
task <- dev_task_row(task_file, index)
loadNamespace("karst")
m <- as.numeric(task$multiple); beta <- 5e-5 * m
out_file <- file.path(out, sprintf("x%02d", as.integer(m)), paste0(gsub("[^A-Za-z0-9_.-]", "_", task$dataset_id), ".rds"))
if (file.exists(out_file)) quit(save = "no")
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
d <- readRDS(v7_path("sims", task$family, "datasets", paste0(task$dataset_id, ".rds")))
key <- list(family = task$family, dataset_id = task$dataset_id, multiple = m, beta = beta, truth_multiple = d$multiple,
            karst = as.character(utils::packageVersion("karst")))
counts <- as.matrix(d$counts)
# split every time point's cells into halves (B gets the odd one)
set.seed(d$seed + 7L)
half <- function(v) {
  cells <- rep.int(seq_along(v), v)
  a <- sample(cells, floor(length(cells) / 2))
  tabulate(a, nbins = length(v))
}
A <- counts
A[, 1] <- half(counts[, 1]); A[, 2] <- half(counts[, 2])
B <- counts - A
keepA <- rowSums(A) > 0
A <- A[keepA, , drop = FALSE]
max_cn <- if (d$family == "e") 8L else as.integer(d$max_cn)
arguments <- if (d$family == "e") list(A, times = c(0, d$interval), beta = beta, birth_rate = 1, min_cn = 1L, max_cn = 8L) else
  list(A, dt = d$interval, beta = beta, birth_rate = 1, min_cn = 1L, max_cn = max_cn)
set.seed(d$seed)
started <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(do.call(karst::fit_sparse_joint, arguments)), error = function(e) e)
seconds <- proc.time()[["elapsed"]] - started
if (inherits(fit, "error")) {
  dev_save_rds(c(key, list(status = "failed", message = conditionMessage(fit), seconds = seconds)), out_file)
  quit(save = "no")
}
s <- fit$summary
expected <- grep("^expected_", names(s), value = TRUE)[1:2]
nll <- 0; inside <- 0; total <- sum(B)
for (j in 1:2) {
  b <- stats::setNames(as.numeric(B[, j]), rownames(B))
  b <- b[b > 0]
  e <- stats::setNames(as.numeric(s[[expected[j]]]), s$karyotype)
  e <- e[is.finite(e) & e > 0]
  p <- e / sum(e)
  hit <- names(b)[names(b) %in% names(p)]
  inside <- inside + sum(b[hit])
  nll <- nll - sum(b[hit] * log(p[hit]))
}
dev_save_rds(c(key, list(status = "ok", seconds = seconds, nll_half = nll, b_cells = total, b_inside = inside,
                         domain = s$karyotype, objective = fit$diagnostics$objective)), out_file)
