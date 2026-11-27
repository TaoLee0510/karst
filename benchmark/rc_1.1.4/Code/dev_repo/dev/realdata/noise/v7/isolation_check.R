#!/usr/bin/env Rscript
# Step 0 of dev/realdata/noise/PLAN_v7.md, isolation.
#   mode "prior": the 20 new landscapes against every earlier landscape. Seeds
#     against the landscape tables of the repo and dev code; centroid matrices
#     against every earlier truth file (read only): files under a /truth/
#     folder in the karst_dev results (v7 excluded), the karst clone, the
#     alfak2-era folders, and the paper's raw simulations (true_landscape).
#   mode "own": every v7 truth file (e sets and C2) has the planned centroids
#     of its landscape, so both families share the landscapes.
#   Rscript dev/realdata/noise/v7/isolation_check.R TASK_FILE ID   (column: mode)
# Output: DEV_RESULTS/realdata/noise/v7/isolation_<mode>.rds and a printed verdict.
source(file.path(Sys.getenv("DEV_REPO", getwd()), "dev", "R", "common.R"))
args <- commandArgs(TRUE)
task <- dev_task_row(args[[1]], args[[2]])
library(karst)
dev_source_validation()
source(file.path(dev_repo, "dev", "realdata", "noise", "v7", "common_v7.R"))
stt_dir <- file.path(dev_repo, "benchmark", "sparse_two_timepoint")
for (part in c("landscape.R", "config.R", "generate.R")) source(file.path(stt_dir, part))
planned <- lapply(seq_len(nrow(v7_landscapes)), function(i) {
  set.seed(v7_landscapes$landscape_seed[i])
  stt_generate_centroids(rep(2L, 22L), 10L, v7_landscapes$wavelength[i])
})
names(planned) <- v7_landscapes$landscape_id
same_matrix <- function(a, b) is.matrix(a) && is.matrix(b) && all(dim(a) == dim(b)) && max(abs(a - b)) < 1e-9
match_planned <- function(mat) {
  if (!is.matrix(mat) || !is.numeric(mat)) return(character())
  names(planned)[vapply(planned, function(p) same_matrix(p, mat), logical(1))]
}

if (task$mode == "prior") {
  ## seeds of the landscape tables in code
  tables <- list(stt_dev_and_heldout = stt_config(tempfile())$landscapes,
                 mt_series = mt_config(tempfile(), "series")$landscapes,
                 mt_heldout = mt_config(tempfile(), "heldout")$landscapes)
  code_seeds <- unlist(lapply(tables, `[[`, "landscape_seed"))
  extra <- c(63000L + 1:20, 64000L + 1:12, 93000L + 1:200, 94000L + 1:200, 95000L + 1:200, 96000L + 1:200,
             99000L + 1:200, 101000L + 1:200, 102000L + 1:200, 103000L + 1:200, 111000L + 1:12, 112000L + 1:12)
  seed_clash <- intersect(v7_landscapes$landscape_seed, c(code_seeds, extra))
  cat("code landscape seeds checked:", length(unique(c(code_seeds, extra))), " clashes:", length(seed_clash), "\n")
  ## centroid matrices of earlier truth files
  roots <- c("/share/lab_crd/taoli/Project/karst_dev/results", "/share/lab_crd/taoli/Project/karst",
             "/share/lab_crd/taoli/Project/alfak2", "/share/lab_crd/taoli/Project/alfak2_results_backup_2026-11-17",
             "/share/lab_crd/taoli/Project/alfak2_results_backup_2026-11-18_extra",
             "/share/lab_crd/taoli/Project/alfak2-multi-timepoint", "/share/lab_crd/taoli/Project/alfak2-paper-benchmark",
             "/share/lab_crd/taoli/Project/alfak2-salehi-variants",
             "/share/lab_crd/taoli/Project/karst_history_choice_untracked_2026-11-19")
  roots <- roots[dir.exists(roots)]
  per_root <- lapply(roots, function(r) {
    out <- system2("find", c(shQuote(r), "-path", shQuote("*/truth/*"), shQuote("("), "-name", shQuote("*.rds"),
                             "-o", "-name", shQuote("*.Rds"), shQuote(")")), stdout = TRUE)
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) stop("find failed on ", r, " (status ", status, ")")
    out
  })
  names(per_root) <- roots
  cat("truth files per root:\n"); print(vapply(per_root, length, integer(1)))
  files <- unlist(per_root)
  files <- files[!grepl("/realdata/noise/v7/", files)]
  if (length(files) < 100L) stop("too few earlier truth files found (", length(files), "); the search is wrong")
  paper <- Sys.getenv("PAPER_DATA", "")
  if (nzchar(paper)) files <- c(files, list.files(file.path(paper, "raw"), pattern = "[.][Rr]ds$", recursive = TRUE,
                                                  full.names = TRUE))
  files <- unique(files)
  cat("earlier files to read:", length(files), "\n")
  hits <- list(); unread <- character(); n_matrices <- 0L
  for (f in files) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) { unread <- c(unread, f); next }
    mats <- if (is.list(x)) list(x$centroids, x$truth$centroids, x$true_landscape) else list()
    for (mat in mats) {
      if (!is.matrix(mat)) next
      n_matrices <- n_matrices + 1L
      h <- match_planned(mat)
      if (length(h)) hits[[length(hits) + 1L]] <- data.frame(file = f, landscape = h)
    }
    rm(x); invisible(gc(FALSE))
  }
  hits <- if (length(hits)) do.call(rbind, hits) else data.frame(file = character(), landscape = character())
  cat("centroid matrices compared:", n_matrices, " unreadable files:", length(unread), " matches:", nrow(hits), "\n")
  verdict <- if (!length(seed_clash) && !nrow(hits)) "isolated" else "NOT isolated"
  dev_save_rds(list(seed_clash = seed_clash, files = length(files), matrices = n_matrices, unread = unread,
                    hits = hits, roots = roots, verdict = verdict), v7_path("isolation_prior.rds"))
} else if (task$mode == "own") {
  tr <- c(list.files(v7_path("sims", "e"), pattern = "_P01[.]rds$", recursive = TRUE, full.names = TRUE),
          list.files(v7_path("sims", "c2", "grid"), pattern = "[.]rds$", recursive = TRUE, full.names = TRUE))
  res <- do.call(rbind, lapply(tr, function(f) {
    x <- readRDS(f)
    data.frame(file = f, landscape = x$landscape_id, ok = same_matrix(x$centroids, planned[[x$landscape_id]]),
               multiple = x$v7$multiple)
  }))
  cat("v7 truth files:", nrow(res), " with the planned centroids:", sum(res$ok), "\n")
  verdict <- if (nrow(res) == 180L && all(res$ok)) "shared and as planned" else "MISMATCH"
  dev_save_rds(list(table = res, verdict = verdict), v7_path("isolation_own.rds"))
} else stop("unknown mode")
cat("VERDICT:", verdict, "\n")
