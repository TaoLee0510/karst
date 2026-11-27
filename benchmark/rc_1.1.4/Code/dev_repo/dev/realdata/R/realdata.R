# Helpers of the real-data plan (dev/realdata/PLAN.md). Sourced after
# dev/R/common.R. Patient data are read in place (PANCAN_DATA, PAPER_DATA,
# both bound read-only on HPC); results go under DEV_RESULTS/realdata.

rd_pancan <- Sys.getenv("PANCAN_DATA", file.path(Sys.getenv("RC_HOME"), "data", "pancan"))  # rc freeze
rd_results <- file.path(dev_results, "realdata")
rd_cancers <- c("GBM", "HNSC", "Ovarian", "Melanoma")
rd_diploid <- paste(rep(2L, 22L), collapse = ".")
rd_beta <- 5e-5
rd_rates <- c(0.02, sqrt(0.02 * 1), 1)
rd_central <- 2L

rd_load_karst <- function() {
  lib <- Sys.getenv("DEV_LIBDIR", "")
  if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(karst))
  invisible(as.character(utils::packageVersion("karst")))
}

# Every NUMBAT-derived input as PANcanKFLs delivers it to alfak.
rd_numbat_inputs <- function() {
  do.call(rbind, lapply(rd_cancers, function(cancer) {
    dir <- file.path(rd_pancan, "processed", "ALFA_K", cancer, "ALFA-K_inputs")
    files <- list.files(dir, pattern = "[.]Rds$", full.names = TRUE)
    data.frame(cancer = cancer, sample = sub("[.]Rds$", "", basename(files)), file = files,
               stringsAsFactors = FALSE)
  }))
}

# The delivered input with the all-diploid karyotype put back. Each delivered
# column is matched to the karyo_in column with the same non-diploid counts.
rd_restore <- function(cancer, sample) {
  file <- file.path(rd_pancan, "processed", "ALFA_K", cancer, "ALFA-K_inputs", paste0(sample, ".Rds"))
  x <- as.matrix(readRDS(file)$x)
  storage.mode(x) <- "integer"
  robj <- file.path(rd_pancan, "processed", "ALFA_K", cancer, "karyotype_count", paste0(sample, ".RObj"))
  out <- list(cancer = cancer, sample = sample, delivered = x, times = as.numeric(colnames(x)),
              status = "ok", matched = NA_character_, diploid = NA_integer_)
  if (!file.exists(robj)) {
    out$status <- "no karyo_in"
    return(out)
  }
  e <- new.env()
  load(robj, envir = e)
  k <- as.matrix(e$karyo_in)
  nd <- k[rownames(k) != rd_diploid, , drop = FALSE]
  if (!all(rownames(x) %in% rownames(nd))) {
    out$status <- "rows differ"
    return(out)
  }
  extra <- setdiff(rownames(nd), rownames(x))
  matched <- vapply(seq_len(ncol(x)), function(j) {
    hit <- which(vapply(seq_len(ncol(nd)), function(c) {
      all(nd[rownames(x), c] == x[, j]) && all(nd[extra, c] == 0)
    }, logical(1)))
    if (length(hit) == 1L) hit else NA_integer_
  }, integer(1))
  if (anyNA(matched) || anyDuplicated(matched)) {
    out$status <- "columns not matched"
    return(out)
  }
  out$matched <- colnames(k)[matched]
  diploid <- if (rd_diploid %in% rownames(k)) as.integer(k[rd_diploid, matched]) else rep(0L, ncol(x))
  out$diploid <- diploid
  restored <- rbind(x, matrix(diploid, nrow = 1L, dimnames = list(rd_diploid, colnames(x))))
  storage.mode(restored) <- "integer"
  out$restored <- restored
  out
}

# First and last time point, karyotypes with no cell there removed.
rd_ends <- function(counts) {
  counts <- counts[, c(1L, ncol(counts)), drop = FALSE]
  counts[rowSums(counts) > 0, , drop = FALSE]
}

# Copy-number bounds of a cancer type (every input of it, restored), as
# real_cohort.R takes them over its input folder: min(1, lowest), highest + 1.
rd_cn_bounds <- function(cancer) {
  inputs <- rd_numbat_inputs()
  labels <- unlist(lapply(inputs$sample[inputs$cancer == cancer], function(s) {
    rownames(readRDS(file.path(rd_pancan, "processed", "ALFA_K", cancer, "ALFA-K_inputs", paste0(s, ".Rds")))$x)
  }))
  cn <- as.integer(unlist(strsplit(unique(c(labels, rd_diploid)), ".", fixed = TRUE)))
  c(min_cn = min(1L, min(cn)), max_cn = max(cn) + 1L)
}

# The compact record of a fit kept for the analyses (the fit object itself is
# not kept).
rd_fit_record <- function(fit, seconds, extra = list()) {
  if (inherits(fit, "error")) {
    return(c(list(status = "failed", condition_class = class(fit)[[1L]], message = conditionMessage(fit),
                  seconds = seconds), extra))
  }
  s <- fit$summary
  keep <- intersect(c("karyotype", "major", "shell", "count_t0", "count_t1", "expected_t0", "expected_t1",
                      "fitness_mean", "relative_fitness", "relative_sd", "relative_sd_model"), names(s))
  c(list(status = "ok", seconds = seconds, summary = s[keep],
         dropped = fit$dropped,
         diagnostics = fit$diagnostics[intersect(c("domain_coverage", "history", "history_profile", "convergence",
                                                   "posterior_objective", "objective", "n_states", "n_major",
                                                   "dropped_karyotypes", "initial_state", "beta"),
                                                 names(fit$diagnostics))]),
    extra)
}

# Pearson dispersion of the counts at the second time point (real_cohort.R).
rd_dispersion <- function(s) {
  seen <- s$count_t1 > 0 | s$expected_t1 >= 1
  if (!any(seen)) return(NA_real_)
  sum(((s$count_t1 - s$expected_t1)^2 / pmax(s$expected_t1, 1e-8))[seen]) / sum(seen)
}

# Non-major domain karyotypes seen at the second time point, over the number
# Poisson(expected) predicts.
rd_zero_ratio <- function(s) {
  minor <- !s$major
  expected <- sum(1 - exp(-s$expected_t1[minor]))
  if (expected <= 0) return(NA_real_)
  sum(s$count_t1[minor] > 0) / expected
}

# Agreement of two fits on the karyotypes both observed. The difference is
# taken after removing its mean over the karyotypes with at least 10 cells
# (the two fits can have different references).
rd_compare <- function(a, b, majors_only = FALSE) {
  a$total <- a$count_t0 + a$count_t1
  at <- match(a$karyotype, b$karyotype)
  use <- !is.na(at) & a$total > 0
  if (majors_only) use <- use & a$major & b$major[ifelse(is.na(at), 1L, at)]
  if (sum(use) < 3L) return(data.frame(n = sum(use), spearman = NA_real_, share_over_1sd = NA_real_))
  fa <- a$relative_fitness[use]; fb <- b$relative_fitness[at[use]]
  d <- fa - fb
  anchor <- a$total[use] >= 10
  d <- d - if (any(anchor)) mean(d[anchor]) else mean(d)
  sd <- sqrt(a$relative_sd[use]^2 + b$relative_sd[at[use]]^2)
  data.frame(n = sum(use), spearman = suppressWarnings(stats::cor(fa, fb, method = "spearman")),
             share_over_1sd = mean(abs(d) > sd))
}

# Two random halves of the cells at each time point (real_cohort.R).
rd_split_half <- function(counts, seed) {
  set.seed(seed)
  halves <- lapply(seq_len(ncol(counts)), function(column) {
    cells <- rep(seq_len(nrow(counts)), counts[, column])
    first <- sample(length(cells), length(cells) %/% 2L)
    cbind(tabulate(cells[first], nrow(counts)), tabulate(cells[-first], nrow(counts)))
  })
  lapply(1:2, function(half) {
    out <- do.call(cbind, lapply(halves, function(h) h[, half]))
    dimnames(out) <- dimnames(counts)
    storage.mode(out) <- "integer"
    out[rowSums(out) > 0, , drop = FALSE]
  })
}

rd_seed <- function(text) sum(utf8ToInt(text) * seq_along(utf8ToInt(text))) %% 100000L + 20261119L
