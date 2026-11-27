#!/usr/bin/env Rscript
# Step 0 of PLAN.md: the dispersion statistics of every existing simulated
# dataset (no fits) and criterion S0 (a threshold tau with at most 5% of the
# x1 datasets above it in every source x wavelength cell, and at least 50% of
# the x16 datasets above it).
#   Rscript validation/sparse_joint/dispersion_gate/step0_dispersion.R \
#     --paper-main=DIR --paper-vm=DIR --v7-e=DIR --v7-c2=DIR --own=DIR --out=DIR [--cores=3]
# DIRs hold the paper's raw ABM .Rds files (ALFA-K data/raw/ABM and
# ABM_varyMisseg), the v7 datasets (dev/realdata/noise/v7 sims/{e,c2}/datasets)
# and the own benchmark round's datasets (sparse_two_timepoint_6ae2a24).
# Output in --out: stats.tsv (one row per dataset and n), cells.tsv (share
# above tau per cell and multiple), s0.tsv (tau and verdict per statistic
# and n), step0.pdf.
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "common_gate.R"))
out <- arg("out")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(arg("cores", "3"))
ns <- c(100L, 300L)
beta_default <- 5e-5

rows_of <- function(meta, cols) {
  do.call(rbind, lapply(ns, function(n) {
    s <- gate_stats(cols$c0, cols$c1, n)
    data.frame(meta, n = n, depth_t0 = sum(cols$c0), depth_t1 = sum(cols$c1), t(s), stringsAsFactors = FALSE)
  }))
}
par_rows <- function(files, f) {
  res <- parallel::mclapply(files, function(p) tryCatch(f(p), error = function(e) {
    message("failed: ", p, ": ", conditionMessage(e)); NULL
  }), mc.cores = cores)
  do.call(rbind, res)
}

# paper: inputs w2, w4 and w8 ends (first and last of the last 2, 4, 8
# sampling times before day 120), 1000 cells stored per time point
paper_rows <- function(dir, source) {
  files <- list.files(dir, pattern = "^w_.*[.]Rds$", full.names = TRUE)
  par_rows(files, function(p) {
    id <- sub("[.]Rds$", "", basename(p))
    w <- as.numeric(sub("p", ".", sub("^w_([0-9p]+)_m_.*", "\\1", id)))
    pm <- as.numeric(sub(".*_m_(.*)_rep_.*", "\\1", id))
    x <- as.matrix(readRDS(p)$abm_output$x)
    times <- as.numeric(colnames(x))
    do.call(rbind, lapply(c(2L, 4L, 8L), function(win) {
      tt <- utils::tail(times[times < 120], win)
      cols <- list(c0 = stats::setNames(as.numeric(x[, times == tt[1]]), rownames(x)),
                   c1 = stats::setNames(as.numeric(x[, times == tt[length(tt)]]), rownames(x)))
      rows_of(data.frame(source = source, dataset_id = paste0(id, "_w", win), design = paste0("w", win, "_ends"),
                         wavelength = w, multiple = round(pm / beta_default), interval = tt[length(tt)] - tt[1]), cols)
    }))
  })
}
v7_rows <- function(dir) {
  files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE)
  par_rows(files, function(p) {
    d <- readRDS(p)
    rows_of(data.frame(source = paste0("v7_", d$design), dataset_id = d$dataset_id, design = d$design,
                       wavelength = d$wavelength, multiple = d$multiple, interval = d$interval),
            gate_columns(d$counts))
  })
}
own_rows <- function(dir) {
  files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE)
  par_rows(files, function(p) {
    d <- readRDS(p)
    rows_of(data.frame(source = "own", dataset_id = d$dataset_id, design = d$scenario,
                       wavelength = as.numeric(sub("^w", "", d$wavelength_class)),
                       multiple = round(d$beta / beta_default), interval = d$interval),
            gate_columns(d$counts))
  })
}

started <- Sys.time()
stats <- rbind(paper_rows(arg("paper-main"), "paper"), paper_rows(arg("paper-vm"), "paper_vm"),
               v7_rows(arg("v7-e")), v7_rows(arg("v7-c2")), own_rows(arg("own")))
stats$cell <- paste0(stats$source, " w", stats$wavelength)
utils::write.table(stats, file.path(out, "stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("datasets:", length(unique(stats$dataset_id)), " rows:", nrow(stats), "\n")
print(table(stats$source, stats$multiple) / length(ns))

# S0 per statistic and n
s0 <- list()
cells <- list()
for (st in c("once", "new")) for (n in ns) {
  d <- stats[stats$n == n & is.finite(stats[[st]]), ]
  x1 <- d[d$multiple == 1, ]
  # smallest tau with <= 5% of every cell's x1 datasets above it
  tau <- max(vapply(split(x1[[st]], x1$cell), function(v) {
    v <- sort(v); v[max(1L, ceiling(0.95 * length(v)))]
  }, numeric(1)))
  hi <- d[d$multiple == 16, ]
  s0[[length(s0) + 1L]] <- data.frame(statistic = st, n = n, tau = tau, x1_datasets = nrow(x1),
                                      x1_cells = length(unique(x1$cell)), x16_datasets = nrow(hi),
                                      x16_above = mean(hi[[st]] > tau), x4_above = mean(d[d$multiple == 4, st] > tau),
                                      worst_x1_cell = names(which.max(tapply(x1[[st]], x1$cell, stats::quantile, 0.95))),
                                      pass = mean(hi[[st]] > tau) >= 0.5)
  a <- stats::aggregate(d[[st]], list(cell = d$cell, multiple = d$multiple), function(v) c(n = length(v), median = stats::median(v),
                                                                                         above = mean(v > tau)))
  a <- data.frame(statistic = st, n = n, tau = tau, a[, c("cell", "multiple")], a$x)
  cells[[length(cells) + 1L]] <- a
}
s0 <- do.call(rbind, s0)
cells <- do.call(rbind, cells)
utils::write.table(s0, file.path(out, "s0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(cells, file.path(out, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n== S0\n"); print(s0, row.names = FALSE)
cat("\nStep 0", if (any(s0$pass)) "PASSES" else "FAILS", "\n")
cat("\n== cells (x1 medians and share above tau; x16 share above tau)\n")
options(width = 220)
for (st in c("once", "new")) for (n in ns) {
  a <- cells[cells$statistic == st & cells$n == n, ]
  w <- reshape(a[, c("cell", "multiple", "median", "above")], idvar = "cell", timevar = "multiple", direction = "wide")
  cat("\n--", st, "n =", n, " tau =", signif(a$tau[1], 3), "\n")
  print(w[order(w$cell), ], row.names = FALSE, digits = 3)
}

# figure: statistic by cell and multiple, tau dashed
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  pd <- do.call(rbind, lapply(c("once", "new"), function(st) data.frame(stats[, c("cell", "source", "wavelength", "multiple", "n")],
                                                                        statistic = st, value = stats[[st]])))
  pd <- pd[is.finite(pd$value), ]
  pd$rate <- factor(paste0("×", pd$multiple), levels = paste0("×", c(1, 2, 4, 8, 16)))
  pd$panel <- paste0(pd$statistic, ", n = ", pd$n)
  td <- data.frame(panel = paste0(s0$statistic, ", n = ", s0$n), tau = s0$tau)
  cols <- c("×1" = "#000000", "×2" = "#56B4E9", "×4" = "#009E73", "×8" = "#E69F00", "×16" = "#D55E00")
  g <- ggplot(pd, aes(x = cell, y = value, colour = rate)) +
    geom_boxplot(outlier.size = 0.3, linewidth = 0.3, position = position_dodge2(preserve = "single")) +
    geom_hline(data = td, aes(yintercept = tau), linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") + scale_colour_manual(values = cols, name = "rate") +
    labs(x = NULL, y = "statistic", caption = "Dashed: tau (at most 5% of every cell's ×1 datasets above it). Rates are multiples of 5e-5.") +
    theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")
  ggsave(file.path(out, "step0.pdf"), g, width = 11, height = 11)
}
cat("\nwritten:", out, " (", format(round(difftime(Sys.time(), started, units = "mins"), 1)), ")\n")
