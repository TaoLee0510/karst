# The paper's Fig. 2d on the Salehi passages, raw and with the sampling noise
# of the passages held fixed (from round 827c4b5 on).
#
# Fig. 2d ranks, across trajectories, the observed change from the parent
# passage to a daughter (Wasserstein distance) by the predicted change (parent
# to forecast). With a few hundred cells per passage much of the observed
# change is sampling noise: splitting a passage's cells at random into two
# halves, the distance between the halves ranks the observed change better
# than any forecast (karst validation/sparse_joint/cloud_level, step 2). This
# script reports, per arm, Spearman's rho of observed and predicted change and
# the partial rho with the noise held fixed, with lineage-type bootstrap
# intervals and differences from alfak, on the rows every arm of a design
# has.
#
#   Rscript scripts/fig2d_noise.R --stage=chunk [--index=N]   trajectories 10(N-1)+1 .. 10N
#   Rscript scripts/fig2d_noise.R --stage=analyse
#
# chunk (index from --index or SLURM_ARRAY_TASK_ID; 28 chunks): runs
#   scripts/fig2d_check.R --stage=salehi for the "ends" design (karst,
#   alfak) and the "all" design (karst, karst_pairwise, karst_state_space,
#   alfak), which writes one row per fit, minobs and daughter passage (base =
#   observed change, pred0 = predicted change, on the landscape's karyotypes,
#   at the paper's evaluation day), then the noise of each row: the mean
#   distance between two random halves (20 splits) of the parent passage and
#   of the daughter passage, on the karyotypes of alfak's landscape of the
#   design (the axis of Fig. 2d), and of the parent on all its cells.
# analyse: BENCH_RESULTS/analysis/fig2d_noise/{fig2d_noise.tsv,
#   fig2d_noise_differences.tsv, figures/fig2d_noise.pdf}.

for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
args <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
stage <- arg("stage")
root <- file.path(paths$results, "analysis", "fig2d_noise")
designs <- list(ends = c("karst", "alfak"), all = c("karst", "karst_pairwise", "karst_state_space", "alfak"))

if (stage == "chunk") {
  index <- as.integer(arg("index", Sys.getenv("SLURM_ARRAY_TASK_ID")))
  # rows of every arm (fig2d_check.R reads the fits, refits nothing)
  for (design in names(designs)) {
    dir <- file.path(root, paste0("rows_", design))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    env <- c(paste0("FIG2D_DESIGN=", design), paste0("FIG2D_METHODS=", paste(designs[[design]], collapse = ",")),
             paste0("FIG2D_OUT=", dir))
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("scripts/fig2d_check.R", "--stage=salehi", paste0("--index=", index)), env = env)
    if (!identical(status, 0L)) stop("fig2d_check.R failed for design ", design, " (status ", status, ")")
  }
  # sampling noise of each row's parent and daughter passage
  meta <- utils::read.csv(paths$salehi_meta)
  inputs <- list.files(paths$salehi_inputs)
  s2m <- function(k) do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
  wd <- function(a, b) {
    a <- a[a > 0]
    b <- b[b > 0]
    if (!length(a) || !length(b)) return(NA_real_)
    transport::wasserstein(transport::wpp(s2m(names(a)), mass = a / sum(a)),
                           transport::wpp(s2m(names(b)), mass = b / sum(b)))
  }
  halves <- function(p, keep = NULL, reps = 20L) {
    if (!is.null(keep)) p <- p[names(p) %in% keep]
    if (sum(p) < 4) return(NA_real_)
    cells <- rep(names(p), p)
    mean(vapply(seq_len(reps), function(r) {
      set.seed(r)
      i <- sample.int(length(cells), length(cells) %/% 2L)
      wd(table(cells[i]), table(cells[-i]))
    }, numeric(1)))
  }
  noise <- list()
  for (design in names(designs)) {
    chunk_file <- file.path(root, paste0("rows_", design), sprintf("chunk_%03d.Rds", index))
    rows <- if (file.exists(chunk_file)) readRDS(chunk_file) else NULL
    if (is.null(rows) || !nrow(rows)) next
    keys <- unique(rows[c("fi", "min_obs", "dj")])
    for (j in seq_len(nrow(keys))) {
      k <- keys[j, ]
      ls_path <- file.path(paths$results, "fits", "salehi", k$fi, design, "alfak", paste0("minobs_", k$min_obs),
                           "landscape.Rds")
      fj <- paste(meta[meta$uid == k$dj, c("datasetname", "timepoint")], collapse = "_")
      hit <- inputs[grepl(fj, inputs)]
      if (!file.exists(ls_path) || !length(hit)) next
      lk <- readRDS(ls_path)$k
      tru <- readRDS(file.path(paths$salehi_inputs, hit[1]))$x
      p <- stats::setNames(tru[, ncol(tru) - 1L], rownames(tru))
      d <- stats::setNames(tru[, ncol(tru)], rownames(tru))
      p <- p[p > 0]
      d <- d[d > 0]
      noise[[length(noise) + 1L]] <- data.frame(design = design, k, cells_parent = sum(p), cells_daughter = sum(d),
                                                noise_parent = halves(p, lk), noise_daughter = halves(d, lk),
                                                noise_parent_all = halves(p), stringsAsFactors = FALSE)
    }
  }
  dir.create(file.path(root, "noise"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(do.call(rbind, noise), file.path(root, "noise", sprintf("chunk_%03d.Rds", index)))
  cat("chunk", index, "done:", length(noise), "noise rows\n")
}

if (stage == "analyse") {
  suppressPackageStartupMessages(library(ggplot2))
  read_chunks <- function(dir) {
    f <- list.files(dir, pattern = "^chunk_.*[.]Rds$", full.names = TRUE)
    bind_rows_fill(lapply(f, readRDS))
  }
  noise <- read_chunks(file.path(root, "noise"))
  noise$key <- paste(noise$design, noise$fi, noise$min_obs, noise$dj)
  noise$noise <- (noise$noise_parent + noise$noise_daughter) / 2
  sp <- function(a, b) {
    k <- is.finite(a) & is.finite(b)
    if (sum(k) < 5) return(NA_real_)
    stats::cor(a[k], b[k], method = "spearman")
  }
  partial <- function(y, a, z) {
    ok <- is.finite(y) & is.finite(a) & is.finite(z)
    if (sum(ok) < 5) return(NA_real_)
    ry <- stats::residuals(stats::lm(rank(y[ok]) ~ rank(z[ok])))
    ra <- stats::residuals(stats::lm(rank(a[ok]) ~ rank(z[ok])))
    stats::cor(ry, ra)
  }
  sets <- list()
  for (design in names(designs)) {
    rows <- read_chunks(file.path(root, paste0("rows_", design)))
    if (is.null(rows)) next
    rows$key <- paste(design, rows$fi, rows$min_obs, rows$dj)
    arms <- intersect(designs[[design]], unique(rows$method))
    ok_keys <- Reduce(intersect, lapply(arms, function(m) rows$key[rows$method == m & is.finite(rows$pred0)]))
    ref <- rows[rows$method == "alfak" & rows$key %in% ok_keys, ]
    w <- data.frame(key = ref$key, fi = ref$fi, type = ref$type, base = ref$base, stringsAsFactors = FALSE)
    w$passages <- as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", w$fi))
    for (m in arms) w[[m]] <- rows$pred0[rows$method == m][match(w$key, rows$key[rows$method == m])]
    w$noise <- noise$noise[match(w$key, noise$key)]
    w <- w[is.finite(w$base) & is.finite(w$noise), ]
    group <- if (design == "ends") "2 time points" else "all time points"
    if (design == "ends") {
      sets[[length(sets) + 1L]] <- list(group = group, subset = "every trajectory", w = w, arms = arms)
      sets[[length(sets) + 1L]] <- list(group = group, subset = "one interval (2 passages)",
                                        w = w[w$passages == 2L, ], arms = arms)
    } else {
      sets[[length(sets) + 1L]] <- list(group = group, subset = "3-10 passages", w = w[w$passages > 2L, ], arms = arms)
    }
  }
  draws <- 4000L
  est <- list()
  dif <- list()
  for (s in sets) {
    w <- s$w
    g <- split(seq_len(nrow(w)), w$type)
    stat <- function(ix) {
      x <- w[ix, ]
      c(unlist(lapply(s$arms, function(m) c(sp(x$base, x[[m]]), partial(x$base, x[[m]], x$noise)))), sp(x$base, x$noise))
    }
    full <- stat(seq_len(nrow(w)))
    set.seed(1)
    boot <- vapply(seq_len(draws), function(i) stat(unlist(g[sample.int(length(g), replace = TRUE)], use.names = FALSE)),
                   numeric(length(full)))
    ci <- apply(boot, 1L, stats::quantile, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    lab <- c(as.vector(t(outer(s$arms, c("raw", "noise held fixed"), paste, sep = "|"))), "noise|raw")
    for (i in seq_along(lab)) {
      parts <- strsplit(lab[i], "|", fixed = TRUE)[[1]]
      est[[length(est) + 1L]] <- data.frame(group = s$group, subset = s$subset, rows = nrow(w),
                                            trajectories = length(unique(w$fi)), lineage_types = length(g),
                                            arm = parts[1], measure = parts[2], rho = round(full[i], 3),
                                            lo = round(ci[1, i], 3), hi = round(ci[2, i], 3), stringsAsFactors = FALSE)
    }
    ia <- match("alfak|raw", lab)
    ip <- match("alfak|noise held fixed", lab)
    for (m in setdiff(s$arms, "alfak")) for (meas in c("raw", "noise held fixed")) {
      i <- match(paste(m, meas, sep = "|"), lab)
      j <- if (meas == "raw") ia else ip
      dv <- boot[i, ] - boot[j, ]
      dif[[length(dif) + 1L]] <- data.frame(group = s$group, subset = s$subset, arm = m, measure = meas,
                                            difference_from_alfak = round(full[i] - full[j], 3),
                                            lo = round(stats::quantile(dv, 0.025, na.rm = TRUE, names = FALSE), 3),
                                            hi = round(stats::quantile(dv, 0.975, na.rm = TRUE, names = FALSE), 3),
                                            stringsAsFactors = FALSE)
    }
  }
  est <- do.call(rbind, est)
  dif <- do.call(rbind, dif)
  utils::write.table(est, file.path(root, "fig2d_noise.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(dif, file.path(root, "fig2d_noise_differences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  options(width = 200)
  cat("== Fig. 2d, Spearman of observed and predicted change (rows every arm has; lineage-type bootstrap)\n")
  print(est, row.names = FALSE)
  cat("\n== Difference from alfak\n")
  print(dif, row.names = FALSE)

  # figure: one panel per design and subset, both measures, every arm
  arm_fit <- function(arm, group) {
    lab <- ifelse(arm == "noise", "sampling noise", paste(unname(method_label[arm]), group))
    lab
  }
  pe <- est[est$arm != "noise", ]
  pe$fit <- factor(arm_fit(pe$arm, pe$group), fit_levels)
  pe$panel <- paste0(pe$group, "\n", pe$subset, " (", pe$rows, " rows)")
  pe$measure <- factor(pe$measure, c("raw", "noise held fixed"))
  pn <- est[est$arm == "noise", ]
  pn$panel <- paste0(pn$group, "\n", pn$subset, " (", pn$rows, " rows)")
  p <- ggplot(pe, aes(x = measure, y = rho, colour = fit)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_hline(data = pn, aes(yintercept = rho), linetype = 2, colour = "grey40") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6), size = 0.25) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = fit_colours, drop = TRUE) +
    labs(x = NULL, y = "Spearman, observed vs predicted change\n(raw; partial, sampling noise held fixed)",
         colour = NULL, caption = paste("Dashed: Spearman of the observed change with the sampling noise of the passages",
                                        "(random halves). Bars: 95% lineage-type bootstrap.")) +
    theme_classic(base_size = 8) + theme(strip.background = element_blank(), legend.position = "top")
  dir.create(file.path(root, "figures"), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(root, "figures", "fig2d_noise.pdf"), p, width = 180, height = 80, units = "mm")
}
