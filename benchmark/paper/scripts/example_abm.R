# The paper's Supp. Fig. 2c-e (S05_plot_abm_results.R): one example
# simulation, w_1p6_m_5e-05_rep_74 (lambda = 1.6), fitted by every method.
#   2c  true against estimated fitness of the frequent karyotypes (both
#       centred), per fit and minobs, windows of 2 and 8 passages
#   2d-e  frequencies of the 8 most common karyotypes over the window:
#       observed, and projected from each bootstrap replicate's frequencies
#       at time 0 and fitness (minobs 20), for the 8- and 2-passage windows
#
#   Rscript scripts/example_abm.R
#
# Writes BENCH_RESULTS/analysis/abm/supp2c_example.tsv (per fit and minobs:
# frequent karyotypes, Pearson r of true and estimated fitness),
# supp2de_example.tsv (per fit: mean absolute error of the median projection
# at the observed times) and figures/supp2c_example.png, supp2de_example.png.

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
out <- file.path(paths$results, "analysis", "abm")
fig <- file.path(paths$results, "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
id <- "w_1p6_m_5e-05_rep_74"
highlight <- c("2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.1.2.2", "2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.3.2.2.2")
sim <- read_abm_sim(id, paths$abm_raw)
w <- abm_sim_info(id)$w

# The fits: label, window, folder of minobs m.
fits <- do.call(rbind, lapply(c(2L, 8L), function(win) {
  ends <- paste0("w", win, "_ends")
  all <- if (win == 2L) ends else paste0("w", win, "_all")
  rbind(data.frame(fit = "KARST 2 time points", window = win, design = ends, method = "karst"),
        data.frame(fit = paste(method_label[c("karst", karst_other_methods)], "all time points"), window = win,
                   design = all, method = if (win == 2L) "karst" else c("karst", karst_other_methods)),
        data.frame(fit = "alfak dev2@f93b2d2 2 time points", window = win, design = ends, method = "alfak"),
        data.frame(fit = "alfak dev2@f93b2d2 all time points", window = win, design = all, method = "alfak"),
        data.frame(fit = "paper (alfak main@4ad82b0)", window = win, design = NA, method = "alfak_paper"))
}))
fit_dir <- function(i, m) {
  if (fits$method[i] == "alfak_paper") {
    return(file.path(paths$paper_abm, id, sprintf("minobs_%s_ntp_%d", m, fits$window[i])))
  }
  file.path(paths$results, "fits", "abm", id, fits$design[i], fits$method[i], paste0("minobs_", m))
}

## 2c ---------------------------------------------------------------------------------
c_rows <- do.call(rbind, lapply(seq_len(nrow(fits)), function(i) do.call(rbind, lapply(c(5, 10, 20), function(m) {
  path <- file.path(fit_dir(i, m), "landscape.Rds")
  if (!file.exists(path)) return(NULL)
  ls <- readRDS(path)
  ls <- ls[ls$fq, ]
  if (!nrow(ls)) return(NULL)
  tru <- grf_fitness(ls$k, w, sim$true_landscape)
  data.frame(fit = fits$fit[i], window = fits$window[i], minobs = m, k = ls$k, est = ls$mean - mean(ls$mean),
             true = tru - mean(tru), hi = ls$k %in% highlight)
}))))
c_tab <- do.call(rbind, lapply(split(c_rows, c_rows[c("fit", "window", "minobs")], drop = TRUE), function(x) data.frame(
  fit = x$fit[1], window = x$window[1], minobs = x$minobs[1], frequent = nrow(x),
  pearson_r = if (nrow(x) > 2) round(stats::cor(x$true, x$est), 3) else NA_real_)))
c_tab <- c_tab[order(c_tab$window, match(c_tab$fit, fit_levels), c_tab$minobs), ]
utils::write.table(c_tab, file.path(out, "supp2c_example.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
c_rows$fit <- factor(c_rows$fit, fit_levels)
p <- ggplot(c_rows, aes(true, est, colour = hi)) +
  facet_grid(fit + paste("window", window) ~ paste("minobs", minobs), scales = "free") +
  geom_abline(linewidth = 0.2) + geom_point(size = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "#FFA500")) +
  labs(x = "true fitness (centred)", y = "estimated fitness (centred)", subtitle = paste(id, "(lambda = 1.6)")) +
  theme_classic(base_size = 7) + theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0))
ggsave(file.path(fig, "supp2c_example.png"), p, width = 170, height = 22 * length(unique(paste(c_rows$fit, c_rows$window))) + 20,
       units = "mm", dpi = 200, bg = "white", limitsize = FALSE)

## 2d-e -------------------------------------------------------------------------------
project <- function(x0, f, tt) {
  lv <- outer(log(x0), rep(1, length(tt))) + outer(f, tt)
  lv <- sweep(lv, 2, apply(lv, 2, max))
  sweep(exp(lv), 2, colSums(exp(lv)), "/")
}
de_rows <- list()
de_tab <- list()
for (win in c(8L, 2L)) {
  yi <- abm_input(sim, abm_designs[abm_designs$design == if (win == 2L) "w2_ends" else "w8_all", ])
  times <- input_times(yi)
  obs <- yi$x[order(rowSums(yi$x), decreasing = TRUE), , drop = FALSE]
  top <- utils::head(rownames(obs), 8)
  obs_f <- sweep(as.matrix(yi$x), 2, colSums(yi$x), "/")[top, , drop = FALSE]
  tt <- seq(min(times), max(times), length.out = 50)
  for (i in which(fits$window == win)) {
    path <- file.path(fit_dir(i, 20), "bootstrap_res.Rds")
    if (!file.exists(path)) next
    b <- readRDS(path)
    reps <- lapply(seq_len(nrow(b$final_frequencies)), function(r) {
      x0 <- b$final_frequencies[r, ]
      ok <- x0 > 0
      pr <- project(x0[ok], b$final_fitness[r, ok], tt)
      rownames(pr) <- colnames(b$final_frequencies)[ok]
      pr[intersect(top, rownames(pr)), , drop = FALSE]
    })
    keep <- Reduce(intersect, lapply(reps, rownames))
    if (!length(keep)) next
    arr <- simplify2array(lapply(reps, function(r) r[keep, , drop = FALSE] / rep(colSums(r[keep, , drop = FALSE]),
                                                                                 each = length(keep))))
    med <- apply(arr, c(1, 2), stats::median)
    at <- vapply(times, function(t) which.min(abs(tt - t)), integer(1))
    obs_k <- sweep(obs_f[keep, , drop = FALSE], 2, colSums(obs_f[keep, , drop = FALSE]), "/")
    de_tab[[length(de_tab) + 1]] <- data.frame(fit = fits$fit[i], window = win, karyotypes = length(keep),
                                               mae = round(mean(abs(med[, at] - obs_k)), 4))
    de_rows[[length(de_rows) + 1]] <- data.frame(fit = fits$fit[i], window = win, k = rep(keep, length(tt)),
                                                 time = rep(tt, each = length(keep)), freq = as.vector(med))
  }
}
de_tab <- do.call(rbind, de_tab)
utils::write.table(de_tab, file.path(out, "supp2de_example.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
de <- do.call(rbind, de_rows)
de$fit <- factor(de$fit, fit_levels)
obs_long <- do.call(rbind, lapply(c(8L, 2L), function(win) {
  yi <- abm_input(sim, abm_designs[abm_designs$design == if (win == 2L) "w2_ends" else "w8_all", ])
  f <- sweep(as.matrix(yi$x), 2, colSums(yi$x), "/")
  top <- utils::head(rownames(yi$x)[order(rowSums(yi$x), decreasing = TRUE)], 8)
  data.frame(window = win, k = rep(top, ncol(f)), time = rep(input_times(yi), each = length(top)),
             freq = as.vector(f[top, ]))
}))
lab <- stats::setNames(LETTERS[seq_along(sort(unique(c(de$k, obs_long$k))))], sort(unique(c(de$k, obs_long$k))))
p <- ggplot(de, aes(time, freq, colour = fit)) +
  facet_grid(paste("window", window) ~ lab[k], scales = "free") +
  geom_line(linewidth = 0.4) +
  geom_point(data = obs_long, aes(time, freq), inherit.aes = FALSE, size = 0.8, colour = "black") +
  scale_colour_manual(NULL, values = fit_colours) +
  labs(x = "time (days)", y = "karyotype frequency", subtitle = "median of bootstrap projections (minobs 20); points observed") +
  theme_classic(base_size = 7) + theme(strip.background = element_blank(), legend.position = "top")
ggsave(file.path(fig, "supp2de_example.png"), p, width = 200, height = 90, units = "mm", dpi = 200, bg = "white")
print(c_tab, row.names = FALSE)
print(de_tab, row.names = FALSE)
