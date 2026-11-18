# karst's choice between the balance prior and free initial abundances
# (fit_sparse_joint(initial_state = "auto"), the default from 827c4b5 on) on
# the two-time-point fits, and what the choice changed against a reference
# round with the balance prior on every fit.
#
#   Rscript scripts/initial_state_summary.R [REFERENCE_RESULTS]
#
# With two time points "auto" first fits the balance prior with its level
# free at the interval; a log level above the threshold (1) makes the fit
# "free", otherwise "balance". run.Rds keeps the check (R/fit.R:
# initial_state, initial_state_check). With more than two time points "auto"
# is "balance", so only the "ends" designs are read.
#
# Writes to BENCH_RESULTS:
#   tables/initial_state.tsv  one row per input and design: status, the state
#                             fitted, the check's level, sd, threshold,
#                             choice and convergence, and the grouping (ABM
#                             wavelength and window, the missegregation
#                             simulations' true rate, the Salehi context and
#                             passages)
#   analysis/initial_state/choice.tsv  share "free" and level quantiles by group
#   analysis/initial_state/effect.tsv  with REFERENCE_RESULTS: paired scores of
#                             the karst two-time-point fits of the two rounds
#                             (same input, design and minobs) by the state
#                             chosen; a fit kept on the balance prior should
#                             score the same
#   analysis/initial_state/figures/initial_state_level.pdf  the check's level by
#                             group, with the threshold

for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(ggplot2))
paths <- bench_paths()
args <- commandArgs(TRUE)
reference <- if (length(args)) normalizePath(args[[1]]) else NULL
out <- file.path(paths$results, "analysis", "initial_state")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(paths$results, "tables"), recursive = TRUE, showWarnings = FALSE)
options(width = 200)
write_tab <- function(d, file) utils::write.table(d, file, sep = "\t", quote = FALSE, row.names = FALSE)

## The choice -----------------------------------------------------------------------------
dirs <- Sys.glob(file.path(paths$results, "fits", "*", "*", "*", "karst"))
dirs <- dirs[grepl("ends$", basename(dirname(dirs)))]
rows <- lapply(dirs, function(d) {
  runs <- Sys.glob(file.path(d, "minobs_*", "run.Rds"))
  parts <- strsplit(d, "/", fixed = TRUE)[[1]]
  n <- length(parts)
  base <- data.frame(source = parts[n - 3], dataset = parts[n - 2], design = parts[n - 1], stringsAsFactors = FALSE)
  if (!length(runs)) return(cbind(base, status = "missing"))
  r <- readRDS(runs[[1]])
  ck <- r$initial_state_check
  cbind(base, status = if (is.null(r$fit_seconds)) "failed" else "ok",
        state = first_or_na(r$initial_state, NA_character_),
        level = if (is.null(ck)) NA_real_ else ck$level, level_sd = if (is.null(ck)) NA_real_ else ck$level_sd,
        threshold = if (is.null(ck)) NA_real_ else ck$threshold,
        chosen = if (is.null(ck)) NA_character_ else as.character(ck$chosen),
        check_convergence = if (is.null(ck)) NA_real_ else ck$check_convergence,
        history = first_or_na(r$history), convergence = first_or_na(r$convergence),
        fit_seconds = first_or_na(r$fit_seconds), stringsAsFactors = FALSE)
})
st <- bind_rows_fill(rows)
st$window <- NA_integer_
st$lambda <- NA_real_
st$rate <- NA_real_
sim <- st$source %in% c("abm", "misseg")
info <- do.call(rbind, lapply(st$dataset[sim], abm_sim_info))
st$lambda[sim] <- info$w
st$rate[sim] <- info$pmis
st$window[sim] <- as.integer(sub("^w([0-9]+)_ends$", "\\1", st$design[sim]))
meta <- utils::read.csv(paths$salehi_meta)
pdx_of <- stats::setNames(meta$PDX_id, paste(meta$datasetname, meta$timepoint, sep = "_"))
sal <- st$source == "salehi"
st$context <- NA_character_
st$context[sal] <- ifelse(pdx_of[sub("_l_.*$", "", st$dataset[sal])] %in% c("SA039", "SA906"), "in vitro", "PDX")
st$passages <- NA_integer_
st$passages[sal] <- as.integer(sub("^.*_l_([0-9]+)_.*$", "\\1", st$dataset[sal]))
st <- st[order(st$source, st$dataset, st$design), ]
write_tab(st, file.path(paths$results, "tables", "initial_state.tsv"))

fitted <- st[st$status == "ok" & !is.na(st$chosen), ]
q <- function(v) stats::quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE, names = FALSE)
choice_rows <- function(d, by, label) {
  do.call(rbind, lapply(split(d, d[by], drop = TRUE), function(x) {
    lv <- round(q(x$level), 2)
    cbind(data.frame(table = label, stringsAsFactors = FALSE), x[1, by, drop = FALSE],
          data.frame(inputs = nrow(x), free = sum(x$chosen == "free"), share_free = round(mean(x$chosen == "free"), 3),
                     level_q05 = lv[1], level_q25 = lv[2], level_median = lv[3], level_q75 = lv[4], level_q95 = lv[5],
                     check_unconverged = sum(x$check_convergence %in% 1)))
  }))
}
fitted$group <- ifelse(fitted$source == "salehi", fitted$context,
                       ifelse(fitted$source == "abm", paste("lambda", fitted$lambda), paste("rate", fitted$rate)))
fitted$passage_group <- ifelse(is.na(fitted$passages), NA_character_,
                               as.character(cut(fitted$passages, c(1, 2, 3, 5, Inf), labels = c("2", "3", "4-5", "6+"))))
choice <- bind_rows_fill(list(
  choice_rows(fitted, c("source", "design"), "source x design"),
  choice_rows(fitted[fitted$source == "abm", ], c("lambda", "window"), "ABM: wavelength x window"),
  choice_rows(fitted[fitted$source == "misseg", ], c("rate", "window"), "misseg: true rate x window"),
  choice_rows(fitted[fitted$source == "misseg", ], "rate", "misseg: true rate"),
  choice_rows(fitted[fitted$source == "salehi", ], "context", "Salehi: context"),
  choice_rows(fitted[fitted$source == "salehi", ], c("context", "passage_group"), "Salehi: context x passages")))
write_tab(choice, file.path(out, "choice.tsv"))
cat("two-time-point inputs:", nrow(st), " fitted with a check:", nrow(fitted), "  status:\n")
print(table(st$source, st$status))
for (t in unique(choice$table)) {
  cat("\n==", t, "\n")
  x <- choice[choice$table == t, ]
  print(x[, colSums(!is.na(x)) > 0 & names(x) != "table"], row.names = FALSE)
}

p <- ggplot(fitted, aes(x = group, y = level)) +
  geom_hline(yintercept = unique(stats::na.omit(fitted$threshold)), linetype = 2, colour = "grey40") +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  facet_grid(~ source, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "log level of the cloud over the balance\n(check fit; dashed: threshold, above it \"free\")") +
  theme_classic(base_size = 8) + theme(strip.background = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "figures", "initial_state_level.pdf"), p, width = 180, height = 70, units = "mm")

## What the choice changed (against REFERENCE_RESULTS) --------------------------------------
if (!is.null(reference)) {
  state_of <- stats::setNames(fitted$chosen, paste(fitted$source, fitted$dataset, fitted$design))
  read_tab <- function(root, name) {
    p <- file.path(root, "tables", paste0(name, ".Rds"))
    if (file.exists(p)) readRDS(p) else NULL
  }
  paired <- function(new, old, keys, value, label, higher = TRUE) {
    if (is.null(new) || is.null(old)) return(NULL)
    m <- merge(new[c(keys, value)], old[c(keys, value)], by = keys, suffixes = c(".new", ".old"))
    m <- m[is.finite(m[[paste0(value, ".new")]]) & is.finite(m[[paste0(value, ".old")]]), ]
    if (!nrow(m)) return(NULL)
    m$state <- state_of[paste(m$source, m[[keys[2]]], m$design)]
    m$state[is.na(m$state)] <- "no check"
    do.call(rbind, lapply(split(m, list(m$source, m$design, m$state), drop = TRUE), function(x) {
      a <- x[[paste0(value, ".new")]]
      b <- x[[paste0(value, ".old")]]
      d <- if (higher) a - b else b - a
      data.frame(score = label, source = x$source[1], design = x$design[1], state = x$state[1], pairs = nrow(x),
                 new = round(mean(a), 3), old = round(mean(b), 3), median_diff = round(stats::median(a - b), 4),
                 better = round(mean(d > 1e-12), 3), worse = round(mean(d < -1e-12), 3),
                 p = if (nrow(x) > 5 && any(abs(d) > 1e-12)) signif(suppressWarnings(stats::wilcox.test(a, b, paired = TRUE)$p.value), 2) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }
  two_tp <- function(d) d[d$method == "karst" & grepl("ends$", d$design), ]
  eff <- list()
  for (src in c("abm", "misseg")) {
    new <- read_tab(paths$results, paste0(src, "_dfl"))
    old <- read_tab(reference, paste0(src, "_dfl"))
    if (is.null(new) || is.null(old)) next
    new <- two_tp(new)
    old <- two_tp(old)
    for (v in intersect(c("r", "rho", "R"), names(new))) {
      eff[[length(eff) + 1L]] <- paired(new, old, c("source", "sim_id", "design", "minobs"), v, paste("landscape", v))
    }
  }
  new <- read_tab(paths$results, "abm_dfp")
  old <- read_tab(reference, "abm_dfp")
  if (!is.null(new) && !is.null(old)) {
    for (mt in c("euclidean", "wasserstein")) {
      a <- two_tp(new[new$passage == 10 & new$metric == mt, ])
      b <- two_tp(old[old$passage == 10 & old$metric == mt, ])
      a$win <- as.numeric(a$win)
      b$win <- as.numeric(b$win)
      eff[[length(eff) + 1L]] <- paired(a, b, c("source", "sim_id", "design", "minobs"), "win",
                                        paste0("forecast beats baseline, passage 10 (", mt, ")"))
    }
  }
  pf <- function(root) {
    p <- file.path(root, "analysis", "salehi", "per_fit.Rds")
    if (!file.exists(p)) return(NULL)
    d <- readRDS(p)
    d <- d[as.character(d$fit) == "KARST 2 time points", ]
    d$source <- "salehi"
    d$design <- "ends"
    d
  }
  a <- pf(paths$results)
  b <- pf(reference)
  if (!is.null(a) && !is.null(b)) {
    for (mt in c("wasserstein", "euclidean")) {
      eff[[length(eff) + 1L]] <- paired(a[a$metric == mt, ], b[b$metric == mt, ], c("source", "fi", "design", "min_obs"),
                                        "win", paste0("Salehi: share of daughters where the forecast beats the parent (",
                                                      mt, ")"))
    }
  }
  eff <- do.call(rbind, eff)
  if (!is.null(eff)) {
    eff$reference <- basename(reference)
    write_tab(eff, file.path(out, "effect.tsv"))
    cat("\n== Paired scores, this round against", basename(reference), "(karst, two time points), by the state chosen\n")
    print(eff[names(eff) != "reference"], row.names = FALSE)
  }
}
