# The analyses of the paper's Fig. 4-6 on the Salehi fits of every method
# (scripts/fig4_cna_effects.R, fig5_wgd.R, fig6_steady_states.R,
# fig456_summary.R), ported from the paper's scripts 06_cna_fitness_effects.R,
# 07_wgd.R and 08_misseg_rate_analysis.R. They read each fit's landscape.Rds
# (alfak's format; karst's from as_alfak()) and the fit's row of 03 (xv,
# ntrain, train_treat, test_treat), and select fits as the paper does.

# Karyotype helpers of the paper (R/utils_karyo.R), unchanged.
paper_s2v <- function(charvec) as.numeric(unlist(strsplit(charvec, split = "[.]")))
paper_gen_all_neighbours <- function(ids, as.strings = TRUE, remove_nullisomes = TRUE) {
  if (as.strings) ids <- lapply(ids, paper_s2v)
  nkern <- do.call(rbind, lapply(1:length(ids[[1]]), function(i) {
    x0 <- rep(0, length(ids[[1]]))
    x1 <- x0
    x0[i] <- -1
    x1[i] <- 1
    rbind(x0, x1)
  }))
  n <- do.call(rbind, lapply(ids, function(ii) t(apply(nkern, 1, function(i) i + ii))))
  n <- unique(n)
  nids <- length(ids)
  n <- rbind(do.call(rbind, ids), n)
  n <- unique(n)
  n <- n[-(1:nids), ]
  if (remove_nullisomes) n <- n[apply(n, 1, function(ni) sum(ni < 1) == 0), ]
  n
}

# The fits compared: each method's fits of the first and last passage ("2
# time points") and of every passage ("all time points"; a trajectory of two passages
# has only the former, which stands for both), and the paper's own fits
# (alfak main@4ad82b0). The fig4 array has 10 tasks per arm (tasks 1-70), the
# fig6 array one task per arm (tasks 1-7); the paper's arm is the seventh.
fig456_arms <- data.frame(
  key = c("karst_2", "karst_all", "karst_pairwise_all", "karst_state_space_all", "alfak_2", "alfak_all", "paper"),
  fit = c("KARST 2 time points", "KARST all time points", "KARST pairwise all time points",
          "KARST state space all time points", "alfak dev2@f93b2d2 2 time points", "alfak dev2@f93b2d2 all time points", "paper (alfak main@4ad82b0)"),
  method = c("karst", "karst", "karst_pairwise", "karst_state_space", "alfak", "alfak", "alfak_paper"),
  all = c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE)

# Folder of one fit of an arm.
fig456_fit_dir <- function(arm, fi, min_obs, paths = bench_paths()) {
  if (arm$method == "alfak_paper") return(file.path(paths$paper_salehi, fi, paste0("minobs_", min_obs)))
  design <- if (arm$all && !grepl("_l_2_", fi)) "all" else "ends"
  method <- if (design == "ends" && arm$method %in% karst_other_methods) "karst" else arm$method
  file.path(paths$results, "fits", "salehi", fi, design, method, paste0("minobs_", min_obs))
}

# The rows of 03 of every fit of an arm, as a list with one data frame per
# fit (as the paper's alfak_outputs_proc/ files): the paper's own files for
# the paper's fits, this round's scores (tables/salehi_scores.Rds) otherwise.
# Trajectories that are not analysable (R/coverage.R) are left out unless
# `exclude` is FALSE.
fig456_fit_rows <- function(arm, paths = bench_paths(), exclude = TRUE) {
  if (arm$method == "alfak_paper") {
    files <- list.files(file.path(paths$data, "processed", "salehi", "alfak_outputs_proc"), full.names = TRUE)
    rows <- lapply(files, readRDS)
  } else {
    s <- readRDS(file.path(paths$results, "tables", "salehi_scores.Rds"))
    two <- grepl("_l_2_", s$fi)
    want_design <- ifelse(arm$all & !two, "all", "ends")
    want_method <- ifelse(want_design == "ends" & arm$method %in% karst_other_methods, "karst", arm$method)
    s <- s[s$design == want_design & s$method == want_method, setdiff(names(s), c("design", "method"))]
    rows <- split(s, paste(s$fi, s$min_obs), drop = TRUE)
    names(rows) <- NULL
  }
  if (exclude) {
    inputs <- readRDS(file.path(paths$results, "tables", "salehi_inputs.Rds"))
    bad <- unique(inputs$fi[inputs$analysable %in% FALSE])
    rows <- Filter(function(r) !r$fi[1] %in% bad, rows)
  }
  rows
}

# Reading a landscape once per fit folder.
fig456_landscape <- local({
  cache <- new.env()
  function(dir) {
    if (!exists(dir, envir = cache, inherits = FALSE)) {
      path <- file.path(dir, "landscape.Rds")
      assign(dir, if (file.exists(path)) readRDS(path) else NULL, envir = cache)
    }
    get(dir, envir = cache, inherits = FALSE)
  }
})
