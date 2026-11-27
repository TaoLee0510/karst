# Whether karst can fit an input, decided from the data before any fit.
#
# karst's model domain is its major karyotypes (share >= major_frequency in
# any time point and >= major_count cells in the time points together) and
# their one-step neighbours. A time point with no cell inside the domain leaves
# nothing to fit: karst stops with "A timepoint has no cell inside the model
# domain." (class karst_domain_coverage_error from 2d86b39 on). The domain is
# built here with the installed karst's own functions and the defaults of
# fit_sparse_joint(), on the two time points fit_karst() fits. Such an input is
# "not analysable" and is left out of every method's analysis; the share of
# each time point inside the domain is kept for a sensitivity analysis.

# The same holds for an input of more than two time points (multi-timepoint
# karst, 21aedca on): by default (`intermediate = "majors"`) the majors are
# taken from the first and the last time point only, the domain is built around
# them, and karst stops when any time point, those in between included, has no
# cell inside the domain. From the commit after 21aedca on, a time point in
# between with no counted cell (none on a major) is left out of the fit, and
# only the first and the last time point must have a cell inside the domain
# (karst_drops_empty_samples()). Each design of an input ("ends", "all") is
# judged on the time points it fits, and a fit is analysed only when its own input
# is analysable, for both methods.
#
# karst's pairwise and state-space methods (fit_sparse_joint(series = ...),
# c7f809c on) give every pair of consecutive time points its own majors and
# domain and stop when a time point of a pair has no cell inside its pair's
# domain (`pair_failures`: such time points among all pairs). From 05a8855 on they
# leave such time points out instead, as anchored does, until every pair is
# covered (karst:::sparse_joint_pair_left_out(); `pair_left_out`); with only
# the first and last left, the two-time-point check applies. `pairs_analysable`
# is whether the fit can go on. Neither is used for the exclusion: the inputs
# analysed stay those of the default ("anchored") method, the time points of the
# round of record, and a fit that stops is a failed fit of that method.

# The count matrix karst fits. Two time points: as every commit reads them.
# More time points: the multi-timepoint prepare_karst_data(), which takes the
# intervals as dt and, if needed, the sampling times from 0.
karst_counts <- function(x, times, pm) {
  prep <- karst:::prepare_karst_data
  first <- tryCatch(prep(x, dt = diff(times), beta = pm)$counts, error = function(e) e)
  if (!inherits(first, "error")) return(first)
  if (ncol(x) > 2L && "times" %in% names(formals(prep))) {
    return(prep(x, dt = diff(times), beta = pm, times = times - times[1])$counts)
  }
  stop(first)
}

# One row: number of time points, cells of the first and last time point and how
# many lie inside the domain, karyotypes shared by the first time point and the
# others, majors, coverage (the lowest share of a time point inside the domain),
# the fewest cells of a time point inside the domain, and analysable (every
# time point has a cell inside the domain). NA where the installed karst cannot
# read the input (more than two time points before the multi-timepoint version).
karst_input_coverage <- function(yi, pm, min_cn = NULL) {
  f <- formals(karst::fit_sparse_joint)
  if (is.null(min_cn)) min_cn <- eval(f$min_cn)
  x <- as.matrix(yi$x)
  storage.mode(x) <- "double"
  times <- input_times(yi)
  k <- ncol(x)
  counts <- tryCatch(karst_counts(x, times, pm), error = function(e) NULL)
  if (is.null(counts)) {
    return(data.frame(samples = k, cells_1 = sum(x[, 1]), cells_2 = sum(x[, k]), in_domain_1 = NA_real_,
                      in_domain_2 = NA_real_, shared = NA_real_, majors = NA_integer_, coverage = NA_real_,
                      min_in_domain = NA_real_, mid_no_domain = NA_integer_, mid_no_major = NA_integer_,
                      mid_neighbours_only = NA_integer_, mid_left_out = NA_integer_, analysable = NA,
                      pair_failures = NA_integer_, pair_left_out = NA_integer_, pairs_analysable = NA))
  }
  max_cn <- max(karst:::parse_karyotypes(rownames(counts))) + 1L
  # The time points the majors come from, as the fit chooses them with its default
  # `intermediate` (every time point in commits without that choice).
  args <- list(counts, eval(f$major_frequency), eval(f$major_count), eval(f$max_states),
               as.integer(min_cn), as.integer(max_cn))
  intermediate <- if ("intermediate" %in% names(f)) eval(f$intermediate)[[1L]] else "all"
  if ("major_columns" %in% names(formals(karst:::sparse_joint_domain))) {
    args$major_columns <- karst:::sparse_joint_major_columns(list(counts = counts),
                                                             list(intermediate = intermediate))
  }
  domain <- do.call(karst:::sparse_joint_domain, args)
  cov <- karst:::sparse_joint_domain_coverage(counts, domain, eval(f$major_frequency), eval(f$major_count))
  # Time points in between with no cell inside the domain, and with none on a
  # major (cells on the majors' neighbours only); with intermediate =
  # "majors" a time point in between counts its cells on the majors only.
  mid <- if (k > 2L) 2:(k - 1L) else integer()
  ends_in <- cov$cells_in_domain[c(1L, k)] > 0
  drops <- karst_drops_empty_samples()
  analysable <- if (drops) all(ends_in) else all(cov$cells_in_domain > 0)
  # the time points in between the fit leaves out (commits that drop them)
  left_out <- if (!drops) NA_integer_ else if (k <= 2L) 0L else
    length(karst:::sparse_joint_uncounted_samples(cov, intermediate))
  # the pair check (pairwise and state space; two time points: the check above)
  pair_args <- list(eval(f$major_frequency), eval(f$major_count), eval(f$max_states), as.integer(min_cn),
                    as.integer(max_cn))
  pair_cov <- function(y) do.call(karst:::sparse_joint_pair_coverage, c(list(y), pair_args))
  pair_failures <- pair_left_out <- NA_integer_
  pairs_analysable <- NA
  if (k <= 2L) {
    pair_failures <- sum(cov$cells_in_domain <= 0)
    pair_left_out <- 0L
    pairs_analysable <- pair_failures == 0L
  } else if (karst_has_series_methods()) {
    pair_failures <- sum(pair_cov(counts)$cells_in_domain <= 0)
    pairs_analysable <- pair_failures == 0L
    if (exists("sparse_joint_pair_left_out", envir = asNamespace("karst"), inherits = FALSE)) {
      out <- do.call(karst:::sparse_joint_pair_left_out, c(list(counts), pair_args))$left_out
      pair_left_out <- length(out)
      kept <- setdiff(seq_len(k), out)
      pairs_analysable <- all(pair_cov(counts[, kept, drop = FALSE])$cells_in_domain > 0)
    }
  }
  data.frame(samples = k, cells_1 = cov$cells[1], cells_2 = cov$cells[k], in_domain_1 = cov$cells_in_domain[1],
             in_domain_2 = cov$cells_in_domain[k], shared = cov$shared[1], majors = domain$n_major,
             coverage = min(cov$fraction), min_in_domain = min(cov$cells_in_domain),
             mid_no_domain = sum(cov$cells_in_domain[mid] <= 0),
             mid_no_major = sum(cov$majors[mid] == 0),
             mid_neighbours_only = sum(cov$majors[mid] == 0 & cov$cells_in_domain[mid] > 0),
             mid_left_out = left_out, analysable = analysable,
             pair_failures = pair_failures, pair_left_out = pair_left_out, pairs_analysable = pairs_analysable)
}

# Whether the installed karst has the pairwise and state-space methods.
karst_has_series_methods <- function() {
  exists("sparse_joint_pair_coverage", envir = asNamespace("karst"), inherits = FALSE)
}

# Whether the installed karst leaves out a time point in between that has no
# counted cell (sparse_joint_uncounted_samples(); recorded in
# diagnostics$dropped_samples) instead of stopping: then only the first and
# the last time point must have a cell inside the domain.
karst_drops_empty_samples <- function() {
  exists("sparse_joint_uncounted_samples", envir = asNamespace("karst"), inherits = FALSE)
}

# Inputs (one row per input and design, with `analysable`) that are not
# analysable, as keys of input_key(); NA (not judged) counts as analysable.
input_key <- function(...) paste(..., sep = "|")
not_analysable_keys <- function(inputs, cols) {
  bad <- inputs[inputs$analysable %in% FALSE, cols, drop = FALSE]
  do.call(input_key, unname(as.list(bad)))
}

# Which candidate history karst chose (af362af on): "dt" or one of the fixed
# times above dt. The candidates are dt followed by the fixed times above it,
# so with k candidates the last k - 1 fixed times are the others; dt itself
# may equal a fixed time (e.g. two 15-day passages, dt = 30). A fit with one
# history and no profile (candidates NA: initial_state "free", where the
# history is dt and does not enter the model; 827c4b5's "auto") is "dt".
karst_chosen_history <- function(history, candidates) {
  fixed <- c(30, 100, 300, 1000, 3000)
  vapply(seq_along(history), function(i) {
    if (is.na(candidates[i]) || candidates[i] < 2L) return("dt")
    others <- utils::tail(fixed, candidates[i] - 1L)
    if (history[i] %in% others) as.character(history[i]) else "dt"
  }, character(1))
}

# The same for a fit of more than two time points, whose candidates are not
# known here: one of the fixed times, or "other".
karst_history_class <- function(history) {
  ifelse(history %in% c(30, 100, 300, 1000, 3000), as.character(history), "other")
}

# Every Salehi trajectory, reduced to its first and last passage ("ends")
# and whole ("all"; for 2-passage trajectories the same two passages).
salehi_input_coverage <- function(paths = bench_paths(), settings = bench_settings$salehi) {
  traj <- salehi_trajectories(paths)
  do.call(rbind, lapply(seq_len(nrow(traj)), function(i) {
    do.call(rbind, lapply(c("ends", "all"), function(design) {
      yi <- salehi_input(traj$traj_id[i], ends = design == "ends", paths)
      cbind(fi = traj$traj_id[i], design = design, gap = traj$length[i] - 1L,
            karst_input_coverage(yi, settings$pm, settings$karst_args$min_cn), stringsAsFactors = FALSE)
    }))
  }))
}

# Every ABM or missegregation simulation and design (window, first and last
# or every time point).
abm_input_coverage <- function(source, paths = bench_paths()) {
  raw <- if (source == "abm") paths$abm_raw else paths$misseg_raw
  sims <- sub("\\.Rds$", "", list.files(raw, pattern = "\\.Rds$"))
  designs <- abm_designs
  do.call(rbind, lapply(sims, function(id) {
    sim <- read_abm_sim(id, raw)
    pm <- if (source == "abm") bench_settings$abm$pm else abm_sim_info(id)$pmis
    do.call(rbind, lapply(seq_len(nrow(designs)), function(j) {
      cbind(source = source, sim_id = id, design = designs$design[j],
            karst_input_coverage(abm_input(sim, designs[j, ]), pm), stringsAsFactors = FALSE)
    }))
  }))
}
