# Names of the fits compared, shared by the analysis scripts.
#
# Each method is fitted to two designs of an input: its first and last
# time point ("2 time points": ABM w2/w4/w8_ends, Salehi ends) and every time point
# ("all time points": ABM w4_all, w8_all, Salehi all). For a window or trajectory
# of two time points the two designs coincide, so its "2 time points" fit is also
# listed as "all time points". "paper (alfak main@4ad82b0)" is the paper's own fits (all
# time points), a reference only. karst's other ways of fitting more than two
# time points (fit_sparse_joint(series = ...), c7f809c on; methods
# karst_pairwise and karst_state_space, R/fit.R) are "all time points" fits of
# their own; "KARST all time points" is the default, "anchored". "KARST auto"
# is, per input, one of the three picked by a rule (scripts/auto_picks.R).
# "KARST beta auto" (ac36fd1 round) is karst with its missegregation rate
# picked per input (beta = "auto"), on both designs; "KARST beta 5e-5" is
# karst at 5e-5 on the missegregation simulations, whose "karst" fits take
# the true rate.

fit_levels <- c("KARST 2 time points", "KARST all time points", "KARST pairwise all time points",
                "KARST state space all time points", "KARST auto all time points",
                "KARST beta auto 2 time points", "KARST beta auto all time points",
                "KARST beta 5e-5 2 time points", "KARST beta 5e-5 all time points", "alfak dev2@f93b2d2 2 time points",
                "alfak dev2@f93b2d2 all time points", "paper (alfak main@4ad82b0)")
fit_colours <- c("KARST 2 time points" = "#D55E00", "KARST all time points" = "#E69F00",
                 "KARST pairwise all time points" = "#CC79A7", "KARST state space all time points" = "#7B3294",
                 "KARST auto all time points" = "#117733",
                 "KARST beta auto 2 time points" = "#1B7837", "KARST beta auto all time points" = "#7FBC41",
                 "KARST beta 5e-5 2 time points" = "#762A83", "KARST beta 5e-5 all time points" = "#C2A5CF",
                 "alfak dev2@f93b2d2 2 time points" = "#0072B2", "alfak dev2@f93b2d2 all time points" = "#56B4E9", "paper (alfak main@4ad82b0)" = "#999999")
method_label <- c(karst = "KARST", karst_pairwise = "KARST pairwise", karst_state_space = "KARST state space",
                  karst_auto = "KARST auto", karst_beta_auto = "KARST beta auto",
                  karst_beta_fixed = "KARST beta 5e-5", karst_free = "KARST free", alfak = "alfak dev2@f93b2d2")
fit_method <- function(fit) sub(" .*$", "", as.character(fit))
fit_group <- function(fit) ifelse(grepl("all time points$", fit), "all time points", "2 time points")

# `all` marks the rows of an "all time points" design, `two` the rows of a
# two-time-point input whose fit also stands for "all time points" (for karst's
# other methods too, which give the same fit with two time points, when the
# round has them).
label_fit_rows <- function(d, all, two) {
  lab <- unname(method_label[d$method])
  lab[is.na(lab)] <- d$method[is.na(lab)]
  d$fit <- ifelse(d$method == "alfak_paper", "paper (alfak main@4ad82b0)", paste(lab, ifelse(all, "all time points", "2 time points")))
  # karst's missegregation-rate arms have two-time-point fits of their own
  dup <- d[two & d$method %in% c("karst", "alfak", names(karst_beta_methods)), , drop = FALSE]
  if (nrow(dup)) {
    dup$fit <- paste(unname(method_label[dup$method]), "all time points")
    series <- lapply(intersect(karst_other_methods, unique(d$method)), function(m) {
      x <- dup[dup$method == "karst", , drop = FALSE]
      x$fit <- rep(paste(method_label[[m]], "all time points"), nrow(x))
      x$method <- rep(m, nrow(x))
      x
    })
    d <- do.call(rbind, c(list(d, dup), series))
  }
  d$fit <- factor(d$fit, fit_levels)
  d
}

# ABM: designs w2_ends, w4_ends, w8_ends, w4_all, w8_all.
label_abm_fits <- function(d) {
  d <- label_fit_rows(d, grepl("_all$", d$design), d$design == "w2_ends")
  d$window <- c(w2_ends = 2L, w4_ends = 4L, w8_ends = 8L, w4_all = 4L, w8_all = 8L)[d$design]
  d$minobs <- as.character(d$minobs)
  d
}

# rbind of data frames whose columns differ; missing columns become NA.
bind_rows_fill <- function(rows) {
  rows <- Filter(function(r) !is.null(r) && nrow(r), rows)
  if (!length(rows)) return(NULL)
  cols <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(r) {
    r[setdiff(cols, names(r))] <- NA
    r[cols]
  }))
}

# Salehi: designs ends and all; trajectories of 2 passages ("_l_2_") have
# "ends" only.
label_salehi_fits <- function(d) {
  d <- label_fit_rows(d, d$design == "all", d$design == "ends" & grepl("_l_2_", d$fi))
  d$min_obs <- as.character(d$min_obs)
  d
}
