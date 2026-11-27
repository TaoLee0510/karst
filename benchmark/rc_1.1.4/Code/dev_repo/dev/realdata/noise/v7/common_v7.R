# Shared settings of dev/realdata/noise/PLAN_v7.md: the v7 results root and its
# path guard, the 20 landscapes (seed base 131000, shared by the e sets and
# C2), the missegregation multiples, the depth pairs and the halves.
v7_root <- file.path(dev_results, "realdata", "noise", "v7")

# Every path a v7 script reads or writes among results must lie in v7_root
# (raw inputs of the paper, Salehi and NUMBAT excepted, read only).
v7_path <- function(...) {
  p <- file.path(v7_root, ...)
  if (any(grepl("(^|/)[.][.](/|$)", p))) stop("v7 path with '..': ", p[1], call. = FALSE)
  if (!all(startsWith(normalizePath(dirname(p), mustWork = FALSE), normalizePath(v7_root, mustWork = FALSE))))
    stop("v7 path outside ", v7_root, ": ", p, call. = FALSE)
  p
}
v7_check_read <- function(path) {
  bad <- grepl("/realdata/noise/v[1-6]/|/realdata/CE/|/abm_cache/|/sparse_two_timepoint", path)
  if (bad) stop("v7 must not read earlier simulations: ", path, call. = FALSE)
  # rc freeze: paths stored in the v7 datasets name the as-run results root;
  # read them from V7_UPSTREAM (dev/hpc/env.sh) instead
  as_run <- "/share/lab_crd/taoli/Project/karst_dev/results/"
  if (startsWith(path, as_run)) path <- file.path(Sys.getenv("V7_UPSTREAM", dev_results), substring(path, nchar(as_run) + 1L))
  path
}

v7_wavelengths <- c(0.2, 0.4, 0.8, 1.2, 1.6)
v7_landscapes <- local({
  w <- rep(v7_wavelengths, each = 4L)
  data.frame(landscape_id = sprintf("v7_L%02d", seq_along(w)), split = "v7",
             wavelength_class = sprintf("w%g", w), wavelength = w,
             landscape_seed = 131000L + seq_along(w),
             half = ifelse(rep(1:4, times = length(v7_wavelengths)) <= 2L, "dev", "confirm"),
             stringsAsFactors = FALSE)
})
v7_multiples <- c(1, 4, 16)
v7_beta_default <- 5e-5
v7_depth_pairs <- data.frame(d0 = c(50, 200, 500, 100, 3000), d1 = c(50, 200, 700, 3000, 100))

# Seeds from text (dataset names), in [1, 2^31 - 2].
v7_seed <- function(text) {
  v <- as.numeric(utf8ToInt(text))
  as.integer((sum(v * seq_along(v) * 7919) + sum(v^2) * 104729) %% 2147483646) + 1L
}

# Version stamp written into every truth file.
v7_stamp <- function(files) {
  list(karst = as.character(utils::packageVersion("karst")), karst_path = find.package("karst"),
       repo_commit = "942182c (dev; simulator code as 1.1.3 6ae2a24 except comments and the held-out list)",
       md5 = tools::md5sum(files), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}
