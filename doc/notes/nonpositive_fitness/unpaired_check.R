# The benchmark's own simulator, one random stream shared by all clones, run
# with the same seed under either rule for clones with fitness <= 0 (one
# held-out landscape, patient 1). Shows why the comparison needs per-clone
# streams (abm_paired.cpp): the first clone removed shifts every later draw.
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/unpaired_check.R [LANDSCAPE]
# Notes: doc/notes/nonpositive_fitness.md.
here <- "doc/notes/nonpositive_fitness"
lid <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "heldout_L13"
for (p in c("landscape.R", "config.R", "generate.R")) source(file.path("benchmark/sparse_two_timepoint", p))
Rcpp::sourceCpp(file.path(here, "unpaired_abm.cpp"))
cfg <- stt_config(output_dir = file.path(tempdir(), "nonpositive_unpaired"))
land <- cfg$landscapes[cfg$landscapes$landscape_id == lid, ]
run <- function(rule) {
  died <- 0
  stt_abm_cpp <<- function(...) { r <- stt_abm_die_cpp(..., nonpositive_dies = rule == "die"); died <<- died + r$died_nonpositive; r }
  truth <- readRDS(stt_generate_truth(cfg, land, 1L, force = TRUE))
  list(truth = truth, died = died)
}
a <- run("persist"); b <- run("die")
fr <- function(x) x / sum(x)
tv <- function(x, y) { k <- union(names(x), names(y)); u <- setNames(numeric(length(k)), k); v <- u; u[names(x)] <- fr(x); v[names(y)] <- fr(y); sum(abs(u - v)) / 2 }
cat(lid, ": cells removed under the paper's rule:", b$died, "\n")
for (day in c("1440", "1620", "1800", "2520")) {
  x <- a$truth$snapshots[[day]]; y <- b$truth$snapshots[[day]]
  cat(sprintf("  day %s  total variation between the rules %.3f; largest clone %.3f / %.3f; same largest clone: %s\n",
              day, tv(x, y), max(fr(x)), max(fr(y)), names(which.max(x)) == names(which.max(y))))
}
