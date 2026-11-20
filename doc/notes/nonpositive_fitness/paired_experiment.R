# Paired comparison of the two rules for f <= 0 on the own benchmark's held-out
# landscapes: same seeds, per-clone random streams (abm_paired.cpp), so that
# the only difference between the runs is the rule and what follows from it.
# Run from the repository root:  Rscript doc/notes/nonpositive_fitness/paired_experiment.R
# Notes: doc/notes/nonpositive_fitness.md. OWN_RESULTS is the own benchmark's held-out
# results (datasets/, nodes/, fits/, dataset_index_heldout.csv); the default is the
# copy made by doc/report/build/fetch_own.sh.
here <- "doc/notes/nonpositive_fitness"
own <- Sys.getenv("OWN_RESULTS", "doc/report/data/own_e7e32f8")
work <- file.path(tempdir(), "nonpositive_paired")   # truths and datasets, not kept
for (p in c("landscape.R", "config.R", "generate.R")) source(file.path("benchmark/sparse_two_timepoint", p))
Rcpp::sourceCpp(file.path(here, "abm_paired.cpp"))
RULE <- "persist"; DIED <- 0
stt_abm_cpp <- function(...) { r <- stt_abm_paired_cpp(..., nonpositive_dies = RULE == "die"); DIED <<- DIED + r$died_nonpositive; r }
base <- stt_config()
lands <- base$landscapes[base$landscapes$split == "heldout", ]
scen <- base$scenarios[base$scenarios$role == "core" | grepl("^turned", base$scenarios$scenario), ]
fr <- function(x) x / sum(x)
tv <- function(a, b) { k <- union(names(a), names(b)); x <- setNames(numeric(length(k)), k); y <- x; x[names(a)] <- fr(a); y[names(b)] <- fr(b); sum(abs(x - y)) / 2 }
job <- expand.grid(l = seq_len(nrow(lands)), patient = 1:8)
one <- function(j) {
  land <- lands[job$l[j], ]; patient <- job$patient[j]
  res <- list()
  for (rule in c("persist", "die")) {
    RULE <<- rule; DIED <<- 0
    config <- stt_config(output_dir = file.path(work, rule))
    tp <- stt_generate_truth(config, land, patient, force = TRUE); truth <- readRDS(tp)
    ds <- lapply(seq_len(nrow(scen)), function(s) {
      seed <- as.integer(land$landscape_seed %% 100000L * 1000L + patient * 20L + match(scen$scenario[s], base$scenarios$scenario))
      readRDS(stt_make_dataset(config, truth, scen[s, ], config$t0, seed, force = TRUE))
    })
    res[[rule]] <- list(truth = truth, ds = ds, died = DIED)
  }
  p <- res$persist; d <- res$die
  pop <- do.call(rbind, lapply(c("1440", "1620", "1800", "2520"), function(day) {
    a <- p$truth$snapshots[[day]]; b <- d$truth$snapshots[[day]]
    fa <- stt_grf_fitness(names(a), p$truth$centroids, p$truth$wavelength)
    data.frame(day = as.numeric(day), tv_rules = tv(a, b), nonpos_persist = sum(fr(a)[fa <= 0]))
  }))
  dsr <- do.call(rbind, lapply(seq_along(p$ds), function(s) {
    a <- p$ds[[s]]; b <- d$ds[[s]]
    same <- identical(a$counts, b$counts)
    ta <- a$truth$frequency_t1; tb <- b$truth$frequency_t1
    data.frame(scenario = a$scenario, counts_identical = same,
               cells_differing = if (same) 0 else { k <- union(rownames(a$counts), rownames(b$counts)); x <- matrix(0, length(k), 2, dimnames = list(k, NULL)); y <- x; x[rownames(a$counts), ] <- a$counts; y[rownames(b$counts), ] <- b$counts; sum(abs(x - y)) / 2 },
               tv_truth_t1 = tv(ta, tb))
  }))
  list(id = data.frame(landscape_id = land$landscape_id, wavelength = land$wavelength, patient = patient, died_total = d$died),
       pop = pop, ds = dsr)
}
t <- system.time(out <- parallel::mclapply(seq_len(nrow(job)), function(j) tryCatch(one(j), error = function(e) conditionMessage(e)), mc.cores = 3))
bad <- vapply(out, is.character, TRUE); if (any(bad)) { print(unique(unlist(out[bad]))); }
saveRDS(out, file.path(here, "paired_experiment.rds"))
cat("runs", sum(!bad), "of", length(out), "seconds", round(t[["elapsed"]]), "\n")
