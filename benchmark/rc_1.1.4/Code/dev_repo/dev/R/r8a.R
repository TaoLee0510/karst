# R8(a): per-dataset parametric bootstrap through the dev ABM on the fitted
# landscape (spec: dev/phase3/r8a_spec.md). Sourced after common.R, eset.R,
# devbench.R and adapters.R, with karst loaded.

r8a_load_abm <- function() {
  if (exists("r8a_abm_cpp", mode = "function")) return(invisible())
  # one cache per R process: tasks compiling into a shared folder on one node collide
  cache <- file.path(tempdir(), "r8a_cpp")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(file.path(dev_repo, "dev", "abm", "abm_table.cpp"), cacheDir = cache, env = globalenv())
  invisible()
}

r8a_classes <- function(distance, count) {
  ifelse(distance == 0L, ifelse(!is.na(count) & count >= 5, "obs5+", "obs<5"),
         ifelse(distance == 1L, "unobs_d1", "unobs_d2"))
}
r8a_class_names <- c("obs<5", "obs5+", "unobs_d1", "unobs_d2")

# 95% quantile with every replicate's total weight within a class equal.
r8a_quantile <- function(s, replicate, p = 0.95) {
  if (!length(s)) return(NA_real_)
  w <- 1 / as.numeric(table(replicate)[as.character(replicate)])
  o <- order(s)
  cw <- cumsum(w[o]) / sum(w)
  s[o][which(cw >= p)[1L]]
}

r8a_q <- function(scores) {
  q <- vapply(r8a_class_names, function(cl) {
    x <- scores[scores$class == cl, ]
    r8a_quantile(x$s, x$replicate)
  }, numeric(1))
  fallback <- is.na(q)
  if (any(fallback)) q[fallback] <- max(q, na.rm = TRUE)
  list(q = q, fallback = names(q)[fallback])
}

r8a_one <- function(kind, id, round, n_rep = 50L, round1_path = NULL) {
  r8a_load_abm()
  dev_devbench_source()
  a <- dev_adapter(kind, id)
  d <- a$dataset
  abm <- stt_config(file.path(dev_results, "abm_cache", "r8a"))$abm
  fit <- suppressWarnings(do.call(karst::fit_sparse_joint, a$args))
  s <- fit$summary
  cn_dom <- karst:::parse_karyotypes(s$karyotype)
  storage.mode(cn_dom) <- "integer"
  b_dom <- exp(fit$posterior$g)
  coef <- as.numeric(fit$trend$coefficients)
  centre <- as.numeric(fit$trend$centre)
  majors <- fit$trend$major
  storage.mode(majors) <- "integer"
  shell_power <- fit$trend$shell_power
  start <- round(2e6 * s$expected_t0 / sum(s$expected_t0))
  names(start) <- s$karyotype
  start <- start[start > 0]
  start_cn <- karst:::parse_karyotypes(names(start))
  storage.mode(start_cn) <- "integer"

  # classes of the domain states in the original data, and round-2 spreads
  observed <- d$counts[rowSums(d$counts) > 0, , drop = FALSE]
  obs_cn <- karst:::parse_karyotypes(rownames(observed))
  dist_obs <- vapply(seq_len(nrow(cn_dom)), function(i) min(rowSums(abs(sweep(obs_cn, 2L, cn_dom[i, ])))), numeric(1))
  dom_class <- r8a_classes(dist_obs, s$count_t0 + s$count_t1)
  noise_dom <- rep(0, length(b_dom)); noise_out <- 0
  q1 <- NULL
  if (round == 2L) {
    q1 <- readRDS(round1_path)$q
    p0 <- stats::predict(fit, a$nodes$karyotype)
    sd_d2 <- stats::median(pmax(p0$relative_sd_model, 0.009 * (1 + p0$shell))[a$nodes$distance == 2L], na.rm = TRUE)
    noise_dom <- q1[dom_class] * pmax(s$relative_sd_model, 0.009 * (1 + s$shell))
    noise_out <- unname(q1[["unobs_d2"]]) * sd_d2
  }
  steps <- as.integer(round(d$interval / abm$delta_t))
  rows <- list(); failures <- character()
  for (r in seq_len(n_rep)) {
    seed <- (as.numeric(d$seed) %% 1e6) * 100 + r + 1000 * round
    set.seed(seed + 7)
    b_tab <- b_dom + if (round == 2L) stats::rnorm(length(b_dom), 0, noise_dom) else 0
    y0 <- stats::rmultinom(1L, d$depth_t0, start / sum(start))[, 1L]
    grow <- start
    if (isTRUE(d$bottleneck > 0)) {
      sv <- stats::rmultinom(1L, d$bottleneck, start / sum(start))[, 1L]
      grow <- sv
    }
    grow_cn <- start_cn[grow > 0, , drop = FALSE]
    raw <- r8a_abm_cpp(grow_cn, as.numeric(grow[grow > 0]), cn_dom, b_tab, coef, centre, majors, shell_power,
                       noise_out, seed, d$beta, abm$delta_t, steps, steps, abm$max_pop, abm$culling_survival,
                       abm$min_cn, abm$max_cn, seed)
    pop1 <- stats::setNames(raw$count, raw$karyotype)
    if (!length(pop1)) { failures <- c(failures, "extinct"); next }
    y1 <- stats::rmultinom(1L, d$depth_t1, pop1 / sum(pop1))[, 1L]
    labels <- union(names(y0)[y0 > 0], names(y1)[y1 > 0])
    counts <- cbind(t0 = unname(y0[labels]), t1 = unname(y1[labels]))
    counts[is.na(counts)] <- 0
    rownames(counts) <- labels
    rf <- tryCatch(suppressWarnings(karst::fit_sparse_joint(counts, dt = d$interval, beta = d$beta, birth_rate = 1,
                                                            min_cn = 1L, max_cn = 8L, initial_state = "balance")),
                   error = function(e) e)
    if (inherits(rf, "error")) { failures <- c(failures, conditionMessage(rf)); next }
    nodes <- stt_scoring_nodes(counts, 1L, 8L, 2000L, 2000L, as.integer(seed %% .Machine$integer.max))
    ncn <- karst:::parse_karyotypes(nodes$karyotype)
    storage.mode(ncn) <- "integer"
    truth <- r8a_birth_cpp(ncn, cn_dom, b_tab, coef, centre, majors, shell_power, noise_out, seed)
    p <- stats::predict(rf, nodes$karyotype)
    pred <- data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
                       model_sd = p$relative_sd_model, shell = p$shell, stringsAsFactors = FALSE)
    sc <- stt_score_predictions(nodes, truth, pred, 5)$nodes
    sc <- sc[sc$distance <= 2L, ]
    model_sd <- pred$model_sd[match(sc$karyotype, pred$karyotype)]
    shell <- pred$shell[match(sc$karyotype, pred$karyotype)]
    scale <- pmax(model_sd, 0.009 * (1 + shell))   # spec amendment 03:44: the package's floor term
    ok <- is.finite(scale) & scale > 0 & is.finite(sc$pred_relative)
    if (!any(ok)) { failures <- c(failures, "no scored node"); next }
    rows[[length(rows) + 1L]] <- data.frame(
      replicate = r, class = r8a_classes(sc$distance[ok], sc$count_t0[ok] + sc$count_t1[ok]),
      s = abs(sc$pred_relative[ok] - sc$truth_relative[ok]) / scale[ok], stringsAsFactors = FALSE)
  }
  scores <- do.call(rbind, rows)
  if (is.null(scores)) return(list(kind = kind, dataset_id = id, round = round, status = "failed",
                                   failures = failures))
  qq <- r8a_q(scores)
  list(kind = kind, dataset_id = id, round = round, status = "ok", q = qq$q, fallback = qq$fallback, q_round1 = q1,
       n_ok = length(unique(scores$replicate)), failures = failures,
       n_by_class = table(factor(scores$class, r8a_class_names)), scores = scores,
       noise_out = noise_out)
}
