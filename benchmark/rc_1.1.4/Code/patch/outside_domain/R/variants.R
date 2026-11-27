# Candidate rules for karyotypes outside a sparse joint fit's domain
# (PLAN.md, step 2), as a wrapper around karst's predict(); the package is not
# changed. `k` is karst's namespace (od_karst()).
#   current  predict() as it is: an observed karyotype left out of the domain
#            gets the observed level (mean log birth of the observed non-major
#            karyotypes inside); any other karyotype the trend.
#   B1       observed karyotypes left out of the domain: relative fitness 0
#            (the count-weighted mean of the majors, the fit's reference).
#   B2       unseen karyotypes outside the domain that are closer to an
#            observed out-of-domain karyotype o than to any major:
#            log birth = log birth(o) + trend(k) - trend(o) (o as in `current`).
#   B1B2     both (o at the B1 level).
# `nearest`: for B2, a list(index = position in `labels` of the anchor(s),
# distance) per label, built with od_anchor_index(); NULL computes it.
# Returns data.frame(karyotype, in_domain, observed_out, fitness_mean,
# relative_fitness, division_rate) for the variant.
od_predict_variant <- function(k, fit, labels, variant = c("current", "B1", "B2", "B1B2"), base = NULL, anchors = NULL) {
  variant <- match.arg(variant)
  if (is.null(base)) base <- k$predict.karst_sparse_joint_fit(fit, labels)
  cn <- od_parse(labels)
  self <- 2 * k$sparse_joint_self_probability(cn, fit$data$beta) - 1
  out <- data.frame(karyotype = labels, in_domain = base$in_domain, observed_out = labels %in% fit$dropped,
                    fitness_mean = base$fitness_mean, relative_fitness = base$relative_fitness, stringsAsFactors = FALSE)
  if (variant != "current") {
    trend <- fit$trend
    outside <- which(!base$in_domain)
    design <- k$sparse_joint_design(cn[outside, , drop = FALSE], trend$centre, trend$major, trend$shell_power)$design
    tr <- rep(NA_real_, length(labels)); tr[outside] <- as.numeric(design %*% trend$coefficients)
    lb <- rep(NA_real_, length(labels)); lb[outside] <- log(base$fitness_mean[outside] / self[outside])
    obs <- out$observed_out
    ref <- trend$reference$fitness
    if (variant %in% c("B1", "B1B2")) lb[obs] <- log(ref / self[obs])
    if (variant %in% c("B2", "B1B2")) {
      if (is.null(anchors)) stop("B2 needs anchors (od_anchor_index())")
      for (i in which(!base$in_domain & !obs & !is.na(anchors$distance))) {
        o <- anchors$index[[i]]
        lb[i] <- mean(lb[o] - tr[o]) + tr[i]
      }
    }
    fm <- self * exp(lb)
    out$fitness_mean[outside] <- fm[outside]
    out$relative_fitness[outside] <- fm[outside] - ref
  }
  out$division_rate <- out$fitness_mean / self
  out
}

# For every label: the observed out-of-domain karyotypes nearest to it among
# those within two steps (through the list's one-step index `eng`,
# engine_build()), if they are closer than every major (distance to the
# majors `dmajor`). Labels with none get distance NA.
od_anchor_index <- function(eng, observed_out, dmajor) {
  n <- length(eng$labels)
  best <- rep(NA_real_, n); idx <- vector("list", n)
  step1 <- function(i) { v <- c(eng$plus[i, ], eng$minus[i, ]); v[!is.na(v)] }
  for (o in which(observed_out)) {
    one <- unique(step1(o))
    two <- setdiff(unique(unlist(lapply(one, step1))), c(one, o))
    for (d in 1:2) {
      targets <- if (d == 1) one else two
      targets <- targets[d < dmajor[targets]]
      for (t in targets) {
        if (is.na(best[t]) || d < best[t]) { best[t] <- d; idx[[t]] <- o }
        else if (d == best[t]) idx[[t]] <- c(idx[[t]], o)
      }
    }
  }
  list(index = idx, distance = best)
}
