# Run from the karst repository root. BENCH_RESULTS may point to a read-only
# mirror of the selected fits, retaining their fits/salehi directory layout.
round_dir <- "benchmark/paper/results_6ae2a24"
fits <- file.path(Sys.getenv("BENCH_RESULTS", round_dir), "fits/salehi")
data_root <- Sys.getenv("PAPER_DATA", "../ALFA-K/data")
panel <- file.path(round_dir, "analysis/fig2c_two_timepoints/shared_full_parent")
out <- "benchmark/paper/docs/fig2c_ploidy_census_6ae2a24.tsv"
pairs <- read.delim(file.path(panel, "daughter_angles.tsv"))
windows <- read.delim(file.path(panel, "windows.tsv"))
lineages <- readRDS(file.path(data_root, "processed/salehi/lineages.Rds"))
meta <- read.csv(file.path(data_root, "raw/salehi/metadata.csv"))
stopifnot(nrow(pairs) == 44L, !anyDuplicated(pairs$daughter_uid), !anyDuplicated(meta$uid),
          nrow(windows) == 35L, !anyDuplicated(windows$window),
          length(unique(pairs$parent_uid)) == 35L, setequal(windows$window, pairs$window))

# Newcombe's difference of two Wilson score intervals (end minus start), in pp.
# These describe independent cell sampling, conditional on this fit's domain.
change <- function(numerator, denominator) {
  if (any(denominator == 0)) {
    return(data.frame(start=if (denominator[1]) numerator[1]/denominator[1] else NA_real_,
                      end=if (denominator[2]) numerator[2]/denominator[2] else NA_real_,
                      change_pp=NA_real_, direction="undefined", lower95_pp=NA_real_, upper95_pp=NA_real_))
  }
  p <- numerator / denominator
  z <- qnorm(0.975)
  scale <- 1 + z^2 / denominator
  center <- (p + z^2/(2*denominator)) / scale
  half <- z * sqrt(p*(1-p)/denominator + z^2/(4*denominator^2)) / scale
  lower <- pmax(0, center-half)
  upper <- pmin(1, center+half)
  delta <- p[2] - p[1]
  direction <- sign(as.double(numerator[2])*denominator[1] - as.double(numerator[1])*denominator[2])
  data.frame(start=p[1], end=p[2], change_pp=100*delta,
             direction=c("down","flat","up")[direction+2],
             lower95_pp=100*(delta - sqrt((p[2]-lower[2])^2 + (upper[1]-p[1])^2)),
             upper95_pp=100*(delta + sqrt((upper[2]-p[2])^2 + (p[1]-lower[1])^2)))
}

one_window <- function(i) {
  w <- windows[i,]
  daughters <- pairs[pairs$window == w$window,]
  m <- unique(daughters[c("parent_uid","lineage","training_start","training_end","n_parent","fit_verdict")])
  ids <- lineages[[w$window]]$ids
  endpoint_ids <- ids[c(1L,length(ids))]
  ix <- match(endpoint_ids,meta$uid)
  stopifnot(nrow(m) == 1L, !anyNA(ix), endpoint_ids[2] == m$parent_uid,
            all(daughters$daughter_uid %in% lineages[[w$window]]$dec1),
            all(daughters$minobs_alfak == w$minobs_alfak))
  family <- meta$PDX_id[ix[2]]
  stopifnot(family == (if (m$lineage == "p53 k.o") "SA906" else m$lineage))
  group <- if (family %in% c("SA609", "SA906")) family else "Other"
  raw <- readRDS(file.path(data_root,"processed/salehi/alfak_inputs",paste0(w$window,".Rds")))
  x <- as.matrix(raw$x[,c(1L,ncol(raw$x)),drop=FALSE])
  stopifnot(identical(paste0("X",colnames(x)),meta$timepoint[ix]),
            all(grepl(paste0(" X",colnames(x)[1],"$"),m$training_start)),
            all(grepl(paste0(" X",colnames(x)[2],"$"),m$training_end)))
  x <- x[rowSums(x)>0,,drop=FALSE]
  tag <- paste0("minobs_",w$minobs_alfak)
  fit <- readRDS(file.path(fits,w$window,"ends/karst",tag,"karst_fit.Rds"))
  labs <- rownames(x)
  stopifnot(fit$data$beta == 5e-5, !anyDuplicated(labs), !anyDuplicated(fit$summary$karyotype),
            setequal(labs,rownames(fit$data$counts)), all(x == fit$data$counts[labs,,drop=FALSE]))
  inside <- labs %in% fit$summary$karyotype
  stopifnot(setequal(labs[!inside],fit$dropped))
  parts <- strsplit(labs,".",fixed=TRUE)
  stopifnot(all(lengths(parts)==22L), all(grepl("^[0-9]+$",unlist(parts))))
  high <- vapply(parts,function(k) sum(as.integer(k)) >= 66L,logical(1))
  masks <- list(in_low=inside & !high, in_high=inside & high, out_low=!inside & !high, out_high=!inside & high)
  counts <- vapply(masks,function(rows) colSums(x[rows,,drop=FALSE]),numeric(2))
  totals <- rowSums(counts)
  outside <- counts[,"out_low"] + counts[,"out_high"]
  stopifnot(all(totals == colSums(x)), totals[2] == m$n_parent,
            all(totals == fit$diagnostics$domain_coverage$cells),
            all(totals-outside == fit$diagnostics$domain_coverage$cells_in_domain))
  absolute <- change(outside,totals)
  conditional <- change(counts[,"out_high"],outside)
  data.frame(group=group, family=family, window=w$window, start_uid=endpoint_ids[1],
             parent_uid=m$parent_uid, training_start=m$training_start, training_end=m$training_end,
             interval_days=diff(as.numeric(colnames(x)))*raw$dt,
             daughter_count=nrow(daughters), daughter_uids=paste(sort(daughters$daughter_uid),collapse=","),
             fit_verdict=m$fit_verdict, start_total=totals[1], end_total=totals[2],
             as.list(setNames(counts[1,],paste0("start_",colnames(counts)))),
             as.list(setNames(counts[2,],paste0("end_",colnames(counts)))),
             setNames(absolute,paste0("outside_share_",names(absolute))),
             setNames(conditional,paste0("high_within_outside_",names(conditional))),row.names=NULL)
}

result <- do.call(rbind,lapply(seq_len(nrow(windows)),one_window))
passage <- function(label) as.integer(sub("^.* X","",label))
result <- result[order(match(result$group,c("SA906","SA609","Other")),result$family,
                       sub(" X[0-9]+$","",result$training_end),passage(result$training_end)),]
stopifnot(sum(result$daughter_count)==44L, !anyDuplicated(result$parent_uid))
write.table(result,out,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
