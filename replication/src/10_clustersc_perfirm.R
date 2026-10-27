# ClusterSC baseline on Geithner Panel B, per-firm version (as in the ClusterSC paper).
#
# Each (pseudo-)treated firm is its own target: it is assigned to a cluster and gets its own SC; the cohort
# effect is the average of the firms' effects. Same six estimators, clusters (r = 10) and placebo design as
# 09_clustersc_panelB.R (cohort-average version), so the two versions can be compared draw by draw only in
# distribution, not pairwise (different draws).
#
# Usage:
#   run:     Rscript --vanilla replication/src/10_clustersc_perfirm.R run <panelB_rds> <crsp_dir> <out_rds> \
#              [n_draws=50] [N1=12] [seed=1] [pools=all,financial] [draw_offset=0] [do_real=1]
#   combine: Rscript --vanilla replication/src/10_clustersc_perfirm.R combine <out_md_prefix> <rds>...
#   Run several "run" processes with different pools / draw_offset to parallelise, then "combine".

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))

if (mode == "combine") {
  prefix <- args[2]
  parts <- lapply(args[-(1:2)], readRDS)
  plac <- do.call(rbind, lapply(parts, `[[`, "placebo"))
  real <- do.call(rbind, lapply(parts, `[[`, "real"))
  summ <- do.call(rbind, lapply(split(plac, list(plac$pool, plac$variant), drop = TRUE), function(d) data.frame(
    pool = d$pool[1], variant = d$variant[1], draws = nrow(d), bias = 100 * mean(d$estimate),
    bias_se = 100 * sd(d$estimate) / sqrt(nrow(d)), rmse = 100 * sqrt(mean(d$estimate^2)), pre_rmse = 100 * mean(d$pre_rmse))))
  cmp <- list(c("simplex_cluster", "simplex_all"), c("simplex_cluster", "simplex_random"),
              c("robust_cluster", "robust_all"), c("robust_cluster", "robust_random"), c("robust_all", "simplex_all"))
  pair <- do.call(rbind, lapply(unique(plac$pool), function(pool) {
    w <- reshape(plac[plac$pool == pool, c("draw", "variant", "estimate")], idvar = "draw", timevar = "variant", direction = "wide")
    do.call(rbind, lapply(cmp, function(c2) {
      dd <- 1e4 * (w[[paste0("estimate.", c2[1])]]^2 - w[[paste0("estimate.", c2[2])]]^2)
      data.frame(pool = pool, comparison = paste(c2, collapse = " vs "), mean_diff = mean(dd),
                 se = sd(dd) / sqrt(length(dd)), t = mean(dd) / (sd(dd) / sqrt(length(dd))))
    }))
  }))
  cat("Placebo (true effect 0), percent per day:\n"); print(summ[order(summ$pool, summ$rmse), ], digits = 3, row.names = FALSE)
  cat("\nPaired differences in squared error (x 1e4); negative = first is better:\n"); print(pair, digits = 3, row.names = FALSE)
  cat("\nReal treated cohorts, percent per day:\n"); print(transform(real, estimate = round(100 * estimate, 3), pre_rmse = round(100 * pre_rmse, 3)), row.names = FALSE)
  saveRDS(list(placebo = plac, real = real, summary = summ, paired = pair), paste0(prefix, ".rds"))
  quit(save = "no")
}

stopifnot(mode == "run")
panel_rds <- args[2]; crsp_dir <- args[3]; out_rds <- args[4]
n_draws     <- if (length(args) >= 5) as.integer(args[5]) else 50L
N1          <- if (length(args) >= 6) as.integer(args[6]) else 12L
seed        <- if (length(args) >= 7) as.integer(args[7]) else 1L
pool_names  <- if (length(args) >= 8) strsplit(args[8], ",")[[1]] else c("all", "financial")
draw_offset <- if (length(args) >= 9) as.integer(args[9]) else 0L
do_real     <- if (length(args) >= 10) as.integer(args[10]) == 1 else TRUE

source(file.path(script_dir, "sc_qp.R"))
source(file.path(script_dir, "clustersc.R"))
t_start <- Sys.time()
log_step <- function(...) { cat(sprintf("[%6.0fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = ""); flush.console() }

R_RANK <- 10L
p <- readRDS(panel_rds)
pre <- p$dif < 0; post <- p$dif >= 0
Y <- p$Y
ctrl <- which(rowSums(p$treat) == 0)
nm <- readRDS(file.path(crsp_dir, "crsp_names.rds"))
nm$namedt <- as.Date(nm$namedt); nm$nameendt <- as.Date(nm$nameendt)
ev <- nm[nm$namedt <= as.Date("2008-11-21") & nm$nameendt >= as.Date("2008-11-21"), ]
sic <- ev$siccd[match(as.integer(p$permno), ev$permno)]
financial <- !is.na(sic) & sic >= 6000 & sic <= 6799

cl <- cluster_donors(Y[ctrl, pre], r = R_RANK, k_max = 10, seed = seed)   # same clusters as 09 (same seed)
log_step(sprintf("clusters r = %d: k = %d, sizes %s", R_RANK, cl$k, paste(cl$sizes, collapse = "/")))

# One cohort: each firm gets its own SC; returns the cohort-average effect per estimator.
estimate_perfirm <- function(firms, donors, rng_seed) {
  lab <- cl$labels[match(donors, ctrl)]
  M_all <- hsvt(Y[donors, ], R_RANK)
  M_cl <- list()
  rows <- list()
  set.seed(rng_seed)
  for (i in firms) {
    y <- Y[i, ]
    g <- unname(assign_cluster(y[pre], cl))
    in_cl <- donors[lab == g]
    key <- as.character(g)
    if (is.null(M_cl[[key]])) M_cl[[key]] <- hsvt(Y[in_cl, ], R_RANK)
    rnd <- sample(donors, length(in_cl))
    res <- list(
      simplex_all     = simplex_sc(Y[donors, ], y, pre, post),
      simplex_cluster = simplex_sc(Y[in_cl, ], y, pre, post),
      simplex_random  = simplex_sc(Y[rnd, ], y, pre, post),
      robust_all      = robust_sc(Y[donors, ], y, pre, post, R_RANK, M = M_all),
      robust_cluster  = robust_sc(Y[in_cl, ], y, pre, post, R_RANK, M = M_cl[[key]]),
      robust_random   = robust_sc(Y[rnd, ], y, pre, post, R_RANK))
    rows[[length(rows) + 1]] <- data.frame(firm = i, cluster = g, variant = names(res),
                                           estimate = sapply(res, `[[`, "estimate"), pre_rmse = sapply(res, `[[`, "pre_rmse"))
  }
  f <- do.call(rbind, rows)
  out <- aggregate(cbind(estimate, pre_rmse) ~ variant, f, mean)
  out$clusters <- paste(names(table(f$cluster[f$variant == "simplex_all"])), table(f$cluster[f$variant == "simplex_all"]), sep = "x", collapse = " ")
  out
}

real <- NULL
if (do_real) {
  real <- do.call(rbind, lapply(names(p$treat), function(cn) {
    tr <- which(p$treat[[cn]] == 1 & p$sample[[cn]])
    cbind(connection = cn, n_treated = length(tr), estimate_perfirm(tr, ctrl, seed))
  }))
  log_step("real treated cohorts done")
}

pools <- list(all = ctrl, financial = ctrl[financial[ctrl]])
plac <- list()
for (pool in pool_names) {
  for (b in draw_offset + seq_len(n_draws)) {
    set.seed(seed + match(pool, names(pools)) * 100000 + b)        # the draw depends only on (pool, b)
    pseudo <- sample(pools[[pool]], N1)
    d <- estimate_perfirm(pseudo, setdiff(ctrl, pseudo), seed + b)
    plac[[length(plac) + 1]] <- cbind(pool = pool, draw = b, d)
    if ((b - draw_offset) %% 10 == 0) log_step(sprintf("placebo %s: %d / %d", pool, b - draw_offset, n_draws))
  }
}
plac <- do.call(rbind, plac)
saveRDS(list(placebo = plac, real = real, clusters = cl, N1 = N1, seed = seed), out_rds)
log_step(paste("saved", out_rds))
