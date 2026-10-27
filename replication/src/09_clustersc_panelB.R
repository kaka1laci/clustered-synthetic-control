# ClusterSC baseline on Geithner Panel B (all CRSP common stocks as donors), cohort-average version.
#
# Question: does restricting donors to the treated cohort's cluster make the counterfactual more accurate?
# Target = the (pseudo-)treated cohort's average daily return, as in GP&L's SC (Definition 5).
# Six estimators = {simplex SC (GP&L), robust SC (HSVT + ridge, ClusterSC paper)} x
#                  {all donors, the assigned cluster, a random donor subset of the cluster's size}.
# Clusters: k-means on SVD features of the donors' pre-period returns, k by silhouette. Daily returns are
# far from low rank (95% of the singular-value sum needs ~208 of 226 components), so the rank r is a
# parameter: main r = 10; ClusterSC sensitivity r = 5 and the paper's 95% rule.
#
# Evaluation: (1) the real treated cohorts (schedule / personal / New York); (2) placebo cohorts of N1 firms
# drawn from untreated firms (all, or financial firms SIC 6000-6799), true effect 0 -> bias and RMSE.
# Clusters are computed once on all untreated firms; for a placebo draw its own firms are removed from the donors.
#
# Usage: Rscript --vanilla replication/src/09_clustersc_panelB.R <panelB_rds> <crsp_dir> <out_rds> [n_draws=100] [N1=12] [seed=1]

args <- commandArgs(trailingOnly = TRUE)
panel_rds <- args[1]; crsp_dir <- args[2]; out_rds <- args[3]
n_draws <- if (length(args) >= 4) as.integer(args[4]) else 100L
N1      <- if (length(args) >= 5) as.integer(args[5]) else 12L
seed    <- if (length(args) >= 6) as.integer(args[6]) else 1L

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(script_dir, "sc_qp.R"))
source(file.path(script_dir, "clustersc.R"))

t_start <- Sys.time()
log_step <- function(...) {
  cat(sprintf("[%6.0fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = ""); flush.console()
}

p <- readRDS(panel_rds)
pre <- p$dif < 0; post <- p$dif >= 0
Y <- p$Y
treated_any <- rowSums(p$treat) > 0
ctrl <- which(!treated_any)

nm <- readRDS(file.path(crsp_dir, "crsp_names.rds"))
nm$namedt <- as.Date(nm$namedt); nm$nameendt <- as.Date(nm$nameendt)
ev <- nm[nm$namedt <= as.Date("2008-11-21") & nm$nameendt >= as.Date("2008-11-21"), ]
sic <- ev$siccd[match(as.integer(p$permno), ev$permno)]
financial <- !is.na(sic) & sic >= 6000 & sic <= 6799
log_step(sprintf("Panel B: %d firms, %d untreated donors (%d financial)", nrow(Y), length(ctrl), sum(financial[ctrl])))

# ---- clusters on all untreated firms ------------------------------------------------
sv <- svd(Y[ctrl, pre], nu = 0, nv = 0)$d
r95 <- which(cumsum(sv) / sum(sv) >= 0.95)[1]
ranks <- c(r10 = 10L, r5 = 5L, r95 = r95)
clusters <- lapply(ranks, function(r) cluster_donors(Y[ctrl, pre], r = r, k_max = 10, seed = seed))
for (nmr in names(clusters)) {
  cl <- clusters[[nmr]]
  fin_share <- tapply(financial[ctrl], cl$labels, mean)
  log_step(sprintf("clusters %s (r = %d): k = %d | sizes %s | financial share %s", nmr, ranks[nmr], cl$k,
                   paste(cl$sizes, collapse = "/"), paste(round(fin_share, 2), collapse = "/")))
}

# ---- one target: all estimators ------------------------------------------------------
estimate_all <- function(y1, donors, rng_seed) {
  out <- list()
  for (nmr in names(clusters)) {
    cl <- clusters[[nmr]]; r <- ranks[[nmr]]
    g <- assign_cluster(y1[pre], cl)
    in_cl <- donors[cl$labels[match(donors, ctrl)] == g]
    out[[paste0("robust_cluster_", nmr)]] <- c(robust_sc(Y[in_cl, ], y1, pre, post, r)[c("estimate", "pre_rmse")], n = length(in_cl), cl = g)
    if (nmr == "r10") {
      set.seed(rng_seed)
      rnd <- sample(donors, length(in_cl))
      out$simplex_all     <- c(simplex_sc(Y[donors, ], y1, pre, post)[c("estimate", "pre_rmse")], n = length(donors), cl = NA)
      out$simplex_cluster <- c(simplex_sc(Y[in_cl, ], y1, pre, post)[c("estimate", "pre_rmse")], n = length(in_cl), cl = g)
      out$simplex_random  <- c(simplex_sc(Y[rnd, ], y1, pre, post)[c("estimate", "pre_rmse")], n = length(rnd), cl = NA)
      out$robust_all      <- c(robust_sc(Y[donors, ], y1, pre, post, r)[c("estimate", "pre_rmse")], n = length(donors), cl = NA)
      out$robust_random   <- c(robust_sc(Y[rnd, ], y1, pre, post, r)[c("estimate", "pre_rmse")], n = length(rnd), cl = NA)
    }
  }
  do.call(rbind, lapply(names(out), function(v) data.frame(variant = v, estimate = out[[v]]$estimate,
                                                           pre_rmse = out[[v]]$pre_rmse, n_donors = out[[v]]$n, cluster = out[[v]]$cl)))
}

# ---- real treated cohorts -----------------------------------------------------------
real <- do.call(rbind, lapply(names(p$treat), function(cn) {
  tr <- which(p$treat[[cn]] == 1 & p$sample[[cn]])
  d <- estimate_all(colMeans(Y[tr, , drop = FALSE]), ctrl, seed)
  cbind(connection = cn, n_treated = length(tr), d)
}))
log_step("real treated cohorts done")

# ---- placebo cohorts ----------------------------------------------------------------
pools <- list(all = ctrl, financial = ctrl[financial[ctrl]])
plac <- list()
for (pool in names(pools)) {
  set.seed(seed + match(pool, names(pools)) * 1000)
  for (b in seq_len(n_draws)) {
    pseudo <- sample(pools[[pool]], N1)
    d <- estimate_all(colMeans(Y[pseudo, , drop = FALSE]), setdiff(ctrl, pseudo), seed + b)
    plac[[length(plac) + 1]] <- cbind(pool = pool, draw = b, d)
    if (b %% 20 == 0) log_step(sprintf("placebo %s: %d / %d", pool, b, n_draws))
  }
}
plac <- do.call(rbind, plac)

summ <- do.call(rbind, lapply(split(plac, list(plac$pool, plac$variant), drop = TRUE), function(d) data.frame(
  pool = d$pool[1], variant = d$variant[1], bias = mean(d$estimate), sd = sd(d$estimate),
  rmse = sqrt(mean(d$estimate^2)), pre_rmse = mean(d$pre_rmse), n_donors = median(d$n_donors))))
summ <- summ[order(summ$pool, summ$rmse), ]
cat("\nPlacebo (true effect 0), average daily effect over days 0..10, in percent:\n")
print(transform(summ, bias = round(100 * bias, 3), sd = round(100 * sd, 3), rmse = round(100 * rmse, 3),
                pre_rmse = round(100 * pre_rmse, 3)), row.names = FALSE)
cat("\nReal treated cohorts (percent):\n")
print(transform(real, estimate = round(100 * estimate, 3), pre_rmse = round(100 * pre_rmse, 3)), row.names = FALSE)

saveRDS(list(real = real, placebo = plac, summary = summ, clusters = clusters, ranks = ranks,
             financial = financial, ctrl = ctrl, N1 = N1, n_draws = n_draws, seed = seed), out_rds)
log_step(paste("saved", out_rds))
