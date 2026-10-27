# Check clustersc.R on the ClusterSC paper's own simulation design (their Section 6.1):
# two groups of sinusoidal signals, n_A = n_B = 500, T = 10, T0 = 8, Gaussian noise s.
# The paper reports lower post-period MSE for ClusterSC than for SC on all donors; our implementation
# should reproduce that direction.
# Usage: Rscript --vanilla replication/tests/test_clustersc.R <src_dir>

args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "sc_qp.R"))
source(file.path(args[1], "clustersc.R"))

set.seed(2025)
make_group <- function(n, a_shape, w_range, r = 3, TT = 10) {
  tt <- seq_len(TT) / TT
  V <- t(sapply(seq_len(r), function(j)
    rbeta(1, a_shape[1], a_shape[2]) * sin(2 * pi * runif(1, w_range[1], w_range[2]) * tt + rnorm(1))))
  matrix(runif(n * r), n, r) %*% V                       # each unit: random mix of the group's r signals
}
TT <- 10; T0 <- 8; s <- 0.3; nA <- 500
S <- rbind(make_group(nA, c(2, 2), c(1, 3)), make_group(nA, c(2, 5), c(3, 6)))
X <- S + matrix(rnorm(length(S), 0, s), nrow(S))
pre <- seq_len(TT) <= T0; post <- !pre
r <- 6

targets <- sample(nA, 30)                                 # placebo targets from group A
res <- t(sapply(targets, function(i) {
  donors <- X[-i, ]
  cl <- cluster_donors(donors[, pre], r = r, k_max = 5)
  g <- assign_cluster(X[i, pre], cl)
  A <- donors[cl$labels == g, , drop = FALSE]
  c(all = mean((S[i, post] - robust_sc(donors, X[i, ], pre, post, r)$fit[post])^2),
    cluster = mean((S[i, post] - robust_sc(A, X[i, ], pre, post, r)$fit[post])^2),
    k = cl$k, n_cluster = nrow(A),
    share_A = mean((setdiff(seq_len(nrow(X)), i)[cl$labels == g]) <= nA))
}))
cat(sprintf("median post MSE: all donors %.4f | ClusterSC %.4f | k chosen: %s | cluster size (median) %d | share of group A in chosen cluster (median) %.2f\n",
            median(res[, "all"]), median(res[, "cluster"]), paste(unique(res[, "k"]), collapse = ","),
            as.integer(median(res[, "n_cluster"])), median(res[, "share_A"])))
stopifnot(median(res[, "cluster"]) < median(res[, "all"]), median(res[, "share_A"]) > 0.9)
cat("PASS ClusterSC beats all-donor SC on the paper's DGP and picks group A\n")
