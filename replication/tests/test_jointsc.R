# Checks for jointsc.R (方案 B). Usage: Rscript --vanilla replication/tests/test_jointsc.R <src_dir>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(gurobi))

# Test 5. Gurobi smoke test: the package loads, the licence works, problems with known solutions solve.
#   QP:   min (x - 1)^2 + (y - 2)^2    s.t. x + y <= 2, x, y >= 0   ->  (0.5, 1.5), objective 0.5
#   MIQP: min (x - 0.8)^2 + (y - 2)^2  same constraints, x integer  ->  (0, 2),    objective 0.64
qp <- list(Q = diag(2), obj = c(-2, -4), objcon = 5, A = matrix(c(1, 1), 1), rhs = 2, sense = "<", lb = c(0, 0))
r <- gurobi(qp, list(OutputFlag = 0))
stopifnot(r$status == "OPTIMAL", max(abs(r$x - c(0.5, 1.5))) < 1e-6, abs(r$objval - 0.5) < 1e-6)
miqp <- modifyList(qp, list(obj = c(-1.6, -4), objcon = 4.64, vtype = c("I", "C")))
r2 <- gurobi(miqp, list(OutputFlag = 0))
stopifnot(r2$status == "OPTIMAL", max(abs(r2$x - c(0, 2))) < 1e-6, abs(r2$objval - 0.64) < 1e-6)
cat("PASS Gurobi smoke test (QP and MIQP)\n")

source(file.path(args[1], "sc_qp.R"))
source(file.path(args[1], "clustersc.R"))
source(file.path(args[1], "jointsc.R"))

# H1. FF loadings: OLS recovers known loadings; an average target's loadings are the average loadings.
set.seed(2)
fac <- matrix(rnorm(100 * 3), 100); rf <- rep(0.01, 100); B <- matrix(runif(15), 5)
R <- 0.2 + B %*% t(fac) + matrix(rf, 5, 100, byrow = TRUE)
stopifnot(max(abs(features_ff3(R, fac, rf) - B)) < 1e-10,
          max(abs(features_ff3(colMeans(R), fac, rf) - colMeans(B))) < 1e-10)
st <- standardize_features(B, B[1:2, ])
stopifnot(max(abs(colMeans(st$donor))) < 1e-12, max(abs(apply(st$donor, 2, sd) - 1)) < 1e-12,
          max(abs(st$target - st$donor[1:2, ])) < 1e-12)
cat("PASS features_ff3 / standardize_features\n")

# H2. Distances, start order, ARI.
stopifnot(isTRUE(all.equal(feature_dist(rbind(c(0, 0), c(3, 4)), c(0, 0)), c(0, 25))),
          identical(start_order_from_clusters(c(5, 1, 4, 2, 3), c(1, 2, 1, 2, 2), 2), c(2L, 4L, 5L, 3L, 1L)),
          isTRUE(all.equal(ari(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)),
          isTRUE(all.equal(ari(c(1, 1, 2, 2), c(1, 2, 1, 2)), -0.5)))
km <- kmeans_features(rbind(matrix(rnorm(40, 0, 0.1), 20), matrix(rnorm(40, 5, 0.1), 20)), k = 2)
stopifnot(km$k == 2, isTRUE(all.equal(ari(km$labels, rep(1:2, each = 20)), 1)))
cat("PASS feature_dist / start_order_from_clusters / ari / kmeans_features\n")

# RF5. A target whose nearest cluster has no donors left goes to the nearest non-empty one;
#      zero-variance donors are ranked last by the pre-screen; n_keep > n keeps everyone.
ctr <- rbind(c(0, 0), c(10, 0), c(0, 10))
stopifnot(nearest_center(c(6, 1), ctr) == 2, nearest_center(c(6, 1), ctr, allowed = c(1, 3)) == 1)
set.seed(5)
Xs <- rbind(rnorm(50), rep(0, 50), rnorm(50)); ys <- Xs[3, ] + rnorm(50, sd = 0.1)
stopifnot(identical(prescreen_corr(Xs, ys, 3), c(3L, 1L, 2L)), length(prescreen_corr(Xs, ys, 10)) == 3,
          identical(prescreen_corr(Xs, rbind(ys, ys), 1), 3L))
cat("PASS nearest_center (allowed groups) / prescreen_corr\n")

# H3. simplex_sc also returns the weights of sc_estimate_qp.
set.seed(6)
Yh <- matrix(rnorm(20 * 40), 20); yh <- colMeans(Yh[1:3, ]) + rnorm(40, sd = 0.1); preh <- seq_len(40) <= 35
stopifnot(identical(simplex_sc(Yh, yh, preh, !preh)$weights, sc_estimate_qp(rbind(Yh, yh), 20, 35)$weights))
cat("PASS simplex_sc returns weights\n")

# Test 1. lambda = 0, m = n: 做法 1 is simplex SC (sc_qp.R).
set.seed(1)
n <- 30; TT <- 70; T0 <- 60
pre <- seq_len(TT) <= T0; post <- !pre
Y0 <- matrix(rnorm(n * TT), n)
y <- as.vector(crossprod(c(0.5, 0.3, 0.2, rep(0, n - 3)), Y0)) + rnorm(TT, sd = 0.5)
d <- runif(n)
EXACT <- list(MIPGap = 1e-9, BarConvTol = 1e-12, OptimalityTol = 1e-9, FeasibilityTol = 1e-9)
f1 <- jsc1_fit(Y0, y, pre, post, d, lambda = 0, m = n, params = EXACT)
sc <- sc_estimate_qp(rbind(Y0, y), n, T0)
cat(sprintf("test 1: |estimate diff| %.2e | max |w diff| %.2e | objective %.8f vs %.8f\n",
            abs(f1$estimate - sc$estimate), max(abs(f1$weights - sc$weights)), f1$objval, mean(sc$gap[pre]^2)))
stopifnot(abs(f1$estimate - sc$estimate) < 1e-6, max(abs(f1$weights - sc$weights)) < 1e-5,
          abs(f1$objval - mean(sc$gap[pre]^2)) < 1e-6, length(f1$selected) == n)
cat("PASS 做法 1 with lambda = 0, m = n equals simplex SC\n")

# Test 2. Very large lambda: the group is the m donors nearest the target.
d2 <- sample(seq_len(n)) / n
f2 <- jsc1_fit(Y0, y, pre, post, d2, 1e5 * lambda_scale(Y0, y, pre, d2), m = 5, params = EXACT)
stopifnot(setequal(f2$selected, order(d2)[1:5]))
cat("PASS 做法 1 with large lambda picks the m nearest donors\n")

# Test 3. n = 8, m = 3: the MIQP optimum equals brute force over all 56 groups.
set.seed(3)
n3 <- 8; TT3 <- 35; pre3 <- seq_len(TT3) <= 30; post3 <- !pre3
Y3 <- matrix(rnorm(n3 * TT3), n3)
y3 <- as.vector(crossprod(c(0.4, 0.4, 0.2, rep(0, 5)), Y3)) + rnorm(TT3, sd = 0.3)
d3 <- runif(n3); lam3 <- lambda_scale(Y3, y3, pre3, d3)
brute <- apply(combn(n3, 3), 2, function(S)
  pre_mse(Y3, y3, pre3, subset_sc_weights(Y3, y3, pre3, S)) + lam3 * sum(d3[S]))
f3 <- jsc1_fit(Y3, y3, pre3, post3, d3, lam3, m = 3, params = EXACT)
cat(sprintf("test 3: MIQP %.10f vs brute force %.10f over %d groups\n", f3$objval, min(brute), length(brute)))
stopifnot(length(brute) == 56, abs(f3$objval - min(brute)) < 1e-6 * max(1, min(brute)))
cat("PASS 做法 1 optimum equals brute force\n")

# RF1. Units: decimal returns (Geithner) and percent returns (simulation) give the same solution.
f3d <- jsc1_fit(Y3 / 100, y3 / 100, pre3, post3, d3, lam3 / 1e4, m = 3, params = EXACT)
stopifnot(setequal(f3d$selected, f3$selected), max(abs(f3d$weights - f3$weights)) < 1e-6,
          abs(f3d$objval * 1e4 - f3$objval) < 1e-6 * max(1, f3$objval))
cat("PASS 做法 1 is invariant to the return units\n")

# RF2. m larger than the candidate set stops with a clear message.
msg <- tryCatch({ jsc1_fit(Y3, y3, pre3, post3, d3, lam3, m = 9); "no error" }, error = conditionMessage)
stopifnot(grepl("m = 9 exceeds the 8 candidate donors", msg, fixed = TRUE))
cat("PASS 做法 1 rejects m > n\n")

# RF3. A solver stop before optimality (here after the first solution) returns that solution and its status.
set.seed(7)
Y7 <- matrix(rnorm(200 * TT), 200); y7 <- colMeans(Y7[1:10, ]) + rnorm(TT, sd = 0.2); d7 <- runif(200)
f7 <- jsc1_fit(Y7, y7, pre, post, d7, lambda_scale(Y7, y7, pre, d7), m = 20, start = order(d7)[1:20],
               params = list(SolutionLimit = 1))
stopifnot(f7$status %in% c("SOLUTION_LIMIT", "OPTIMAL"), abs(sum(f7$weights) - 1) < 1e-12,
          length(f7$selected) == 20, all(f7$weights[-f7$selected] < 1e-6), is.numeric(f7$mipgap))
cat(sprintf("PASS 做法 1 returns an early-stopped solution (status %s, gap %.3f)\n", f7$status, f7$mipgap))

# H4. Tuning: the chosen (m, c) minimises the hold-out error; the final fit uses it on the whole pre-period.
ft <- jsc1_tune(Y0, y, pre, post, d, order(d), m_grid = c(3, 5), c_grid = c(0, 1), params = list(MIPGap = 1e-6))
best <- which.min(ft$tuning$val_mse)
stopifnot(nrow(ft$tuning) == 4, ft$m == ft$tuning$m[best], ft$c == ft$tuning$c[best],
          length(ft$selected) == ft$m, abs(sum(ft$weights) - 1) < 1e-12,
          isTRUE(all.equal(ft$lambda, ft$c * lambda_scale(Y0, y, pre, d))), ft$runtime_total >= ft$runtime)
cat(sprintf("PASS jsc1_tune picks m = %d, c = %g\n", ft$m, ft$c))

# Test 4. 做法 2 on three well-separated groups, started with 10% of donors in a wrong group:
#              recovers the donor groups and the target groups (ARI = 1); the objective never goes up.
set.seed(4)
K <- 3; n_per <- 20; TT4 <- 60; pre4 <- seq_len(TT4) <= 50; post4 <- !pre4
truth <- rep(seq_len(K), each = n_per)
G <- matrix(rnorm(K * TT4), K)                                              # one factor per group
Y4 <- G[truth, ] * runif(K * n_per, 0.8, 1.2) + matrix(rnorm(K * n_per * TT4, sd = 0.3), K * n_per)
tg <- rep(seq_len(K), each = 2)
Yt4 <- G[tg, ] + matrix(rnorm(length(tg) * TT4, sd = 0.3), length(tg))
F4 <- rbind(c(-3, 0), c(0, 3), c(3, 0))[truth, ] + matrix(rnorm(2 * K * n_per, sd = 0.2), ncol = 2)
wrong <- truth; wrong[c(1, 5, 25, 30, 45, 50)] <- c(2, 3, 3, 1, 2, 1)
init <- list(labels = wrong, centers = t(sapply(seq_len(K), function(k) colMeans(F4[wrong == k, ]))),
             target_groups = tg)
gam <- jsc2_gamma(Y4, Yt4, pre4, F4, init, 1)
f4 <- jsc2_fit(Y4, Yt4, pre4, post4, F4, K, gam, init, params = list(MIPGap = 1e-6))
print(f4$history)
path <- c(f4$obj_init, as.vector(rbind(f4$history$obj_step1, f4$history$obj_step2)))
stopifnot(isTRUE(all.equal(ari(f4$labels, truth), 1)), isTRUE(all.equal(ari(f4$target_groups, tg), 1)),
          f4$converged, all(diff(path) <= 1e-6 * abs(path[1])),
          all(abs(rowSums(f4$weights) - 1) < 1e-12),
          all(f4$weights[outer(f4$target_groups, f4$labels, `!=`)] == 0))
cat(sprintf("PASS 做法 2 recovers the groups (ARI 1) in %d iterations; objective monotone\n", f4$n_iter))

# RF4. An empty group keeps its centre (no NaN) and no target is assigned to it.
init5 <- list(labels = truth, centers = rbind(init$centers, c(100, 100)), target_groups = tg)
f5 <- jsc2_fit(Y4, Yt4, pre4, post4, F4, 4, gam, init5, params = list(MIPGap = 1e-6))
stopifnot(all(is.finite(f5$centers)), all(f5$centers[4, ] == c(100, 100)), !any(f5$labels == 4),
          !any(f5$target_groups == 4), f5$converged)
cat("PASS 做法 2 keeps an empty group's centre\n")

# Review I1. Convergence is a fixed point of the alternation: when the start centres are not the means of the start
# labels (Geithner: centres from all untreated firms, labels on the pre-screened candidates), 做法 2 must re-solve
# against the updated centres before it stops.
set.seed(1)
n6 <- 30; TT6 <- 40; pre6 <- seq_len(TT6) <= 30; post6 <- !pre6
Y6 <- matrix(rnorm(n6 * TT6), n6)
Yt6 <- rbind(colMeans(Y6[1:5, ]) + rnorm(TT6, sd = 0.3), colMeans(Y6[20:25, ]) + rnorm(TT6, sd = 0.3))
F6 <- matrix(c(seq(0, 1, length.out = 14), 1.3, seq(2.5, 3.5, length.out = 15)), ncol = 1)
init6 <- list(labels = c(rep(1, 14), rep(2, 16)), centers = matrix(c(0.5, 2.0), ncol = 1), target_groups = c(1, 2))
gam6 <- jsc2_gamma(Y6, Yt6, pre6, F6, init6, 1)
f6 <- jsc2_fit(Y6, Yt6, pre6, post6, F6, 2, gam6, init6, params = list(MIPGap = 1e-6))
f6b <- jsc2_fit(Y6, Yt6, pre6, post6, F6, 2, gam6,
                list(labels = f6$labels, centers = f6$centers, target_groups = f6$target_groups), params = list(MIPGap = 1e-6))
cat(sprintf("I1: first run %d iterations; restart from its state moves %d donors\n", f6$n_iter, sum(f6b$labels != f6$labels)))
stopifnot(f6$converged, f6$n_iter >= 2, all(f6b$labels == f6$labels), all(f6b$target_groups == f6$target_groups),
          f6b$n_iter == 1)
cat("PASS 做法 2 stops only at a fixed point (start centres off the group means)\n")
