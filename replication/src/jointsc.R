# 方案 B: joint donor-group selection and synthetic control weights, solved with Gurobi.
#
# 做法 1 (jsc1_fit): one target; choose a group of m donors and simplex weights on it jointly:
#   min_{w,z} (1/T0) ||y_pre - X_pre' w||^2 + lambda * sum_j z_j d_j
#   s.t. w >= 0, sum(w) = 1, w_j <= z_j, z_j in {0, 1}, sum(z) = m                       (convex MIQP)
#   d_j = squared feature distance between donor j and the target.
# 做法 2 (jsc2_fit): several targets share one partition of the donors into K groups; each target belongs
#   to one group and may only use that group's donors. Given centres c_k, step 1 solves
#   min_{a,b,w} sum_t (1/T0) ||y_t,pre - X_pre' w_t||^2 + gamma * sum_{j,k} a_jk ||f_j - c_k||^2
#   s.t. sum_k a_jk = 1, sum_k b_tk = 1, sum_j w_tj = 1, w_tj <= a_jk + 1 - b_tk          (convex MIQP)
#   and step 2 moves each centre to the mean of its group's features; the two steps alternate as in k-means.
# Both MIQPs are solved in units of the donors' pre-period return sd, so decimal and percent returns give
# the same solution; objective values are reported in the input units.
# Needs sc_qp.R (sc_weights_qp) and clustersc.R (choose_k) sourced first.

suppressMessages({ library(gurobi); library(Matrix) })

# ---- features and helpers ---------------------------------------------------------------

# Factor loadings: OLS of (R - rf) on an intercept and the factors over the pre-period; n x p betas.
features_ff3 <- function(R_pre, fac_pre, rf_pre = 0) {
  if (is.null(dim(R_pre))) R_pre <- matrix(R_pre, 1)
  X <- cbind(1, as.matrix(fac_pre))
  ex <- R_pre - matrix(rf_pre, nrow(R_pre), ncol(R_pre), byrow = TRUE)
  B <- solve(crossprod(X), crossprod(X, t(ex)))
  t(B[-1, , drop = FALSE])
}

# Standardize donor and target features by the donors' mean and sd.
standardize_features <- function(F_donor, F_target) {
  mu <- colMeans(F_donor); s <- apply(F_donor, 2, sd)
  sc <- function(M) sweep(sweep(rbind(M), 2, mu), 2, s, "/")
  list(donor = sc(F_donor), target = sc(F_target))
}

feature_dist <- function(F_donor, f) colSums((t(F_donor) - as.vector(f))^2)

kmeans_features <- function(F, k = NULL, k_max = 10, seed = 1) {
  if (is.null(k)) k <- choose_k(F, k_max, seed = seed)$k
  set.seed(seed)
  km <- kmeans(F, centers = k, nstart = 20, iter.max = 100)
  list(labels = km$cluster, centers = km$centers, k = k)
}

# Nearest centre among the allowed groups (e.g. the groups that still have donors).
nearest_center <- function(f, centers, allowed = seq_len(nrow(centers))) {
  dd <- colSums((t(centers[allowed, , drop = FALSE]) - as.vector(f))^2)
  allowed[which.min(dd)]
}

# Candidates ordered for a 做法 1 start: the target's cluster g by distance, then the rest by distance;
# the start group for size m is the first m.
start_order_from_clusters <- function(d, labels, g) {
  inside <- which(labels == g); outside <- which(labels != g)
  c(inside[order(d[inside])], outside[order(d[outside])])
}

# Pre-screen: the n_keep donors with the highest pre-period return correlation with the target (with several
# targets, the highest average correlation). Zero-variance donors have no correlation and come last.
prescreen_corr <- function(X_pre, Ytarget_pre, n_keep) {
  r <- suppressWarnings(cor(t(X_pre), t(rbind(Ytarget_pre))))
  head(order(rowMeans(r), decreasing = TRUE, na.last = TRUE), n_keep)
}

# Adjusted Rand index of two labelings.
ari <- function(x, y) {
  tab <- table(x, y)
  comb2 <- function(v) sum(v * (v - 1) / 2)
  s_ij <- comb2(tab); s_a <- comb2(rowSums(tab)); s_b <- comb2(colSums(tab))
  expected <- s_a * s_b / comb2(length(x))
  (s_ij - expected) / ((s_a + s_b) / 2 - expected)
}

# ---- Gurobi and SC helpers ---------------------------------------------------------------

GRB_DEFAULTS <- list(OutputFlag = 0, MIPGap = 0.01, Threads = 2)
grb_params <- function(time_limit, params) modifyList(c(GRB_DEFAULTS, list(TimeLimit = time_limit)), params)
grb_value <- function(res, field) if (is.null(res[[field]])) NA_real_ else res[[field]]

# Simplex SC weights on a subset of donors, returned as a full-length weight vector.
subset_sc_weights <- function(Y0, y, pre, members) {
  stopifnot(length(members) >= 1)
  w <- numeric(nrow(Y0))
  w[members] <- sc_weights_qp(Y0[members, pre, drop = FALSE], y[pre], 0)
  w
}

pre_mse <- function(Y0, y, pre, w) mean((y[pre] - as.vector(crossprod(Y0[, pre, drop = FALSE], w)))^2)

# lambda = c * lambda_scale: SC's pre-period MSE on all candidates over the median distance, so c = 1
# makes the two terms of the objective comparable.
lambda_scale <- function(Y0, y, pre, d) pre_mse(Y0, y, pre, subset_sc_weights(Y0, y, pre, seq_len(nrow(Y0)))) / median(d)

# ---- 做法 1 ----------------------------------------------------------------------------------

# Y0: n x T candidate donors; y: length-T target; pre / post: logical over columns; d: length-n distances;
# start: optional m donor rows for the MIP start (weights = simplex SC on them).
jsc1_fit <- function(Y0, y, pre, post, d, lambda, m, start = NULL, params = list()) {
  n <- nrow(Y0); T0 <- sum(pre)
  if (m > n) stop(sprintf("jsc1_fit: m = %d exceeds the %d candidate donors", m, n))
  stopifnot(m >= 1, length(d) == n, length(y) == ncol(Y0), all(is.finite(d)), lambda >= 0,
            is.null(start) || length(start) == m)
  s <- sd(as.vector(Y0[, pre]))                                    # work in units of the donors' sd
  X <- Y0[, pre, drop = FALSE] / s
  yp <- y[pre] / s
  # variables: w (n) | z (n) | e (T0), with e = y_pre - X'w
  iw <- seq_len(n); iz <- n + iw; ie <- 2 * n + seq_len(T0); nv <- 2 * n + T0
  r_fit <- 2 + n + seq_len(T0)
  # rows: sum(w) = 1 | sum(z) = m | w_j - z_j <= 0 (n) | e_t + sum_j X_jt w_j = y_t (T0)
  A <- sparseMatrix(
    i = c(rep(1, n), rep(2, n), 2 + iw, 2 + iw, r_fit, rep(r_fit, each = n)),
    j = c(iw, iz, iw, iz, ie, rep(iw, T0)),
    x = c(rep(1, 3 * n), rep(-1, n), rep(1, T0), as.vector(X)),
    dims = c(2 + n + T0, nv))
  model <- list(A = A, rhs = c(1, m, rep(0, n), yp), sense = c("=", "=", rep("<", n), rep("=", T0)),
                obj = c(rep(0, n), lambda / s^2 * d, rep(0, T0)),
                Q = sparseMatrix(i = ie, j = ie, x = rep(1 / T0, T0), dims = c(nv, nv)),
                lb = c(rep(0, 2 * n), rep(-Inf, T0)), ub = c(rep(1, 2 * n), rep(Inf, T0)),
                vtype = c(rep("C", n), rep("B", n), rep("C", T0)), modelsense = "min")
  if (!is.null(start)) {
    w0 <- subset_sc_weights(Y0, y, pre, start)
    z0 <- numeric(n); z0[start] <- 1
    model$start <- c(w0, z0, yp - as.vector(crossprod(X, w0)))
  }
  res <- gurobi(model, grb_params(60, params))
  if (is.null(res$x)) stop(sprintf("jsc1_fit: Gurobi returned no solution (status %s)", res$status))
  w <- pmax(res$x[iw], 0); w <- w / sum(w)
  gap <- y - as.vector(crossprod(Y0, w))
  list(estimate = mean(gap[post]), gap = gap, weights = w, selected = which(res$x[iz] > 0.5),
       pre_rmse = sqrt(mean(gap[pre]^2)), objval = res$objval * s^2, status = res$status,
       mipgap = grb_value(res, "mipgap"), runtime = res$runtime, lambda = lambda, m = m)
}

# Tune (m, c) on a hold-out: fit on the first (1 - holdout) of the pre-period, score on the rest, refit the
# best pair on the whole pre-period. lambda = c * lambda_scale on the window being fitted.
# start_order: candidate rows in start order (start_order_from_clusters); the start for size m is its first m.
jsc1_tune <- function(Y0, y, pre, post, d, start_order, m_grid = c(20, 50, 100), c_grid = c(0, 0.1, 1, 10),
                      holdout = 0.2, params = list()) {
  pidx <- which(pre)
  n_tr <- floor((1 - holdout) * length(pidx))
  tr <- seq_along(pre) %in% pidx[seq_len(n_tr)]
  va <- seq_along(pre) %in% pidx[-seq_len(n_tr)]
  scale_tr <- lambda_scale(Y0, y, tr, d)
  grid <- expand.grid(m = m_grid, c = c_grid)
  fits <- lapply(seq_len(nrow(grid)), function(g)
    jsc1_fit(Y0, y, tr, va, d, grid$c[g] * scale_tr, grid$m[g], start_order[seq_len(grid$m[g])], params))
  grid$val_mse <- sapply(fits, function(f) mean(f$gap[va]^2))
  grid$status <- sapply(fits, `[[`, "status")
  grid$mipgap <- sapply(fits, `[[`, "mipgap")
  grid$runtime <- sapply(fits, `[[`, "runtime")
  best <- which.min(grid$val_mse)
  fit <- jsc1_fit(Y0, y, pre, post, d, grid$c[best] * lambda_scale(Y0, y, pre, d), grid$m[best],
                  start_order[seq_len(grid$m[best])], params)
  fit$c <- grid$c[best]
  fit$tuning <- grid
  fit$runtime_total <- sum(grid$runtime) + fit$runtime
  fit
}

# ---- 做法 2 ----------------------------------------------------------------------------------

# gamma = c * (mean pre-period MSE of the targets' SC within their initial groups)
#           / (mean k-means cost per donor under the initial partition)
jsc2_gamma <- function(Y0, Ytargets, pre, F, init, c = 1) {
  mse <- sapply(seq_len(nrow(Ytargets)), function(t) {
    w <- subset_sc_weights(Y0, Ytargets[t, ], pre, which(init$labels == init$target_groups[t]))
    pre_mse(Y0, Ytargets[t, ], pre, w)
  })
  c * mean(mse) / mean(rowSums((F - init$centers[init$labels, , drop = FALSE])^2))
}

jsc2_objective <- function(Y0, Ytargets, pre, F, W, labels, centers, gamma) {
  fit <- sum(sapply(seq_len(nrow(Ytargets)), function(t) pre_mse(Y0, Ytargets[t, ], pre, W[t, ])))
  fit + gamma * sum((F - centers[labels, , drop = FALSE])^2)
}

# Y0: n x T donors; Ytargets: N1 x T targets; F: n x p donor features;
# init: list(labels = length n in 1..K, centers = K x p, target_groups = length N1).
jsc2_fit <- function(Y0, Ytargets, pre, post, F, K, gamma, init, max_iter = 10, params = list()) {
  n <- nrow(Y0); N1 <- nrow(Ytargets); T0 <- sum(pre)
  stopifnot(nrow(F) == n, length(init$labels) == n, nrow(init$centers) == K,
            length(init$target_groups) == N1, gamma >= 0)
  if (!all(init$target_groups %in% init$labels))
    stop("jsc2_fit: every target's initial group must contain at least one donor")
  s <- sd(as.vector(Y0[, pre]))
  X <- Y0[, pre, drop = FALSE] / s
  Yt <- Ytargets[, pre, drop = FALSE] / s
  # variables: w (t, j) | a (j, k) | b (t, k) | e (t, u), with e_t = y_t,pre - X'w_t
  nw <- N1 * n; na <- n * K; nb <- N1 * K; ne <- N1 * T0; nv <- nw + na + nb + ne
  iw <- function(t, j) (t - 1) * n + j
  ia <- function(j, k) nw + (k - 1) * n + j
  ib <- function(t, k) nw + na + (k - 1) * N1 + t
  ie <- function(t, u) nw + na + nb + (t - 1) * T0 + u
  # rows: sum_k a_jk = 1 (n) | sum_k b_tk = 1 (N1) | sum_j w_tj = 1 (N1) |
  #       w_tj - a_jk + b_tk <= 1 (N1 n K) | e_tu + sum_j X_ju w_tj = y_tu (N1 T0)
  r2 <- n; r3 <- n + N1; r4 <- n + 2 * N1; r5 <- r4 + N1 * n * K
  jj <- rep(seq_len(n), K); kk <- rep(seq_len(K), each = n)
  tb <- rep(seq_len(N1), K); kb <- rep(seq_len(K), each = N1)
  tw <- rep(seq_len(N1), each = n); jw <- rep(seq_len(n), N1)
  tl <- rep(rep(seq_len(N1), each = n), K); jl <- rep(seq_len(n), N1 * K); kl <- rep(seq_len(K), each = N1 * n)
  rl <- r4 + seq_len(N1 * n * K)
  tf <- rep(seq_len(N1), each = n * T0); uf <- rep(rep(seq_len(T0), each = n), N1); jf <- rep(seq_len(n), T0 * N1)
  te <- rep(seq_len(N1), each = T0); ue <- rep(seq_len(T0), N1)
  A <- sparseMatrix(
    i = c(jj, r2 + tb, r3 + tw, rl, rl, rl, r5 + (tf - 1) * T0 + uf, r5 + (te - 1) * T0 + ue),
    j = c(ia(jj, kk), ib(tb, kb), iw(tw, jw), iw(tl, jl), ia(jl, kl), ib(tl, kl), iw(tf, jf), ie(te, ue)),
    x = c(rep(1, na + nb + nw), rep(1, length(rl)), rep(-1, length(rl)), rep(1, length(rl)),
          rep(as.vector(X), N1), rep(1, ne)),
    dims = c(r5 + ne, nv))
  model <- list(A = A, rhs = c(rep(1, n + 2 * N1 + N1 * n * K), as.vector(t(Yt))),
                sense = c(rep("=", n + 2 * N1), rep("<", N1 * n * K), rep("=", ne)),
                Q = sparseMatrix(i = nw + na + nb + seq_len(ne), j = nw + na + nb + seq_len(ne),
                                 x = rep(1 / T0, ne), dims = c(nv, nv)),
                lb = c(rep(0, nw + na + nb), rep(-Inf, ne)), ub = c(rep(1, nw + na + nb), rep(Inf, ne)),
                vtype = c(rep("C", nw), rep("B", na + nb), rep("C", ne)), modelsense = "min")
  start_vec <- function(W, labels, groups) {
    Am <- matrix(0, n, K); Am[cbind(seq_len(n), labels)] <- 1
    Bm <- matrix(0, N1, K); Bm[cbind(seq_len(N1), groups)] <- 1
    c(as.vector(t(W)), as.vector(Am), as.vector(Bm), as.vector(t(Yt - W %*% X)))
  }

  labels <- init$labels; groups <- init$target_groups; centers <- init$centers
  W <- t(sapply(seq_len(N1), function(t) subset_sc_weights(Y0, Ytargets[t, ], pre, which(labels == groups[t]))))
  obj_init <- jsc2_objective(Y0, Ytargets, pre, F, W, labels, centers, gamma)
  history <- data.frame()
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    # step 1: given the centres, the partition and all targets' weights jointly (warm start = current state)
    D <- sapply(seq_len(K), function(k) colSums((t(F) - centers[k, ])^2))           # n x K
    model$obj <- c(rep(0, nw), gamma / s^2 * as.vector(D), rep(0, nb + ne))
    model$start <- start_vec(W, labels, groups)
    res <- gurobi(model, grb_params(300, params))
    if (is.null(res$x)) stop(sprintf("jsc2_fit: Gurobi returned no solution in iteration %d (status %s)", it, res$status))
    new_labels <- apply(matrix(res$x[nw + seq_len(na)], n, K), 1, which.max)
    new_groups <- apply(matrix(res$x[nw + na + seq_len(nb)], N1, K), 1, which.max)
    W <- matrix(pmax(res$x[seq_len(nw)], 0), N1, n, byrow = TRUE)
    W[outer(new_groups, new_labels, `!=`)] <- 0                        # only the target's own group
    W <- W / rowSums(W)
    obj1 <- jsc2_objective(Y0, Ytargets, pre, F, W, new_labels, centers, gamma)
    # step 2: centres = group means of the features; an empty group keeps its centre
    old_centers <- centers
    for (k in seq_len(K)) if (any(new_labels == k)) centers[k, ] <- colMeans(F[new_labels == k, , drop = FALSE])
    obj2 <- jsc2_objective(Y0, Ytargets, pre, F, W, new_labels, centers, gamma)
    # stop only at a fixed point: moved centres need another step 1 even when the partition did not change
    # (e.g. start centres from a larger set of firms than the candidates)
    changed <- any(new_labels != labels) || any(new_groups != groups) ||
               max(abs(centers - old_centers)) > 1e-10 * max(1, abs(old_centers))
    history <- rbind(history, data.frame(iter = it, obj_step1 = obj1, obj_step2 = obj2, status = res$status,
                                         mipgap = grb_value(res, "mipgap"), runtime = res$runtime,
                                         donors_moved = sum(new_labels != labels)))
    labels <- new_labels; groups <- new_groups
    if (!changed) { converged <- TRUE; break }
  }
  gaps <- Ytargets - W %*% Y0
  list(estimates = rowMeans(gaps[, post, drop = FALSE]), gaps = gaps, weights = W, labels = labels,
       target_groups = groups, centers = centers, history = history, obj_init = obj_init,
       converged = converged, n_iter = nrow(history), gamma = gamma, runtime = sum(history$runtime),
       mipgap = max(history$mipgap))
}
