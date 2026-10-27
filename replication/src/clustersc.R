# ClusterSC (Rho, Tang, Bergam, Cummings & Misra, AISTATS 2025), adapted to a single target series.
#
#   1. SVD of the donors' pre-period matrix; features = rows of U_r Sigma_r (= X V_r)
#   2. k-means on the features (k by silhouette unless given)
#   3. Project the target's pre-period series on V_r and assign it to the nearest centre
#   4. SC on the donors of that cluster only:
#        robust SC  = hard singular value thresholding (HSVT) + ridge, as in the paper; or
#        simplex SC = the synthdid sc objective (sc_qp.R), as in GP&L
#
# Differences from the paper (see paper-card-clustersc.md): GP&L match the treated cohort's average
# return, so the "target" here is that average series; and daily returns are far from low rank, so the
# rank r is an explicit argument instead of the paper's 95%-of-singular-values rule.

suppressMessages(library(cluster))

svd_features <- function(X_pre, r) {
  sv <- svd(X_pre, nu = 0, nv = r)
  list(V = sv$v, features = X_pre %*% sv$v, d = sv$d)
}

choose_k <- function(features, k_max = 10, n_sil = 2000, seed = 1) {
  set.seed(seed)
  idx <- if (nrow(features) > n_sil) sample(nrow(features), n_sil) else seq_len(nrow(features))
  dm <- dist(features[idx, , drop = FALSE])
  sil <- sapply(2:k_max, function(k) {
    km <- kmeans(features, centers = k, nstart = 10, iter.max = 100)
    mean(silhouette(km$cluster[idx], dm)[, "sil_width"])
  })
  list(k = (2:k_max)[which.max(sil)], silhouette = setNames(sil, 2:k_max))
}

cluster_donors <- function(X_pre, r, k = NULL, k_max = 10, seed = 1) {
  f <- svd_features(X_pre, r)
  ks <- if (is.null(k)) choose_k(f$features, k_max, seed = seed) else list(k = k, silhouette = NA)
  set.seed(seed)
  km <- kmeans(f$features, centers = ks$k, nstart = 20, iter.max = 100)
  list(labels = km$cluster, centers = km$centers, V = f$V, k = ks$k, silhouette = ks$silhouette,
       sizes = as.vector(table(km$cluster)))
}

assign_cluster <- function(x_pre, cl) {
  u <- as.vector(crossprod(cl$V, x_pre))
  which.min(colSums((t(cl$centers) - u)^2))
}

hsvt <- function(M, r) {
  sv <- svd(M, nu = r, nv = r)
  sv$u %*% (sv$d[seq_len(r)] * t(sv$v))
}

# Robust SC: HSVT on the donors' full window, ridge weights fitted on the pre-period, lambda chosen on the
# last `holdout` share of the pre-period. Donors are untreated, so using their post-period data is allowed.
# M: optional precomputed hsvt(Y0, r), so several targets can share one donor set's SVD.
robust_sc <- function(Y0, y1, pre, post, r, lambda_grid = 10^seq(-4, 1), holdout = 0.2, M = NULL) {
  if (is.null(M)) M <- hsvt(Y0, min(r, nrow(Y0), ncol(Y0)))
  Mp <- M[, pre, drop = FALSE]; yp <- y1[pre]
  # (M M' + lam I_n)^{-1} M y = M (M'M + lam I_T)^{-1} y: solve the T x T system, since n >> T
  ridge <- function(Mt, yt, lam) Mt %*% solve(crossprod(Mt) + lam * diag(ncol(Mt)), yt)
  scale <- mean(Mp^2) * ncol(Mp)                          # lambda grid relative to the donors' scale
  n_tr <- floor((1 - holdout) * ncol(Mp))
  tr <- seq_len(n_tr); va <- (n_tr + 1):ncol(Mp)
  val_err <- sapply(lambda_grid, function(l) {
    w <- ridge(Mp[, tr, drop = FALSE], yp[tr], l * scale)
    mean((yp[va] - crossprod(Mp[, va, drop = FALSE], w))^2)
  })
  lam <- lambda_grid[which.min(val_err)] * scale
  w <- ridge(Mp, yp, lam)
  fit <- as.vector(crossprod(M, w))
  list(estimate = mean(y1[post] - fit[post]), pre_rmse = sqrt(mean((yp - fit[pre])^2)), lambda = lam, fit = fit)
}

# Simplex SC on a donor set for a single target series (wraps sc_qp.R)
simplex_sc <- function(Y0, y1, pre, post) {
  fit <- sc_estimate_qp(rbind(Y0, y1), nrow(Y0), sum(pre))
  list(estimate = fit$estimate, pre_rmse = sqrt(mean(fit$gap[pre]^2)), fit = y1 - fit$gap, weights = fit$weights)
}
