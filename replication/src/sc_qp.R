# Synthetic control by quadratic programming.
#
# Solves the same problem as synthdid::sc_estimate (simplex weights, no intercept,
# ridge zeta = eta * noise.level with eta = 1e-6):
#   min_w  (1/T0) * || Y0_pre' w - ybar1_pre ||^2 + zeta^2 * ||w||^2,   w >= 0, sum(w) = 1
# synthdid uses Frank-Wolfe, which takes ~15 minutes per fit on the Geithner panel because
# the near-zero ridge makes it converge slowly; quadprog gets the exact optimum in well under a second.
# This makes the 100-replication placebo standard error feasible.

suppressMessages(library(quadprog))

qp_simplex <- function(Y0_pre, y1_pre, ridge) {
  T0 <- ncol(Y0_pre)
  N0 <- nrow(Y0_pre)
  Dmat <- 2 * (tcrossprod(Y0_pre) / T0 + diag(ridge, N0))
  dvec <- 2 * as.vector(Y0_pre %*% y1_pre) / T0
  Amat <- cbind(rep(1, N0), diag(N0))                # sum(w) = 1, then w >= 0
  w <- solve.QP(Dmat, dvec, Amat, c(1, rep(0, N0)), meq = 1)$solution
  w[w < 0] <- 0
  w / sum(w)
}

sc_weights_qp <- function(Y0_pre, y1_pre, zeta, min_ridge = 1e-10,
                          working_set = nrow(Y0_pre) > 600, ws_init = 300, ws_add = 200) {
  # Y0_pre: N0 x T0 control pre-period outcomes; y1_pre: length-T0 treated average.
  T0 <- ncol(Y0_pre)
  N0 <- nrow(Y0_pre)
  # quadprog needs a strictly positive definite matrix; N0 > T0 makes the Gram matrix singular,
  # so the ridge is floored at min_ridge * mean(diag(Gram)) (negligible next to the fit term).
  ridge <- max(zeta^2, min_ridge * mean(rowSums(Y0_pre^2)) / T0)
  if (!working_set) return(qp_simplex(Y0_pre, y1_pre, ridge))

  # Working-set method for large donor pools: the solution is sparse, so solve the QP on a subset
  # of donors and add any donor that violates the KKT conditions, until none does. The result is
  # the optimum of the full problem; the full N0 x N0 QP is O(N0^3) and too slow for N0 ~ 4,000.
  score <- as.vector(Y0_pre %*% y1_pre) / sqrt(rowSums(Y0_pre^2))
  S <- order(score, decreasing = TRUE)[seq_len(min(N0, ws_init))]
  repeat {
    wS <- qp_simplex(Y0_pre[S, , drop = FALSE], y1_pre, ridge)
    w <- numeric(N0); w[S] <- wS
    resid <- as.vector(crossprod(Y0_pre, w)) - y1_pre
    grad <- 2 * as.vector(Y0_pre %*% resid) / T0 + 2 * ridge * w
    nu <- mean(grad[w > 0])                          # multiplier on sum(w) = 1
    tol <- 1e-9 * max(abs(grad))
    viol <- setdiff(which(grad < nu - tol), S)
    if (length(viol) == 0) return(w)
    S <- c(S, viol[order(grad[viol])][seq_len(min(length(viol), ws_add))])
  }
}

# Y: N x T with controls in rows 1..N0, treated after; columns 1..T0 are pre-period.
sc_estimate_qp <- function(Y, N0, T0, eta = 1e-6) {
  pre  <- seq_len(T0)
  post <- (T0 + 1):ncol(Y)
  Y0   <- Y[seq_len(N0), , drop = FALSE]
  Y1   <- Y[(N0 + 1):nrow(Y), , drop = FALSE]
  noise_level <- sd(apply(Y0[, pre, drop = FALSE], 1, diff))   # synthdid's default
  w <- sc_weights_qp(Y0[, pre, drop = FALSE], colMeans(Y1[, pre, drop = FALSE]), eta * noise_level)
  gap <- colMeans(Y1) - as.vector(crossprod(w, Y0))
  list(estimate = mean(gap[post]), weights = w, gap = gap)
}

# Placebo SE as in synthdid::vcov(method = "placebo"): reassign N1 random controls as
# pseudo-treated, refit on the remaining controls, take the sd of the placebo estimates.
sc_placebo_se_qp <- function(Y, N0, T0, replications = 100, eta = 1e-6) {
  N1 <- nrow(Y) - N0
  stopifnot(N0 > N1)
  est <- replicate(replications, {
    ind <- sample(seq_len(N0))
    sc_estimate_qp(Y[ind, , drop = FALSE], N0 - N1, T0, eta)$estimate
  })
  sqrt((replications - 1) / replications) * sd(est)
}
