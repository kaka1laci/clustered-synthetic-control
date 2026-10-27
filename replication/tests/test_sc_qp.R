# Checks for sc_qp.R.  Usage: Rscript --vanilla replication/tests/test_sc_qp.R <src_dir> <panel_rds>
args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "sc_qp.R"))

# 1. Exact convex combination: weights and a known effect are recovered.
set.seed(1)
N0 <- 60; T0 <- 40; T1 <- 5
Y0 <- matrix(rnorm(N0 * (T0 + T1), sd = 0.02), N0)
w_true <- c(0.5, 0.3, 0.2, rep(0, N0 - 3))
y1 <- as.vector(crossprod(w_true, Y0))
y1[(T0 + 1):(T0 + T1)] <- y1[(T0 + 1):(T0 + T1)] + 0.03
fit <- sc_estimate_qp(rbind(Y0, y1), N0, T0)
stopifnot(abs(fit$estimate - 0.03) < 1e-6, max(abs(fit$weights - w_true)) < 1e-4,
          abs(sum(fit$weights) - 1) < 1e-12, all(fit$weights >= 0))
cat("PASS exact convex combination: estimate", fit$estimate, "\n")

# 2. Geithner, schedule connections: compare with synthdid::sc_estimate (0.0161 from our run).
p <- readRDS(args[2])
D <- p$treat$schedule
Ys <- p$Y[c(which(D == 0), which(D == 1)), ]
tm <- system.time(fit <- sc_estimate_qp(Ys, sum(D == 0), sum(p$dif < 0)))
pre_rmse <- sqrt(mean(fit$gap[p$dif < 0]^2))
cat(sprintf("Geithner schedule SC (QP): %.4f | pre-period RMSE %.5f | donors with w > 0.01: %d | %.2fs\n",
            fit$estimate, pre_rmse, sum(fit$weights > 0.01), tm["elapsed"]))
stopifnot(abs(fit$estimate - 0.0161) < 0.001)
cat("PASS matches synthdid sc_estimate within 0.001\n")

# 3. The working-set solver reaches the same optimum as the full QP.
pre <- p$dif < 0; N0g <- sum(D == 0)
Y0p <- Ys[seq_len(N0g), pre]; y1p <- colMeans(Ys[-seq_len(N0g), pre])
zeta <- 1e-6 * sd(apply(Y0p, 1, diff))
w_full <- sc_weights_qp(Y0p, y1p, zeta, working_set = FALSE)
tm <- system.time(w_ws <- sc_weights_qp(Y0p, y1p, zeta, working_set = TRUE, ws_init = 50, ws_add = 25))
obj <- function(w) mean((crossprod(Y0p, w) - y1p)^2)
cat(sprintf("Working set vs full QP: max |w diff| %.2e | objective %.10f vs %.10f | %.2fs\n",
            max(abs(w_ws - w_full)), obj(w_ws), obj(w_full), tm["elapsed"]))
stopifnot(max(abs(w_ws - w_full)) < 1e-6, abs(obj(w_ws) - obj(w_full)) < 1e-12)
cat("PASS working-set solution equals full QP\n")
