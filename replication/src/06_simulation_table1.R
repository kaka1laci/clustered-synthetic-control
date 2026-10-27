# Simulation of GP&L Table 1 (Section 4): bias and rejection rates of event-study estimators when
# returns follow a two-factor model (market + SMB) and the "CAPM" estimator omits SMB.
#
# DGP as described in the paper, with the ambiguous parts pinned down (assumptions marked [A]):
#   r_it = rf_t + b_mkt_i * (mkt_t - rf_t) + b_smb_i * smb_t + e_it      (percent units)
#   b_mkt, b_smb ~ N(1, 0.3^2) independent; 500 firms; t = 1..250, event at t = 240
#   (239 pre periods, event day, 10 post periods); treatment effect = 3 on the event day only.
#   Factors and rf: a block of 250 consecutive days drawn from FF daily data, July 1926 - Dec 2022.
#   [A] e_it ~ N(0, sigma^2) with sigma set so the residual variance is half of the average total
#       variance (the paper states this target; its literal "0.1" contradicts it in either unit).
#   Assignment: random 10% of firms, or P(treated) = logistic(delta * b_smb), delta = log(0.1)/E(b_smb).
#   Timing: event day = t = 240 of the random block, or timing selection:
#     rank SMB over the 250 days (rank 1 = highest), P(selected_t) = logistic(delta * rank_t),
#     delta = log(1/250)/E(rank); draw each day; among drawn days take the highest SMB.
#     [A] "Rank2" in the paper is read as rank_t. [A] The selected day's factor row is swapped into t = 240.
# Estimators (per event-time period t = 240..250):
#   Simple Means : mean(treated) - mean(control)                         two-sample t-test
#   CAPM         : mean over treated of AR = excess r - a - b * (mkt-rf)  one-sample t-test across treated
#   Correct      : same with mkt-rf and smb                               one-sample t-test across treated
#   Gsynth       : gsynth 1.2.1, force = "unit", CV r in 0..5           parametric-bootstrap CI
# Metrics (reverse-engineered from the paper's Table 1, see README):
#   Treated   : event-day estimate vs 3          -> E(Bias), MAD = mean |error|, rejection rate
#   Untreated : average of the 10 post-day estimates vs 0 -> E(Bias), MAD; rejection = per-day rate
#   All       : average of the 11 estimates vs 3/11 -> E(Bias), MAD; RMSE of the cumulative (sum) error
#   [A] The paper's "Coverage" column is read as the rejection rate (1 = power on the event day).
#
# Usage: Rscript --vanilla replication/src/06_simulation_table1.R <ff_csv> <out_rds> <panel A|B|C|D> \
#          [n_sims=50] [gsynth_boots=200] [seed=2025]

args <- commandArgs(trailingOnly = TRUE)
ff_csv       <- args[1]
out_rds      <- args[2]
panel        <- args[3]
n_sims       <- if (length(args) >= 4) as.integer(args[4]) else 50L
gsynth_boots <- if (length(args) >= 5) as.integer(args[5]) else 200L
seed         <- if (length(args) >= 6) as.integer(args[6]) else 2025L
stopifnot(panel %in% c("A", "B", "C", "D"))
assign_select <- panel %in% c("B", "D")
timing_select <- panel %in% c("C", "D")

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
gsynth_lib <- normalizePath(file.path(script_dir, "..", "rlib", "gsynth121"), mustWork = FALSE)
if (dir.exists(gsynth_lib)) .libPaths(c(gsynth_lib, .libPaths()))
suppressMessages(library(gsynth))

N <- 500; TT <- 250; T_EVENT <- 240; EFFECT <- 3; N_TREAT <- 50
post_t <- T_EVENT:TT

# ---- Fama-French daily factors, percent ------------------------------------
ff_lines <- readLines(ff_csv)
ff <- read.csv(text = ff_lines[grep("^,Mkt-RF", ff_lines):length(ff_lines)], strip.white = TRUE)
names(ff) <- c("date", "mkt_rf", "smb", "hml", "rf")
ff <- ff[grepl("^[0-9]{8}$", ff$date), ]
ff$date <- as.Date(as.character(ff$date), "%Y%m%d")
ff[, -1] <- lapply(ff[, -1], as.numeric)
ff <- ff[ff$date >= as.Date("1926-07-01") & ff$date <= as.Date("2022-12-31"), ]

# [A] sigma: residual variance = half of average total variance, i.e. sigma^2 = E[b' Sigma_F b]
Sig <- cov(ff[, c("mkt_rf", "smb")])
factor_var <- sum(diag(Sig)) * (1 + 0.3^2) + 2 * Sig[1, 2]
SIGMA <- sqrt(factor_var)
cat(sprintf("Panel %s | sims %d | gsynth boots %d | sigma_e = %.3f pp\n", panel, n_sims, gsynth_boots, SIGMA))

logistic <- function(x) exp(x) / (1 + exp(x))

one_sim <- function(s) {
  set.seed(seed + 1000 * match(panel, LETTERS) + s)
  start <- sample.int(nrow(ff) - TT + 1, 1)
  F <- ff[start:(start + TT - 1), c("mkt_rf", "smb", "rf")]
  if (timing_select) {
    rk <- rank(-F$smb, ties.method = "first")
    delta <- log(1 / TT) / mean(rk)
    drawn <- which(runif(TT) < logistic(delta * rk))
    sel <- if (length(drawn)) drawn[which.max(F$smb[drawn])] else which.max(F$smb)
    F[c(sel, T_EVENT), ] <- F[c(T_EVENT, sel), ]
  }
  b_mkt <- rnorm(N, 1, 0.3); b_smb <- rnorm(N, 1, 0.3)
  if (assign_select) {
    treat <- runif(N) < logistic(log(0.1) / 1 * b_smb)      # E(b_smb) = 1
    if (!any(treat)) treat[which.min(b_smb)] <- TRUE
  } else {
    treat <- seq_len(N) %in% sample.int(N, N_TREAT)
  }
  R <- outer(rep(1, N), F$rf) + outer(b_mkt, F$mkt_rf) + outer(b_smb, F$smb) +
       matrix(rnorm(N * TT, 0, SIGMA), N, TT)
  R[treat, T_EVENT] <- R[treat, T_EVENT] + EFFECT
  truth <- c(EFFECT, rep(0, length(post_t) - 1))

  est <- list(); rej <- list()
  # Simple means
  est$SimpleMeans <- sapply(post_t, function(t) mean(R[treat, t]) - mean(R[!treat, t]))
  rej$SimpleMeans <- sapply(seq_along(post_t), function(k) {
    t <- post_t[k]; t.test(R[treat, t], R[!treat, t], var.equal = TRUE)$p.value < 0.05 })
  # Abnormal-return estimators on the treated firms (Definition 3)
  ar_est <- function(X) {
    Xp <- cbind(1, X[1:(T_EVENT - 1), , drop = FALSE])
    Yex <- R[treat, , drop = FALSE] - outer(rep(1, sum(treat)), F$rf)
    B <- solve(crossprod(Xp), crossprod(Xp, t(Yex[, 1:(T_EVENT - 1), drop = FALSE])))
    AR <- Yex[, post_t, drop = FALSE] - t(cbind(1, X[post_t, , drop = FALSE]) %*% B)
    list(est = colMeans(AR),
         rej = apply(AR, 2, function(a) length(a) > 1 && t.test(a)$p.value < 0.05))
  }
  capm <- ar_est(as.matrix(F[, "mkt_rf", drop = FALSE]))
  corr <- ar_est(as.matrix(F[, c("mkt_rf", "smb")]))
  est$CAPM <- capm$est; rej$CAPM <- capm$rej
  est$Correct <- corr$est; rej$Correct <- corr$rej
  # Gsynth
  long <- data.frame(id = rep(seq_len(N), times = TT), time = rep(seq_len(TT), each = N),
                     y = as.vector(R), d = as.vector(outer(treat, seq_len(TT) >= T_EVENT)))
  invisible(capture.output(gs <- suppressMessages(gsynth(y ~ d, data = long, index = c("id", "time"),
              force = "unit", CV = TRUE, r = c(0, 5), se = gsynth_boots > 0, inference = "parametric",
              nboots = max(gsynth_boots, 1L), parallel = FALSE))))
  att <- if (gsynth_boots > 0) gs$est.att else cbind(ATT = gs$att)
  post_rows <- as.numeric(rownames(att)) >= 1
  est$Gsynth <- unname(att[post_rows, "ATT"])
  rej$Gsynth <- if (gsynth_boots > 0) unname(att[post_rows, "CI.lower"] > 0 | att[post_rows, "CI.upper"] < 0) else rep(NA, length(post_t))
  stopifnot(length(est$Gsynth) == length(post_t))

  list(est = est, rej = rej, truth = truth, n_treat = sum(treat), smb_event = F$smb[T_EVENT],
       gsynth_r = gs$r.cv)
}

t0 <- Sys.time()
sims <- vector("list", n_sims)
for (s in seq_len(n_sims)) {
  sims[[s]] <- one_sim(s)
  if (s %% 10 == 0 || s == n_sims)
    cat(sprintf("[%5.0fs] sim %d/%d\n", as.numeric(difftime(Sys.time(), t0, units = "secs")), s, n_sims))
}

summarise <- function(sims) {
  out <- lapply(names(sims[[1]]$est), function(m) {
    E <- t(sapply(sims, function(x) x$est[[m]] - x$truth))          # sims x 11 errors
    J <- t(sapply(sims, function(x) x$rej[[m]]))
    all_avg <- rowMeans(E); untr_avg <- rowMeans(E[, -1, drop = FALSE])
    data.frame(panel = panel, model = m,
               all_bias = mean(all_avg), all_mad = mean(abs(all_avg)), all_rmse = sqrt(mean(rowSums(E)^2)),
               tr_bias = mean(E[, 1]), tr_mad = mean(abs(E[, 1])), tr_reject = mean(J[, 1]),
               un_bias = mean(untr_avg), un_mad = mean(abs(untr_avg)), un_reject = mean(J[, -1]))
  })
  do.call(rbind, out)
}
res <- summarise(sims)
saveRDS(list(summary = res, sims = sims, sigma = SIGMA, panel = panel, n_sims = n_sims,
             gsynth_boots = gsynth_boots, seed = seed), out_rds)
print(res, digits = 2, row.names = FALSE)
cat(sprintf("mean N treated %.1f | mean event-day SMB %.2f | gsynth r chosen: %s\n",
            mean(sapply(sims, `[[`, "n_treat")), mean(sapply(sims, `[[`, "smb_event")),
            paste(names(table(sapply(sims, `[[`, "gsynth_r"))), table(sapply(sims, `[[`, "gsynth_r")), sep = "x", collapse = " ")))
cat("Saved", out_rds, "\n")
