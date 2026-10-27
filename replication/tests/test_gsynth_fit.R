# Gsynth wrapper (gsynth_fit.R): reproduces the replicated Table 2 Gsynth exactly, recovers known effects on toy
# data, prints nothing and rejects bad input.
# Usage: Rscript --vanilla replication/tests/test_gsynth_fit.R <src_dir> <output_dir>
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; outd <- args[2]
source(file.path(src, "gsynth_fit.R")); use_gsynth121(src)

# 1. anchors: Table 2's treated firms and controls (every other firm), schedule connection
pa <- readRDS(file.path(outd, "geithner_panel.rds"))
D <- pa$treat$schedule == 1
ga <- gsynth_att(pa$Y[!D, ], pa$Y[D, ], pa$dif < 0, pa$dif >= 0)
ta <- read.csv(file.path(outd, "table2_panelA_schedule.csv")); ta <- ta$est[ta$estimator == "Gsynth"]
pb <- readRDS(file.path(outd, "geithner_panelB_complete.rds")); rows <- pb$sample$schedule
Db <- pb$treat$schedule[rows] == 1; Yb <- pb$Y[rows, ]
gb <- gsynth_att(Yb[!Db, ], Yb[Db, ], pb$dif < 0, pb$dif >= 0)
tb <- read.csv(file.path(outd, "table2_panelB_complete_schedule.csv")); tb <- tb$est[tb$estimator == "Gsynth"]
cat(sprintf("Panel A %.16f (Table 2 %.16f, r %d) | Panel B %.16f (Table 2 %.16f, r %d)\n",
            ga$estimate, ta, ga$r, gb$estimate, tb, gb$r))
stopifnot(length(ta) == 1, length(tb) == 1, abs(ga$estimate - ta) < 1e-10, ga$r == 3,
          abs(gb$estimate - tb) < 1e-10, gb$r == 5,
          identical(dim(ga$gaps), c(sum(D), ncol(pa$Y))),
          abs(mean(ga$gaps[, pa$dif >= 0]) - ga$estimate) < 1e-12, ga$pre_rmse > 0)
cat("PASS anchors: Table 2 Gsynth reproduced (Panel A and Panel B, schedule)\n")

# 2. toy factor model: two factors and unit effects; treated effect 1 (unit 1: 2) after period 40
set.seed(7)
N0 <- 80; N1 <- 5; TT <- 50; T0 <- 40
Fm <- matrix(rnorm(TT * 2), TT); L <- matrix(rnorm((N0 + N1) * 2), N0 + N1)
Y <- L %*% t(Fm) + rnorm(N0 + N1, 0, 1) + matrix(rnorm((N0 + N1) * TT, 0, 0.1), N0 + N1)
pre <- seq_len(TT) <= T0; post <- !pre
eff <- c(2, rep(1, N1 - 1))
Y1 <- Y[N0 + 1:N1, ]; Y1[, post] <- Y1[, post] + eff
g1 <- gsynth_att(Y[1:N0, ], Y1, pre, post)
g0 <- gsynth_att(Y[1:N0, ], Y[N0 + 1:N1, ], pre, post)
cat(sprintf("toy: estimate %.4f (truth %.2f), unit 1 %.4f (truth 2), null %.4f, r %d, pre RMSE %.4f\n",
            g1$estimate, mean(eff), mean(g1$gaps[1, post]), g0$estimate, g1$r, g1$pre_rmse))
stopifnot(abs(g1$estimate - mean(eff)) < 0.1, abs(mean(g1$gaps[1, post]) - 2) < 0.15, abs(g0$estimate) < 0.1,
          g1$r >= 1, g1$pre_rmse < 0.2, isTRUE(all.equal(g1$att, colMeans(g1$gaps))))
cat("PASS toy data: effects recovered, rows in the order of Y1\n")

# 3. silent, and bad input stops with an error
noise <- capture.output(invisible(gsynth_att(Y[1:N0, ], Y1, pre, post)))
bad <- Y1; bad[1, 3] <- NA
fails <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(length(noise) == 0, fails(gsynth_att(Y[1:N0, ], bad, pre, post)),
          fails(gsynth_att(Y[1:N0, ], Y1[, -1], pre, post)), fails(gsynth_att(Y[1:N0, ], Y1, pre, pre)))
cat("PASS no printed output; NA, column mismatch and overlapping pre/post are rejected\n")
