# Does synthdid's own Frank-Wolfe sc_estimate give the paper's Panel B SC (schedule 0.004)? -> no, 0.0074, same as sc_qp.R
# Usage: Rscript --vanilla replication/src/check_panelB_sc_synthdid.R replication/output/geithner_panelB_complete.rds
# Result: replication/output/diagnostics/panelB_sc_synthdid.txt
suppressMessages(library(synthdid))
p <- readRDS(commandArgs(TRUE)[1]); D <- p$treat$schedule
Ys <- p$Y[c(which(D == 0), which(D == 1)), ]
t0 <- Sys.time()
fit <- sc_estimate(Ys, sum(D == 0), sum(p$dif < 0))
w <- attr(fit, "weights")$omega
cat(sprintf("synthdid sc_estimate (schedule, Panel B): %.4f | donors w>0.01: %d | %.0f s\n",
            as.numeric(fit), sum(w > 0.01), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
