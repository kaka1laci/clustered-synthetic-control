# Which gsynth 1.2.1 setting gives the paper's Panel B Gsynth (schedule 0.001)? -> none tried; ATT falls with r (r=15: 0.0051)
# Usage: Rscript --vanilla replication/src/check_panelB_gsynth_specs.R replication/output/geithner_panelB_complete.rds replication/rlib/gsynth121
# Result: replication/output/diagnostics/panelB_gsynth_specs.txt
.libPaths(c(commandArgs(TRUE)[2], .libPaths())); suppressMessages(library(gsynth))
p <- readRDS(commandArgs(TRUE)[1]); post <- p$dif >= 0; D <- p$treat$schedule; Y <- p$Y
long <- data.frame(id = rep(seq_len(nrow(Y)), times = ncol(Y)), time = rep(seq_len(ncol(Y)), each = nrow(Y)),
                   y = as.vector(Y), d = as.vector(outer(D, as.integer(post))))
specs <- list(list("unit", FALSE, 5), list("unit", FALSE, 8), list("unit", FALSE, 10), list("unit", FALSE, 15),
              list("two-way", FALSE, 5), list("unit", TRUE, c(0, 10)), list("two-way", TRUE, c(0, 10)),
              list("unit", TRUE, c(0, 20)))
for (s in specs) {
  tm <- system.time(fit <- suppressMessages(capture.output(g <- gsynth(y ~ d, data = long, index = c("id", "time"),
          force = s[[1]], CV = s[[2]], r = s[[3]], se = FALSE))))
  cat(sprintf("force=%-7s CV=%-5s r=%-5s -> r*=%2d ATT=%.4f (%.0fs)\n", s[[1]], s[[2]], paste(s[[3]], collapse = ":"),
              g$r.cv, g$att.avg, tm["elapsed"])); flush.console()
}
