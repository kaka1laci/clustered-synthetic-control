# Which gsynth version/specification reproduces the paper's Gsynth = 0.012 (schedule connections)?
# Usage: R_LIBS=<lib> Rscript --vanilla replication/src/check_gsynth_versions.R <panel_rds> [<lib with gsynth 1.2.1>]
# Results: replication/output/diagnostics/gsynth_v121.txt and gsynth_v140.txt
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 2 && nzchar(args[2])) args[2] else NULL
suppressMessages(library(gsynth, lib.loc = lib))
cat("gsynth version:", as.character(packageVersion("gsynth", lib.loc = lib)), "\n")
p <- readRDS(args[1])
Y <- p$Y; post <- p$dif >= 0
D <- p$treat$schedule
long <- data.frame(id = rep(seq_len(nrow(Y)), times = ncol(Y)),
                   time = rep(seq_len(ncol(Y)), each = nrow(Y)),
                   y = as.vector(Y), d = as.vector(outer(D, as.integer(post))))
specs <- list(
  list(force = "unit",    CV = FALSE, r = 1), list(force = "unit",    CV = FALSE, r = 2),
  list(force = "unit",    CV = FALSE, r = 3), list(force = "two-way", CV = FALSE, r = 0),
  list(force = "two-way", CV = FALSE, r = 1), list(force = "two-way", CV = FALSE, r = 2),
  list(force = "two-way", CV = TRUE,  r = c(0, 5)), list(force = "unit", CV = TRUE, r = c(0, 5))
)
for (s in specs) {
  tm <- system.time(fit <- suppressMessages(suppressWarnings(
    gsynth(y ~ d, data = long, index = c("id", "time"), force = s$force, CV = s$CV, r = s$r, se = FALSE))))
  cat(sprintf("force=%-7s CV=%-5s r_in=%-5s -> r=%d  ATT=%.4f  (%.0fs)\n", s$force, s$CV,
              paste(s$r, collapse = ":"), fit$r.cv, fit$att.avg, tm["elapsed"]))
  flush.console()
}
