# Gsynth (Xu 2017) as used in GP&L Table 2 (Definition 6) and 02_estimate_table2.R: interactive fixed effects with
# unit fixed effects, the number of factors r chosen by cross-validation over 0..5, gsynth 1.2.1 (the version the
# paper cites; gsynth >= 1.4 picks other r on these panels). Point estimates only (no bootstrap).

# Put replication/rlib/gsynth121 first on the library path and load gsynth 1.2.1.
use_gsynth121 <- function(src_dir) {
  lib <- normalizePath(file.path(src_dir, "..", "rlib", "gsynth121"), mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(gsynth))
  if (getNamespaceVersion("gsynth") != "1.2.1") stop("gsynth 1.2.1 is required, loaded ", getNamespaceVersion("gsynth"))
  invisible(lib)
}

# Y0: donors x periods; Y1: treated units x periods; pre / post: logical over periods (treatment = Y1 rows x post).
# Returns the average effect over the post periods, the treated units' gaps (rows in the order of Y1), the
# per-period effect, the mean per-unit pre-period RMSE, the chosen r and the wall time.
gsynth_att <- function(Y0, Y1, pre, post) {
  stopifnot(is.matrix(Y0), is.matrix(Y1), nrow(Y0) > 0, nrow(Y1) > 0, ncol(Y0) == ncol(Y1),
            length(pre) == ncol(Y0), length(post) == ncol(Y0), !any(pre & post), any(pre), any(post),
            !anyNA(Y0), !anyNA(Y1), "gsynth" %in% loadedNamespaces(), getNamespaceVersion("gsynth") == "1.2.1")
  Y <- rbind(Y0, Y1); D <- rep(0:1, c(nrow(Y0), nrow(Y1)))
  long <- data.frame(id = rep(seq_len(nrow(Y)), times = ncol(Y)), time = rep(seq_len(ncol(Y)), each = nrow(Y)),
                     y = as.vector(Y), d = as.vector(outer(D, as.integer(post))))
  t0 <- proc.time()[["elapsed"]]
  invisible(utils::capture.output(                                    # the cross-validation prints a table
    fit <- suppressMessages(gsynth::gsynth(y ~ d, data = long, index = c("id", "time"), force = "unit",
                                           CV = TRUE, r = c(0, 5), se = FALSE, parallel = FALSE))))
  gaps <- t(fit$Y.tr - fit$Y.ct)                                      # treated ids ascending = the rows of Y1
  list(estimate = fit$att.avg, gaps = gaps, att = colMeans(gaps),
       pre_rmse = mean(sqrt(rowMeans(gaps[, pre, drop = FALSE]^2))), r = fit$r.cv,
       elapsed = proc.time()[["elapsed"]] - t0)
}
