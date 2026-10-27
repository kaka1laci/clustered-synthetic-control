# Estimate GP&L Table 2 (Geithner): estimators x 3 connection types, for Panel A or Panel B.
#   Panel A (bank controls, 01_build_geithner_panel.R): all 8 estimators.
#   Panel B (all CRSP common stocks, 05_build_panelB.R): Average, DID, SC, SDID, Gsynth. Market /
#     CAPM / FF3F use treated firms only (Definition 3), so they do not depend on the control
#     group and the paper's Panel B values equal Panel A's.
#
# ATT = average over event days 0..10 of the per-day effect (the synthdid / gsynth convention).
# Where the paper leaves a choice unstated, several variants are computed and saved so the
# matching one can be identified (see the "variant" column).
#
# Usage: Rscript --vanilla replication/src/02_estimate_table2.R \
#          <panel_rds> <out_csv> [connections=all] [placebo_reps=100] [gsynth_boots=1000]
#   connections: "all" or a comma list of schedule,personal,newyork
#   placebo_reps / gsynth_boots = 0 skips that standard error (point estimates only)
#
# Software choices (each checked against the paper's point estimates):
#   DID, SDID : synthdid package
#   SC        : sc_qp.R, the synthdid::sc_estimate objective solved by quadratic programming
#               (synthdid's Frank-Wolfe needs ~15 min per fit here, too slow for placebo SEs)
#   Gsynth    : gsynth 1.2.1, the version the paper cites, from replication/rlib/gsynth121.
#               gsynth >= 1.4 (a wrapper of fect) picks r = 0 or 1 by CV on this panel and does
#               not reproduce the paper; 1.2.1 picks r = 3 and does.

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
gsynth_lib <- normalizePath(file.path(script_dir, "..", "rlib", "gsynth121"), mustWork = FALSE)
if (dir.exists(gsynth_lib)) .libPaths(c(gsynth_lib, .libPaths()))
suppressMessages({
  library(synthdid)
  library(gsynth)
})
source(file.path(script_dir, "sc_qp.R"))
cat("gsynth version:", as.character(packageVersion("gsynth")), "\n")

args <- commandArgs(trailingOnly = TRUE)
panel_rds    <- args[1]
out_csv      <- args[2]
connections  <- if (length(args) >= 3 && args[3] != "all") strsplit(args[3], ",")[[1]] else c("schedule", "personal", "newyork")
PLACEBO_REPS <- if (length(args) >= 4) as.integer(args[4]) else 100L
gsynth_boots <- if (length(args) >= 5) as.integer(args[5]) else 1000L
set.seed(20081121)

t_start <- Sys.time()
log_step <- function(...) {
  cat(sprintf("[%6.1fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")
  flush.console()
}

p    <- readRDS(panel_rds)
Y    <- p$Y
pre  <- p$dif < 0
post <- p$dif >= 0
T0   <- sum(pre)

panel_label <- if (is.null(p$label)) "A" else p$label
# Paper values, Table 2 (identical in v1 and the Aug-2026 version)
paper <- read.csv(file.path(script_dir, "paper_table2.csv"))
paper <- paper[paper$panel == panel_label, setdiff(names(paper), "panel")]
estimator_order <- c("Average", "DID", "Market", "CAPM", "FF3F", "SC", "SDID", "Gsynth")
log_step(sprintf("Panel %s: %d firms x %d days", panel_label, nrow(Y), ncol(Y)))

# ---- abnormal returns from firm-level factor regressions on the pre window ----
abnormal_returns <- function(Y, X, pre, post) {
  # X: T x K factor matrix (no intercept). Returns N x sum(post) matrix of AR.
  Xp <- cbind(1, X[pre, , drop = FALSE])
  coefs <- solve(crossprod(Xp), crossprod(Xp, t(Y[, pre, drop = FALSE])))   # (K+1) x N
  fitted_post <- t(cbind(1, X[post, , drop = FALSE]) %*% coefs)             # N x T1
  Y[, post, drop = FALSE] - fitted_post
}

two_sample <- function(x_treat, x_ctrl) {
  tt <- t.test(as.vector(x_treat), as.vector(x_ctrl), var.equal = TRUE)
  c(est = unname(diff(rev(tt$estimate))), se = unname(tt$stderr))
}
treated_only <- function(x_treat) {
  firm_avg <- rowMeans(x_treat)
  c(est = mean(x_treat), se = sd(firm_avg) / sqrt(length(firm_avg)))
}

ar_sets <- list()
if (!is.null(p$mkt)) {
  rf      <- p$ff$rf
  mkt_sp  <- matrix(p$mkt, ncol = 1)
  ff_full <- as.matrix(p$ff[, c("mkt_rf", "smb", "hml")])
  ff_day0 <- ff_full
  ff_day0[p$dif == 0, "mkt_rf"] <- p$mkt[p$dif == 0] - rf[p$dif == 0]   # 3-4pm S&P in place of full-day Mkt-RF
  Yex <- sweep(Y, 2, rf)                                                  # excess returns for FF3F

  ar_sets <- list(
    Market     = list(ar = Y[, post] - matrix(p$mkt[post], nrow(Y), sum(post), byrow = TRUE),
                      note = "R - R_m (beta = 1, no alpha)"),
    Market_alpha = list(ar = abnormal_returns(Y - matrix(p$mkt, nrow(Y), ncol(Y), byrow = TRUE),
                                              matrix(0, ncol(Y), 0), pre, post),
                      note = "R - alpha_hat - R_m"),
    CAPM       = list(ar = abnormal_returns(Y, mkt_sp, pre, post),
                      note = "R - alpha_hat - beta_hat * S&P500"),
    FF3F       = list(ar = abnormal_returns(Yex, ff_full, pre, post),
                      note = "excess R on FF3, full-day factors on day 0"),
    FF3F_day0sp = list(ar = abnormal_returns(Yex, ff_day0, pre, post),
                      note = "excess R on FF3, day-0 Mkt-RF = 3-4pm S&P")
  )

  # CAPM on the Acemoglu et al. estimation window, days -280..-31 (a firm may miss a few days).
  # GP&L use this window for betas in Section 5.1.2 / Table 3, and it reproduces their CAPM column
  # (0.016 / 0.016 / 0.011); the -256..-31 window gives 0.0154 / 0.0151 / 0.0108.
  capm_window_ar <- function(a, window) {
    est_win <- a$dif >= window[1] & a$dif <= window[2] & !is.na(a$mkt)
    ev <- a$dif >= 0 & a$dif <= 10
    t(sapply(seq_len(nrow(a$Y)), function(i) {
      ok <- est_win & !is.na(a$Y[i, ])
      b <- coef(lm(a$Y[i, ok] ~ a$mkt[ok]))
      a$Y[i, ev] - (b[1] + b[2] * a$mkt[ev])
    }))
  }
  if (!is.null(p$all)) {
    ar_sets$CAPM_w280 <- list(ar = capm_window_ar(p$all, c(-280L, -31L)),
                              note = "R - alpha_hat - beta_hat * S&P500, betas on days -280..-31")
  }
}

results <- list()
add <- function(connection, estimator, variant, est, se, note = "") {
  results[[length(results) + 1]] <<- data.frame(connection, estimator, variant,
                                                est = unname(est), se = unname(se), note)
}

Y_all <- Y
for (conn in connections) {
  # Panel B can use a different sample per connection (treated firms that miss days stay in only
  # for the connection they are treated for); Panel A uses every row.
  rows <- if (is.null(p$sample)) rep(TRUE, nrow(Y_all)) else p$sample[[conn]]
  Y <- Y_all[rows, , drop = FALSE]
  D <- p$treat[[conn]][rows]
  log_step(sprintf("=== %s connections: N1 = %d, N0 = %d ===", conn, sum(D), sum(D == 0)))

  # (1) Average: difference in mean daily returns over days 0..10
  r <- two_sample(Y[D == 1, post], Y[D == 0, post])
  add(conn, "Average", "diff", r["est"], r["se"], "two-sample t-test, stacked firm-days")

  # (3)-(5) abnormal-return estimators: treated-only (Def. 3) and treated-minus-control
  for (nm in names(ar_sets)) {
    est_name <- sub("_.*", "", nm)
    ar <- ar_sets[[nm]]$ar[rows, , drop = FALSE]
    r1 <- treated_only(ar[D == 1, , drop = FALSE])
    add(conn, est_name, paste0(nm, ":treated_only"), r1["est"], r1["se"], ar_sets[[nm]]$note)
    r2 <- two_sample(ar[D == 1, , drop = FALSE], ar[D == 0, , drop = FALSE])
    add(conn, est_name, paste0(nm, ":diff"), r2["est"], r2["se"], ar_sets[[nm]]$note)
  }

  # (2), (6), (7): synthdid package; controls first, treated last
  ord <- c(which(D == 0), which(D == 1))
  Ys  <- Y[ord, , drop = FALSE]
  N0  <- sum(D == 0)
  for (e in list(list("DID", did_estimate), list("SDID", synthdid_estimate))) {
    fit <- e[[2]](Ys, N0, T0)
    log_step(sprintf("%-5s point estimate %.4f", e[[1]], as.numeric(fit)))
    se <- if (PLACEBO_REPS > 0) sqrt(vcov(fit, method = "placebo", replications = PLACEBO_REPS)) else NA_real_
    add(conn, e[[1]], "synthdid", as.numeric(fit), se,
        sprintf("synthdid pkg, placebo SE (%d reps)", PLACEBO_REPS))
    if (PLACEBO_REPS > 0) log_step(sprintf("%-5s placebo SE %.4f", e[[1]], se))
  }
  sc <- sc_estimate_qp(Ys, N0, T0)
  sc_se <- if (PLACEBO_REPS > 0) sc_placebo_se_qp(Ys, N0, T0, PLACEBO_REPS) else NA_real_
  add(conn, "SC", "sc_qp", sc$estimate, sc_se,
      sprintf("synthdid sc objective via quadprog, placebo SE (%d reps)", PLACEBO_REPS))
  log_step(sprintf("SC    %.4f (placebo SE %s)", sc$estimate, format(round(sc_se, 4))))

  # (8) Gsynth: interactive fixed effects with unit FE (Def. 6), r by CV, parametric bootstrap
  long <- data.frame(
    id   = rep(seq_len(nrow(Y)), times = ncol(Y)),
    time = rep(seq_len(ncol(Y)), each = nrow(Y)),
    y    = as.vector(Y),
    d    = as.vector(outer(D, as.integer(post)))
  )
  gs <- suppressMessages(gsynth(y ~ d, data = long, index = c("id", "time"), force = "unit",
               CV = TRUE, r = c(0, 5), se = gsynth_boots > 0, inference = "parametric",
               nboots = max(gsynth_boots, 1L), parallel = FALSE))
  gs_se <- if (gsynth_boots > 0) gs$est.avg[1, "S.E."] else NA_real_
  add(conn, "Gsynth", sprintf("gsynth force=unit r_cv=%d", gs$r.cv),
      gs$att.avg, gs_se, sprintf("gsynth pkg, parametric bootstrap (%d)", gsynth_boots))
  log_step(sprintf("Gsynth %.4f (SE %s), r = %d", gs$att.avg, format(round(gs_se, 4)), gs$r.cv))
}

log_step("all estimators done")
res <- do.call(rbind, results)
res <- merge(res, paper, by = c("connection", "estimator"), all.x = TRUE)
res$gap <- res$est - res$paper_est
res$panel <- panel_label
res <- res[order(match(res$connection, c("schedule", "personal", "newyork")),
                 match(res$estimator, estimator_order), res$variant), ]
write.csv(res, out_csv, row.names = FALSE)

cat("\n")
print(transform(res[, c("connection", "estimator", "variant", "est", "se", "paper_est", "paper_se", "gap")],
                est = round(est, 4), se = round(se, 4), gap = round(gap, 4)), row.names = FALSE)
cat(sprintf("\nSaved %s\n", out_csv))
