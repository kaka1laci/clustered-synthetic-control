# S&P 500 inclusions, 1980 - Sep 1989: GP&L Table 5 (announcement-day effects, 1980-89 row) and
# Table 4 Panel A (betas of included firms).
#
# Event dates: in this period S&P announced additions after the Wednesday close, effective Thursday;
# GP&L infer announcement dates from that rule (Siblis has none before Oct 1989), so CRSP's effective
# dates (crsp.dsp500list.start) give the same events. Two day-0 definitions are computed:
#   effective  = the CRSP start date (229 of 241 are Thursdays; the price reaction is on this day)
#   prior      = the trading day before it (the Wednesday "announcement date" taken literally)
# Only "effective" reproduces Table 5 (DiM 3.25% vs 3.27%); "prior" gives ~0.1%, so Gsynth runs on "effective" only.
# Oct-Dec 1989 additions (5) are dropped: S&P had started pre-announcing, so the rule no longer holds.
#
# Universe on each event date: CRSP share code 10/11, exchange 1-3. Controls: universe firms with a
# return on every day of the window, excluding firms added to the S&P 500 within the window.
# Loadings: days -250..-101 (Table 5 caption). Table 4 betas: days -250..-100.
#
# Usage: Rscript --vanilla replication/src/08_sp500_1980s.R <crsp_dir> <ff_csv> <out_rds> [gsynth=1]

args <- commandArgs(trailingOnly = TRUE)
crsp_dir <- args[1]; ff_csv <- args[2]; out_rds <- args[3]
run_gsynth <- if (length(args) >= 4) as.integer(args[4]) == 1 else TRUE

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
gsynth_lib <- normalizePath(file.path(script_dir, "..", "rlib", "gsynth121"), mustWork = FALSE)
if (dir.exists(gsynth_lib)) .libPaths(c(gsynth_lib, .libPaths()))
if (run_gsynth) suppressMessages(library(gsynth))

t_start <- Sys.time()
log_step <- function(...) {
  cat(sprintf("[%6.0fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = "")
  flush.console()
}

PRE_EST <- c(-250L, -101L)     # loading estimation window
GS_POST <- c(-100L, 20L)       # gsynth "treated" periods (gives the drift path and day 0)

# ---- data ---------------------------------------------------------------------
dsf <- readRDS(file.path(crsp_dir, "crsp_dsf.rds"))
dsi <- readRDS(file.path(crsp_dir, "crsp_dsi.rds"))
nm  <- readRDS(file.path(crsp_dir, "crsp_names.rds"))
spl <- readRDS(file.path(crsp_dir, "crsp_sp500list.rds"))
dsf$date <- as.Date(dsf$date); dsi$date <- as.Date(dsi$date)
nm$namedt <- as.Date(nm$namedt); nm$nameendt <- as.Date(nm$nameendt)
spl$start <- as.Date(spl$start)

cal <- sort(dsi$date)
permnos <- sort(unique(dsf$permno))
R <- matrix(NA_real_, length(permnos), length(cal), dimnames = list(permnos, NULL))
R[cbind(match(dsf$permno, permnos), match(dsf$date, cal))] <- dsf$ret
rm(dsf); invisible(gc())
sp  <- dsi$sprtrn[match(cal, dsi$date)]
vw  <- dsi$vwretd[match(cal, dsi$date)]

ff_lines <- readLines(ff_csv)
ff <- read.csv(text = ff_lines[grep("^,Mkt-RF", ff_lines):length(ff_lines)], strip.white = TRUE)
names(ff) <- c("date", "mkt_rf", "smb", "hml", "rf")
ff <- ff[grepl("^[0-9]{8}$", ff$date), ]
ff$date <- as.Date(as.character(ff$date), "%Y%m%d")
ff <- ff[match(cal, ff$date), ]
FF <- as.matrix(sapply(ff[, c("mkt_rf", "smb", "hml", "rf")], as.numeric)) / 100
log_step(sprintf("returns matrix %d x %d", nrow(R), ncol(R)))

# ---- events -------------------------------------------------------------------
ev <- spl[spl$start >= as.Date("1980-01-01") & spl$start <= as.Date("1989-09-30"), c("permno", "start")]
ev$i_eff <- match(ev$start, cal)
stopifnot(!anyNA(ev$i_eff))
universe_on <- function(d) {
  x <- nm[nm$namedt <= d & nm$nameendt >= d & nm$shrcd %in% c(10, 11) & nm$exchcd %in% 1:3, "permno"]
  intersect(as.character(unique(x)), rownames(R))
}
ev$in_universe <- mapply(function(p, d) as.character(p) %in% universe_on(d), ev$permno, ev$start)
log_step(sprintf("events 1980..1989-09: %d (in common-stock universe: %d), dates: %d",
                 nrow(ev), sum(ev$in_universe), length(unique(ev$start))))
ev <- ev[ev$in_universe, ]

added_within <- function(i0, lo, hi) {
  d_lo <- cal[max(1, i0 + lo)]; d_hi <- cal[min(length(cal), i0 + hi)]
  as.character(spl$permno[spl$start >= d_lo & spl$start <= d_hi])
}

ols_ar <- function(y, X, est, day) {            # abnormal return on `day` with intercept
  ok <- est[!is.na(y[est]) & stats::complete.cases(X[est, , drop = FALSE])]
  if (length(ok) < 100 || is.na(y[day])) return(c(ar = NA, ar_noalpha = NA, beta = rep(NA, ncol(X))))
  b <- coef(lm.fit(cbind(1, X[ok, , drop = FALSE]), y[ok]))
  pred <- sum(c(1, X[day, ]) * b)
  c(ar = y[day] - pred, ar_noalpha = y[day] - (pred - b[1]), beta = b[-1])
}

# ---- per-firm estimators (Diff-in-Means, Market, CAPM, FF3F) -------------------
firm_rows <- list(); beta_rows <- list(); mkt_paths <- list()
for (rule in c("effective", "prior")) {
  for (j in seq_len(nrow(ev))) {
    p <- as.character(ev$permno[j]); i0 <- ev$i_eff[j] - (rule == "prior")
    if (i0 + PRE_EST[1] < 1) next
    est <- (i0 + PRE_EST[1]):(i0 + PRE_EST[2])
    y <- R[p, ]
    if (is.na(y[i0])) next
    univ <- universe_on(cal[i0])
    ctrl <- setdiff(univ, c(p, added_within(i0, PRE_EST[1], 20)))
    ctrl_day0 <- R[ctrl, i0]; ctrl_day0 <- ctrl_day0[!is.na(ctrl_day0)]
    capm_sp <- ols_ar(y, cbind(sp), est, i0)
    capm_vw <- ols_ar(y, cbind(vw), est, i0)
    capm_ff <- ols_ar(y - FF[, "rf"], FF[, "mkt_rf", drop = FALSE], est, i0)
    ff3     <- ols_ar(y - FF[, "rf"], FF[, c("mkt_rf", "smb", "hml")], est, i0)
    firm_rows[[length(firm_rows) + 1]] <- data.frame(
      rule, permno = p, date = cal[i0],
      DiM = y[i0] - mean(ctrl_day0), Market_sp = y[i0] - sp[i0], Market_vw = y[i0] - vw[i0],
      CAPM_sp = capm_sp["ar"], CAPM_sp_noalpha = capm_sp["ar_noalpha"], CAPM_vw = capm_vw["ar"],
      CAPM_ff = capm_ff["ar"], FF3F = ff3["ar"], FF3F_noalpha = ff3["ar_noalpha"])
    if (rule == "effective" && i0 + GS_POST[2] <= length(cal)) {
      w <- (i0 + GS_POST[1]):(i0 + GS_POST[2])
      mkt_paths[[length(mkt_paths) + 1]] <- data.frame(permno = p, k = w - i0, ar = y[w] - sp[w])
    }
    if (rule == "effective") {
      bw <- (i0 - 250):(i0 - 100)
      beta_rows[[length(beta_rows) + 1]] <- data.frame(
        permno = p, capm_sp = ols_ar(y, cbind(sp), bw, i0)[3], capm_vw = ols_ar(y, cbind(vw), bw, i0)[3],
        capm_ff = ols_ar(y - FF[, "rf"], FF[, "mkt_rf", drop = FALSE], bw, i0)[3],
        ff_mkt = ols_ar(y - FF[, "rf"], FF[, 1:3], bw, i0)[3],
        ff_smb = ols_ar(y - FF[, "rf"], FF[, 1:3], bw, i0)[4],
        ff_hml = ols_ar(y - FF[, "rf"], FF[, 1:3], bw, i0)[5])
    }
  }
  log_step(sprintf("per-firm estimators done (%s)", rule))
}
firms <- do.call(rbind, firm_rows); betas <- do.call(rbind, beta_rows)

# ---- Gsynth, one model per event date -------------------------------------------
gs_rows <- list(); gs_paths <- list()
if (run_gsynth) {
  for (rule in "effective") {            # "prior" gives ~0.1% for every per-firm estimator: not the paper's day 0
    dates <- sort(unique(ev$start))
    for (d in seq_along(dates)) {
      i_e <- match(dates[d], cal); i0 <- i_e - (rule == "prior")
      cols <- (i0 + PRE_EST[1]):(i0 + GS_POST[2])
      if (min(cols) < 1 || max(cols) > length(cal)) next
      tr <- intersect(as.character(ev$permno[ev$start == dates[d]]), rownames(R))
      tr <- tr[rowSums(is.na(R[tr, cols, drop = FALSE])) == 0]
      if (!length(tr)) next
      univ <- universe_on(cal[i0])
      ctrl <- setdiff(univ, c(tr, added_within(i0, PRE_EST[1], GS_POST[2])))
      ctrl <- ctrl[rowSums(is.na(R[ctrl, cols, drop = FALSE])) == 0]
      keep <- c(intersect(cols, (i0 + PRE_EST[1]):(i0 + PRE_EST[2])), (i0 + GS_POST[1]):(i0 + GS_POST[2]))
      Y <- R[c(ctrl, tr), keep, drop = FALSE]
      k <- keep - i0
      D <- c(rep(0L, length(ctrl)), rep(1L, length(tr)))
      long <- data.frame(id = rep(seq_len(nrow(Y)), times = ncol(Y)), time = rep(seq_len(ncol(Y)), each = nrow(Y)),
                         y = as.vector(Y), d = as.vector(outer(D, as.integer(k >= GS_POST[1]))))
      fit <- tryCatch({
        invisible(capture.output(g <- suppressMessages(gsynth(y ~ d, data = long, index = c("id", "time"),
                    force = "unit", CV = TRUE, r = c(0, 5), se = FALSE))))
        g }, error = function(e) NULL)
      if (is.null(fit)) next
      # gsynth 1.2.1 with se = FALSE names `att` by the panel's time index (1..ncol(Y)), not by time
      # relative to treatment (checked on simulated data: a jump placed at index t shows up at name t).
      att <- fit$att
      stopifnot(length(att) == ncol(Y), all(as.numeric(names(att)) == seq_len(ncol(Y))))
      post <- k >= GS_POST[1]
      post_k <- k[post]
      path <- setNames(att[post], post_k)
      gs_rows[[length(gs_rows) + 1]] <- data.frame(rule, date = cal[i0], n_treat = length(tr), n_ctrl = length(ctrl),
                                                   r = fit$r.cv, day0 = unname(path["0"]))
      gs_paths[[length(gs_paths) + 1]] <- data.frame(rule, date = cal[i0], n_treat = length(tr), k = post_k, att = unname(path))
      if (d %% 25 == 0) log_step(sprintf("gsynth %s: %d / %d dates", rule, d, length(dates)))
    }
  }
}
gsy <- if (length(gs_rows)) do.call(rbind, gs_rows) else NULL
paths <- if (length(gs_paths)) do.call(rbind, gs_paths) else NULL

# ---- summaries ------------------------------------------------------------------
pct <- function(x) 100 * mean(x, na.rm = TRUE)
summ <- do.call(rbind, lapply(split(firms, firms$rule), function(f) data.frame(
  rule = f$rule[1], n = nrow(f), DiM = pct(f$DiM), Market_sp = pct(f$Market_sp), Market_vw = pct(f$Market_vw),
  CAPM_sp = pct(f$CAPM_sp), CAPM_sp_noalpha = pct(f$CAPM_sp_noalpha), CAPM_vw = pct(f$CAPM_vw), CAPM_ff = pct(f$CAPM_ff),
  FF3F = pct(f$FF3F), FF3F_noalpha = pct(f$FF3F_noalpha),
  Gsynth = if (is.null(gsy)) NA else with(gsy[gsy$rule == f$rule[1], ], 100 * weighted.mean(day0, n_treat)))))
cat("\nTable 5, 1980-89 row (percent). Paper: DiM 3.27 | Market 3.25 | CAPM 3.15 | FF3F 3.05 | Gsynth 3.06\n")
print(summ, digits = 3, row.names = FALSE)
bsum <- sapply(betas[, -1], function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
cat("\nTable 4 Panel A, treated betas (paper mean/sd: CAPM 0.961/0.523, FF Mkt 1.108/0.539, SMB 0.558/0.604, HML -0.148/0.987)\n")
print(round(bsum, 3))
# Pre-inclusion drift: cumulative from day -100 (market-adjusted vs Gsynth); missing days count as 0
mp <- do.call(rbind, mkt_paths); mp$ar[is.na(mp$ar)] <- 0
mkt_car <- cumsum(tapply(mp$ar, mp$k, mean))
drift <- data.frame(k = as.numeric(names(mkt_car)), market_car = unname(mkt_car))
if (!is.null(paths)) {
  w <- tapply(paths$att * paths$n_treat, paths$k, sum) / tapply(paths$n_treat, paths$k, sum)
  drift$gsynth_car <- unname(cumsum(w[order(as.numeric(names(w)))]))[match(drift$k, sort(as.numeric(names(w))))]
}
cat("\nPre-inclusion drift, cumulative from day -100 (percent):\n")
print(transform(drift[drift$k %in% c(-80, -50, -15, -1, 0, 5, 20), ], market_car = round(100 * market_car, 2),
                gsynth_car = round(100 * gsynth_car, 2)), row.names = FALSE)
saveRDS(list(firms = firms, betas = betas, gsynth = gsy, gsynth_paths = paths, summary = summ, beta_summary = bsum,
             drift = drift), out_rds)
log_step(paste("saved", out_rds))
