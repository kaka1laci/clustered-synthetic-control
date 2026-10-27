# Build the Panel B panel (Geithner, all NYSE/AMEX/NASDAQ common stocks as controls) from CRSP.
#
# Sample (reproduces the paper's 4,078 firms x 237 days):
#   CRSP securities with share code 10/11 and exchange code 1-3 on Nov 21, 2008, with a return on
#   every day of the window (days -256..-31 and 0..10). Returns are CRSP close-to-close for every
#   firm, treated ones included (no 3pm-4pm day-0 return in Panel B).
# Treated: the Acemoglu et al. connected firms, matched to PERMNO by ticker on the event date and
#   kept only if they are in the common-stock universe (drops AB, BX, IEP, LAZ, OZM).
# treated_rule:
#   "complete"      treated firms also need all 237 returns (main; matches N = 4,078)
#   "allow_missing" treated firms may miss days (filled with 0); CARV and CMSB each miss 2008-03-10.
#                   Reproduces the paper's personal-connection Average/DID but not New York's.
#
# Usage: Rscript --vanilla replication/src/05_build_panelB.R <crsp_dir> <panelA_rds> <out_rds> [treated_rule]

args <- commandArgs(trailingOnly = TRUE)
crsp_dir     <- args[1]
panelA_rds   <- args[2]
out_rds      <- args[3]
treated_rule <- if (length(args) >= 4) args[4] else "complete"
stopifnot(treated_rule %in% c("complete", "allow_missing"))

EVENT <- as.Date("2008-11-21")
crsp_names <- readRDS(file.path(crsp_dir, "crsp_names.rds"))
dsf <- readRDS(file.path(crsp_dir, "crsp_dsf.rds"))
dsi <- readRDS(file.path(crsp_dir, "crsp_dsi.rds"))
pa  <- readRDS(panelA_rds)

cal <- sort(as.Date(dsi$date))
dif <- seq_along(cal) - which(cal == EVENT)
in_window <- (dif >= -256 & dif <= -31) | (dif >= 0 & dif <= 10)
win <- cal[in_window]

dsf$date <- as.Date(dsf$date)
permnos <- sort(unique(dsf$permno))
k <- match(dsf$date, win)
R <- matrix(NA_real_, length(permnos), length(win), dimnames = list(permnos, format(win)))
R[cbind(match(dsf$permno[!is.na(k)], permnos), k[!is.na(k)])] <- dsf$ret[!is.na(k)]

crsp_names$namedt   <- as.Date(crsp_names$namedt)
crsp_names$nameendt <- as.Date(crsp_names$nameendt)
on_event <- crsp_names[crsp_names$namedt <= EVENT & crsp_names$nameendt >= EVENT, ]
universe <- as.character(unique(on_event$permno[on_event$shrcd %in% c(10, 11) & on_event$exchcd %in% 1:3]))
universe <- intersect(universe, rownames(R))
complete <- universe[rowSums(is.na(R[universe, , drop = FALSE])) == 0]

# Acemoglu connected firms -> PERMNO (ticker on the event date)
acem_permno <- as.character(on_event$permno[match(pa$ticker, on_event$ticker)])
stopifnot(!anyNA(acem_permno[rowSums(pa$treat) > 0]))
treated <- lapply(pa$treat, function(d) {
  tp <- intersect(unique(acem_permno[d == 1]), universe)
  if (treated_rule == "complete") intersect(tp, complete) else tp
})

rows <- union(complete, unlist(treated))
Y <- R[rows, , drop = FALSE]
filled <- rowSums(is.na(Y))
Y[is.na(Y)] <- 0

treat <- as.data.frame(lapply(treated, function(tp) as.integer(rows %in% tp)))
# Controls for each connection = complete firms not treated for it; treated firms that miss days
# never serve as controls.
sample <- lapply(treated, function(tp) rows %in% union(complete, tp))

cat(sprintf("Universe on event date: %d | complete over 237 days: %d | rows: %d | rule: %s\n",
            length(universe), length(complete), length(rows), treated_rule))
for (cn in names(treated)) {
  cat(sprintf("  %-8s treated %2d | sample N %d | treated with filled days: %s\n", cn, length(treated[[cn]]),
              sum(sample[[cn]]), paste(names(filled)[filled > 0 & rows %in% treated[[cn]]], collapse = " ")))
}

ticker_map <- on_event[match(as.integer(rows), on_event$permno), c("permno", "ticker", "comnam")]
saveRDS(list(label = "B", Y = Y, dates = win, dif = dif[in_window], permno = rows,
             ticker = ticker_map$ticker, treat = treat, sample = sample, filled = filled,
             treated_rule = treated_rule), out_rds)
cat("Saved", out_rds, "\n")
