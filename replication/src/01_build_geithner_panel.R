# Build the Geithner event-study panel used in Goldsmith-Pinkham & Lyu, Table 2 Panel A.
#
# Follows the sample construction of Acemoglu et al. (2016) replication files
# ("Geithner event study.do", "create correlations.do", "Table2.do"), then keeps
# the balanced window used by GP&L: pre days -256..-31 and event days 0..10.
#
# Usage: Rscript --vanilla replication/src/01_build_geithner_panel.R <kermani_stata_dir> <ff_csv> <out_rds>

suppressMessages({
  library(haven)
})

args <- commandArgs(trailingOnly = TRUE)
stata_dir <- args[1]
ff_csv    <- args[2]
out_rds   <- args[3]

EVENT_DATE    <- as.Date("2008-11-21")
MKT_RET_DAY0  <- 0.06876978          # S&P 500 return 3pm-4pm on Nov 21, 2008 (Acemoglu et al.)
PRE_WINDOW    <- c(-256L, -31L)      # GP&L balanced pre-window
EVENT_WINDOW  <- c(0L, 10L)

stata_date <- function(x) as.Date(x, origin = "1960-01-01")
holidays <- stata_date(c(17351, 17412, 17492, 17525, 17532, 17552, 17580, 17612, 17678,
                         17717, 17776, 17863, 17891, 17898, 17916, 17944, 17997))

# ---- returns (long) --------------------------------------------------------
ri <- read_dta(file.path(stata_dir, "Return Data.dta"))
ri_cols <- grep("^ri", names(ri), value = TRUE)
dates <- as.Date(sub("^ri", "", ri_cols), format = "%d%b%y")
stopifnot(!anyNA(dates))

keep_dates <- !(dates %in% holidays)
dates  <- dates[keep_dates]
ri_mat <- as.matrix(ri[, ri_cols[keep_dates]])
ord    <- order(dates)
dates  <- dates[ord]
ri_mat <- ri_mat[, ord, drop = FALSE]

ret_mat <- cbind(NA, ri_mat[, -1] / ri_mat[, -ncol(ri_mat)] - 1)
colnames(ret_mat) <- format(dates)
rownames(ret_mat) <- ri$index

dif <- seq_along(dates) - which(dates == EVENT_DATE)
stopifnot(sum(dif == 0) == 1)

market_index <- 679
mkt <- ret_mat[as.character(market_index), ]
firm_rows <- ri$index != market_index
ret <- ret_mat[firm_rows, , drop = FALSE]
tickers <- ri$ticker[firm_rows]

# ---- Acemoglu et al. sample filters ---------------------------------------
est_obs   <- rowSums(!is.na(ret[, dif >= -280 & dif < -30, drop = FALSE]))
event_obs <- rowSums(!is.na(ret[, dif >= 0 & dif <= 10, drop = FALSE]))
keep_ace  <- est_obs >= 225 & event_obs >= 11

# day-0 return is the 3pm-4pm return
nov21 <- read_dta(file.path(stata_dir, "Return Nov. 21 Data.dta"))
ret0  <- nov21$ret3to4[match(tickers, nov21$ticker)]
ret[, dif == 0] <- ret0
mkt[dif == 0]   <- MKT_RET_DAY0
keep_taq <- !is.na(ret0)

# connections and characteristics
conn  <- read_dta(file.path(stata_dir, "Connections Data.dta"))
sched <- read_dta(file.path(stata_dir, "Schedule Data.dta"))
fin   <- read_dta(file.path(stata_dir, "Financial Data.dta"))
personal_num <- conn$geithner2num[match(tickers, conn$ticker)]
sched_num    <- sched$geithnerschednum[match(tickers, sched$ticker)]
city         <- fin$addresscity_ds[match(tickers, fin$ticker)]
fullsample   <- fin$fullsample[match(tickers, fin$ticker)]

# correlation with Citigroup over the 250-day window (Stata dates 17454..17814)
corr_window <- dates >= stata_date(17454) & dates <= stata_date(17814)
citi <- ret[which(tickers == "C"), corr_window]
citicorr <- apply(ret[, corr_window, drop = FALSE], 1, function(r) {
  ok <- !is.na(r) & !is.na(citi)
  if (sum(ok) < 3) NA_real_ else cor(r[ok], citi[ok])
})
in_full <- !is.na(fullsample) & fullsample == 1
# Stata egen pctile (default definition) == R quantile type 2
citicorr_p90 <- quantile(tapply(citicorr[in_full], tickers[in_full], mean),
                         0.90, type = 2, na.rm = TRUE)
# The cutoff is computed over fullsample firms (as in "create correlations.do"), but the
# cutoff is applied to every firm without requiring fullsample == 1. This reproduces the
# N = 545 firms in GP&L Table 2 Panel A (requiring fullsample too gives 525).
keep_citi <- !is.na(citicorr) & citicorr <= citicorr_p90

keep <- keep_ace & keep_taq & keep_citi
cat(sprintf("Firms: %d total | Acemoglu window filter %d | + TAQ day-0 %d | + Citi-corr <= p90 %d\n",
            length(tickers), sum(keep_ace), sum(keep_ace & keep_taq), sum(keep)))

# ---- GP&L balanced window ---------------------------------------------------
in_window <- (dif >= PRE_WINDOW[1] & dif <= PRE_WINDOW[2]) |
             (dif >= EVENT_WINDOW[1] & dif <= EVENT_WINDOW[2])
Y <- ret[keep, in_window, drop = FALSE]
complete <- rowSums(is.na(Y)) == 0
cat(sprintf("Window: %d days (%d pre, %d post) | firms complete in window: %d of %d\n",
            sum(in_window), sum(in_window & dif < 0), sum(in_window & dif >= 0),
            sum(complete), nrow(Y)))

# ---- Fama-French daily factors ---------------------------------------------
ff_lines <- readLines(ff_csv)
start <- grep("^,Mkt-RF", ff_lines)
ff <- read.csv(text = ff_lines[start:length(ff_lines)], header = TRUE, strip.white = TRUE)
names(ff) <- c("date", "mkt_rf", "smb", "hml", "rf")
ff <- ff[grepl("^[0-9]{8}$", ff$date), ]
ff$date <- as.Date(as.character(ff$date), format = "%Y%m%d")
ff[, -1] <- lapply(ff[, -1], function(x) as.numeric(x) / 100)

win_dates <- dates[in_window]
ffw <- ff[match(win_dates, ff$date), ]
stopifnot(!anyNA(ffw$date))

panel <- list(
  Y          = Y,                                  # firms x days, daily simple returns
  dates      = win_dates,
  dif        = dif[in_window],
  ticker     = tickers[keep],
  complete   = complete,
  treat      = data.frame(
    schedule = as.integer(!is.na(sched_num[keep]) & sched_num[keep] > 0),
    personal = as.integer(!is.na(personal_num[keep]) & personal_num[keep] > 0),
    newyork  = as.integer(!is.na(city[keep]) & city[keep] == "New York")
  ),
  mkt        = mkt[in_window],                     # S&P 500 return (Datastream row 679), day 0 = 3-4pm
  ff         = ffw,                                # full-day FF factors (no intraday version exists)
  all        = list(Y = ret[keep, , drop = FALSE], # every trading day, NA where missing (robustness checks)
                    dif = dif, dates = dates, mkt = mkt)
)
saveRDS(panel, out_rds)

cat("Treated counts:\n"); print(colSums(panel$treat[complete, ]))
cat(sprintf("Saved %s\n", out_rds))
