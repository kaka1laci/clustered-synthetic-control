# Where does the S&P 500 inclusion jump sit relative to the CRSP effective date (crsp.dsp500list.start)?
# Decides whether effective dates can stand in for the Siblis announcement dates used by GP&L. Needs WRDS (~/.pgpass).
# Usage: Rscript --vanilla replication/src/check_sp500_effective_dates.R replication/data/raw/crsp/sp500_effective_windows.rds
# Result: replication/output/diagnostics/sp500_effective_dates.txt
suppressMessages({library(DBI); library(RPostgres)})
out_rds <- commandArgs(TRUE)[1]
u <- strsplit(grep("^wrds-pgdata", readLines(path.expand("~/.pgpass")), value = TRUE)[1], ":")[[1]][4]
con <- dbConnect(Postgres(), host = "wrds-pgdata.wharton.upenn.edu", port = 9737, dbname = "wrds", user = u, sslmode = "require", bigint = "numeric")
ev <- dbGetQuery(con, "
  select l.permno, l.start, d.date, d.ret, i.vwretd, i.sprtrn
    from crsp.dsp500list l
    join crsp.dsf d on d.permno = l.permno and d.date between l.start - 60 and l.start + 30
    join crsp.dsi i on i.date = d.date
   where l.start between '1976-01-01' and '2020-12-31'")
dbDisconnect(con)
saveRDS(ev, out_rds)
ev$start <- as.Date(ev$start); ev$date <- as.Date(ev$date)
cal <- sort(unique(ev$date))
k <- match(ev$date, cal) - match(ev$start, cal)       # trading-day offset; NA if start not a trading day in sample
ev$k <- k; ev$ar <- ev$ret - ev$sprtrn
ev$era <- cut(as.integer(format(ev$start, "%Y")), c(1975, 1989, 1999, 2009, 2020), labels = c("1976-89", "1990-99", "2000-09", "2010-20"))
ev$era[ev$start >= as.Date("1989-10-01") & ev$era == "1976-89"] <- NA   # S&P began pre-announcing in Oct 1989
cat("events:", length(unique(paste(ev$permno, ev$start))), "\n\n")
tab <- with(ev[!is.na(ev$k) & ev$k >= -10 & ev$k <= 3 & !is.na(ev$era), ], tapply(ar, list(k, era), mean, na.rm = TRUE))
cat("Mean market-adjusted return (pp) by trading day relative to CRSP effective date:\n")
print(round(100 * tab, 2))
win <- ev[!is.na(ev$k) & ev$k >= -10 & ev$k <= 0 & !is.na(ev$era), ]
share <- sapply(split(win, win$era), function(d) { m <- tapply(d$ar, d$k, mean, na.rm = TRUE); c(total_m10_0 = sum(m), largest_day = as.numeric(names(which.max(m))), share_largest = max(m) / sum(m)) })
cat("\nCumulative AR days -10..0 (pp), day with the largest AR, and its share of the cumulative:\n"); print(round(share * c(100, 1, 1), 2))
