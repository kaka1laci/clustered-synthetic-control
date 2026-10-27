# Pull CRSP daily data from WRDS for a date range.
#   Table 2 Panel B (Geithner, all common stocks as controls): default range 2007-09-01 .. 2008-12-31
#   S&P 500 inclusions, 1980s (Tables 4-5):                    1978-12-01 .. 1990-03-31
#
# Needs a WRDS account with the password in ~/.pgpass:
#   wrds-pgdata.wharton.upenn.edu:9737:wrds:<username>:<password>   (chmod 600)
# CRSP is licensed data: the output stays local under replication/data/raw/.
#
# Usage: Rscript --vanilla replication/src/04_pull_crsp.R <out_dir> [start=2007-09-01] [end=2008-12-31]

suppressMessages({
  library(DBI)
  library(RPostgres)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
START <- if (length(args) >= 2) args[2] else "2007-09-01"   # Panel B: covers day -280 before Nov 21, 2008
END   <- if (length(args) >= 3) args[3] else "2008-12-31"   # Panel B: covers day +10 (Dec 8, 2008)
stopifnot(!is.na(as.Date(START)), !is.na(as.Date(END)))

pgpass <- readLines(path.expand("~/.pgpass"))
entry <- grep("^wrds-pgdata\\.wharton\\.upenn\\.edu:", pgpass, value = TRUE)
if (length(entry) == 0) stop("No WRDS entry in ~/.pgpass")
user <- strsplit(entry[1], ":")[[1]][4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Connecting to WRDS as", user, "(approve the Duo push if one arrives)...\n")
con <- dbConnect(Postgres(), host = "wrds-pgdata.wharton.upenn.edu", port = 9737,
                 dbname = "wrds", user = user, sslmode = "require", bigint = "numeric")
on.exit(dbDisconnect(con))

pull <- function(name, sql) {
  t0 <- Sys.time()
  d <- dbGetQuery(con, sql)
  saveRDS(d, file.path(out_dir, paste0(name, ".rds")))
  cat(sprintf("%-15s %9d rows  (%.0fs)\n", name, nrow(d),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(d)
}

# Name history (share code, exchange code, ticker) for every security alive in the window
pull("crsp_names", sprintf(
  "select permno, namedt, nameendt, shrcd, exchcd, ticker, comnam, ncusip, siccd
     from crsp.dsenames
    where nameendt >= '%s' and namedt <= '%s'", START, END))

# Daily stock returns (legacy SIZ daily stock file); shrout for market capitalisation
pull("crsp_dsf", sprintf(
  "select permno, date, ret, prc, shrout
     from crsp.dsf
    where date between '%s' and '%s'", START, END))

# Daily market indices (S&P 500 and CRSP value-weighted)
pull("crsp_dsi", sprintf(
  "select date, vwretd, ewretd, sprtrn
     from crsp.dsi
    where date between '%s' and '%s'", START, END))

# S&P 500 membership spells (start = date the security enters the index)
pull("crsp_sp500list", "select permno, start, ending from crsp.dsp500list")

cat("Saved to", out_dir, "\n")
