# Regression check after an R / package upgrade: Panel A point estimates equal the stored results
# (which match the paper in 22 of 24 cells, output/table2_panelA_comparison.md).
# Usage: Rscript --vanilla replication/tests/test_panelA_regression.R <new_csv_dir> <stored_csv_dir>
args <- commandArgs(trailingOnly = TRUE)
for (cn in c("schedule", "personal", "newyork")) {
  f <- paste0("table2_panelA_", cn, ".csv")
  new <- read.csv(file.path(args[1], f)); old <- read.csv(file.path(args[2], f))
  m <- merge(old, new, by = c("connection", "estimator", "variant"), suffixes = c(".old", ".new"))
  stopifnot(nrow(m) == nrow(old))
  diff <- max(abs(m$est.new - m$est.old))
  cat(sprintf("%-8s %2d estimates, max |new - stored| = %.2e\n", cn, nrow(m), diff))
  stopifnot(diff < 1e-6)
}
cat("PASS Panel A point estimates unchanged\n")
