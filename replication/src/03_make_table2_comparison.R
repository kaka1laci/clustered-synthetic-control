# Combine per-connection results and write a replication-vs-paper table (markdown + csv).
#
# Usage: Rscript --vanilla replication/src/03_make_table2_comparison.R <out_md> <in_csv>...

args <- commandArgs(trailingOnly = TRUE)
out_md <- args[1]
res <- do.call(rbind, lapply(args[-1], read.csv))

# The variant that matches the paper's definition of each estimator:
#   abnormal-return estimators = mean abnormal return of treated firms only (Definition 3),
#   CAPM betas on days -280..-31, FF3F with the published full-day factors on day 0.
chosen <- c(Average = "^diff$", DID = "^synthdid$", Market = "^Market:treated_only$",
            CAPM = "^CAPM_w280:treated_only$", FF3F = "^FF3F:treated_only$", SC = "^sc_qp$",
            SDID = "^synthdid$", Gsynth = "^gsynth")
keep <- mapply(function(e, v) res$estimator == e & grepl(v, res$variant), names(chosen), chosen)
tab <- res[rowSums(keep) > 0, ]
tab$match <- round(tab$est, 3) == round(tab$paper_est, 3)

conns <- c(schedule = "Schedule connections", personal = "Personal connections",
           newyork = "New York connections")
fmt <- function(est, se) sprintf("%.4f (%s)", est, ifelse(is.na(se), "–", sprintf("%.4f", se)))
# Significance stars from |est / se| (two-sided normal: 10% *, 5% **, 1% ***). The paper's stars are
# recomputed from its rounded numbers, so a borderline case can differ from the printed table.
stars <- function(est, se) {
  if (is.na(se) || se <= 0) return("–")
  z <- abs(est / se)
  if (z > 2.576) "\\*\\*\\*" else if (z > 1.96) "\\*\\*" else if (z > 1.645) "\\*" else "n.s."
}
tab$sig <- mapply(stars, tab$est, tab$se)
tab$paper_sig <- mapply(stars, tab$paper_est, tab$paper_se)
lines <- c(
  "| Connection | Estimator | Ours: est (SE) | Paper: est (SE) | Gap | Same at 3 dp | Sig: ours / paper |",
  "|---|---|---|---|---|---|---|")
for (cn in names(conns)) for (e in names(chosen)) {
  r <- tab[tab$connection == cn & tab$estimator == e, ]
  if (nrow(r) != 1) next
  lines <- c(lines, sprintf("| %s | %s | %s | %.3f (%.3f) | %+.4f | %s | %s / %s |", conns[cn], e,
                            fmt(r$est, r$se), r$paper_est, r$paper_se, r$gap,
                            if (r$match) "yes" else "**no**", r$sig, r$paper_sig))
}
writeLines(lines, out_md)
write.csv(tab, sub("\\.md$", ".csv", out_md), row.names = FALSE)
cat(lines, sep = "\n")
cat(sprintf("\n%d of %d point estimates equal the paper at 3 decimals\n", sum(tab$match), nrow(tab)))
