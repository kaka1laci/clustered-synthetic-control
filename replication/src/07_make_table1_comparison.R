# Compare the simulated Table 1 (06_simulation_table1.R output) with the paper's Table 1.
#
# Usage: Rscript --vanilla replication/src/07_make_table1_comparison.R <out_md> <sim_rds>...

args <- commandArgs(trailingOnly = TRUE)
out_md <- args[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
sims  <- lapply(args[-1], readRDS)
ours  <- do.call(rbind, lapply(sims, `[[`, "summary"))
paper <- read.csv(file.path(script_dir, "paper_table1.csv"))

titles <- c(A = "Random Assignment + Random Timing", B = "Assignment Selection + Random Timing",
            C = "Random Assignment + Timing Selection", D = "Assignment Selection + Timing Selection")
labels <- c(SimpleMeans = "Simple Means", CAPM = "CAPM", Correct = "Correct Factor Structure",
            Gsynth = "Gsynth (PCA)")
cols <- c("all_bias", "all_mad", "all_rmse", "tr_bias", "tr_mad", "tr_reject", "un_bias", "un_mad", "un_reject")
pair <- function(o, p) sprintf("%.2f / %.2f", o, p)

info <- sims[[1]]
lines <- c(sprintf("Cells: **ours / paper**. %d simulations per panel; sigma_e = %.3f pp; gsynth bootstrap %d.",
                   info$n_sims, info$sigma, info$gsynth_boots),
           "Bias and MAD in percentage points; Reject = rejection rate of a 5% test (the paper's \"Coverage\").", "")
for (pn in names(titles)) {
  o <- ours[ours$panel == pn, ]
  if (!nrow(o)) next
  p <- paper[paper$panel == pn, ]
  lines <- c(lines, sprintf("### Panel %s: %s", pn, titles[pn]), "",
             "| Model | All: E(Bias) | All: MAD | All: RMSE | Treated: E(Bias) | Treated: MAD | Treated: Reject | Untreated: E(Bias) | Untreated: MAD | Untreated: Reject |",
             "|---|---|---|---|---|---|---|---|---|---|")
  for (m in names(labels)) {
    oo <- o[o$model == m, cols]; pp <- p[p$model == m, cols]
    lines <- c(lines, paste0("| ", labels[m], " | ", paste(mapply(pair, oo, pp), collapse = " | "), " |"))
  }
  lines <- c(lines, "")
}
writeLines(lines, out_md)
cat(lines, sep = "\n")
