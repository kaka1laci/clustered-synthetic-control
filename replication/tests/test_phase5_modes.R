# Phase 5 run modes: methods= runs only the listed methods on the same data / draws as the earlier full runs,
# Panel A support, Gsynth rows. Needs the local (gitignored) outputs of the 方案 B runs:
# jointsc_sim_p2.rds, jointsc_geithner_{all,financial}_v2.rds, and the Panel A / B panels.
# Usage: Rscript --vanilla replication/tests/test_phase5_modes.R <src_dir> <output_dir> <ff_csv> <crsp_dir> <scratch_dir> [sim|geithner|all]
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; outd <- args[2]; ff <- args[3]; crsp <- args[4]; scr <- args[5]
part <- if (length(args) >= 6) args[6] else "all"
dir.create(scr, showWarnings = FALSE, recursive = TRUE)
run <- function(script, ...) {
  res <- suppressWarnings(system2("Rscript", c("--vanilla", file.path(src, script), ...), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(res, "status"))) 0L else attr(res, "status"), out = res)
}

if (part %in% c("sim", "all")) {
  f <- file.path(scr, "sim_check.rds")
  r <- run("11_jointsc_sim.R", "run", paste0("ff=", ff), paste0("out=", f), "scenarios=weak_blur", "reps=1-2",
           "pf_reps=0", "methods=sc,gsynth")
  if (r$status != 0) { cat(tail(r$out, 20), sep = "\n"); stop("sim run failed") }
  x <- readRDS(f); m <- x$metrics
  old <- readRDS(file.path(outd, "jointsc_sim_p2.rds"))$metrics
  old <- old[old$scenario == "weak_blur" & old$rep %in% 1:2 & old$method == "sc" & old$version == "avg", ]
  new <- m[m$method == "sc", ]
  stopifnot(nrow(m) == 4, setequal(paste(m$method, m$version), c("sc avg", "gsynth pf")),
            identical(new$err0[order(new$rep)], old$err0[order(old$rep)]),
            all(is.finite(m$err0[m$method == "gsynth"])), all(is.na(m$w_g1[m$method == "gsynth"])),
            nrow(x$gsynth) == 2, all(x$gsynth$r %in% 0:5), is.null(x$choices))
  cat("PASS sim methods=sc,gsynth: same data as the full run (sc/avg identical), Gsynth rows written\n")
  r <- run("11_jointsc_sim.R", "run", paste0("ff=", ff), paste0("out=", file.path(scr, "sim_bad.rds")),
           "scenarios=weak_blur", "reps=1-1", "methods=sc,gsynthh")
  stopifnot(r$status != 0, any(grepl("unknown method", r$out)))
  cat("PASS sim: an unknown method name stops the run\n")
}
