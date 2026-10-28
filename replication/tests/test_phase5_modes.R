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

if (part %in% c("geithner", "all")) {
  source(file.path(src, "sc_qp.R")); source(file.path(src, "clustersc.R"))
  source(file.path(src, "gsynth_fit.R")); use_gsynth121(src)
  pa_rds <- file.path(outd, "geithner_panel.rds"); pb_rds <- file.path(outd, "geithner_panelB_complete.rds")
  # Panel A: the real cohorts and 2 placebo draws with SC, SC on the candidates and Gsynth
  f <- file.path(scr, "geiA_check.rds")
  r <- run("12_jointsc_geithner.R", "run", paste0("panel=", pa_rds), "panel_type=A", paste0("ff=", ff), paste0("out=", f),
           "draws=1-2", "pf_draws=0", "do_real=1", "methods=sc,sc_cand,gsynth")
  if (r$status != 0) { cat(tail(r$out, 20), sep = "\n"); stop("Panel A run failed") }
  stopifnot(any(grepl("497 untreated firms", r$out)))
  x <- readRDS(f); p <- readRDS(pa_rds)
  pre <- p$dif < 0; post <- p$dif >= 0; ctrl <- which(rowSums(p$treat) == 0); sched <- which(p$treat$schedule == 1)
  rl <- x$real
  n_tr <- tapply(rl$n_treated, rl$connection, function(v) v[1])
  stopifnot(identical(sort(unique(rl$connection)), c("newyork", "personal", "schedule")),
            all(n_tr[c("schedule", "personal", "newyork")] == c(15, 9, 38)),
            setequal(unique(paste(x$placebo$method, x$placebo$version)), c("sc avg", "sc_cand avg", "gsynth pf")),
            identical(unique(x$placebo$pool), "all"), all(is.na(x$placebo$fin)), nrow(x$gsynth) == 5)
  sc_direct <- simplex_sc(p$Y[ctrl, ], colMeans(p$Y[sched, ]), pre, post)$estimate
  gs_direct <- gsynth_att(p$Y[ctrl, ], p$Y[sched, , drop = FALSE], pre, post)$estimate
  stopifnot(identical(rl$estimate[rl$connection == "schedule" & rl$method == "sc" & rl$version == "avg"], sc_direct),
            identical(rl$estimate[rl$connection == "schedule" & rl$method == "gsynth"], gs_direct))
  a <- x$placebo[x$placebo$method == "sc", ]; b <- x$placebo[x$placebo$method == "sc_cand", ]
  stopifnot(identical(a$estimate[order(a$draw)], b$estimate[order(b$draw)]))      # n_keep > donors: same donors
  cat("PASS Panel A: 497 donors, cohorts 15 / 9 / 38, SC and Gsynth equal direct fits, sc_cand identical to sc\n")
  r <- run("12_jointsc_geithner.R", "run", paste0("panel=", pa_rds), "panel_type=A", paste0("ff=", ff),
           paste0("out=", file.path(scr, "geiA_bad.rds")), "pools=all,financial", "draws=1-1", "do_real=0", "methods=sc")
  stopifnot(r$status != 0, any(grepl("one pool", r$out)))
  cat("PASS Panel A: asking for a financial pool stops the run\n")
  # Panel B: methods=sc gives the earlier runs' sc/avg on draws 1-2 of both pools (same cohorts)
  f <- file.path(scr, "geiB_check.rds")
  r <- run("12_jointsc_geithner.R", "run", paste0("panel=", pb_rds), paste0("crsp=", crsp), paste0("ff=", ff),
           paste0("out=", f), "pools=all,financial", "draws=1-2", "pf_draws=0", "do_real=0", "methods=sc")
  if (r$status != 0) { cat(tail(r$out, 20), sep = "\n"); stop("Panel B run failed") }
  new <- readRDS(f)$placebo
  for (pool in c("all", "financial")) {
    old <- readRDS(file.path(outd, sprintf("jointsc_geithner_%s_v2.rds", pool)))$placebo
    o <- old[old$method == "sc" & old$version == "avg" & old$draw %in% 1:2, ]; n <- new[new$pool == pool, ]
    stopifnot(nrow(n) == 2, identical(n$estimate[order(n$draw)], o$estimate[order(o$draw)]))
  }
  cat("PASS Panel B methods=sc: same cohorts as the earlier runs (sc/avg identical on draws 1-2)\n")
  r <- run("12_jointsc_geithner.R", "run", paste0("panel=", pb_rds), paste0("crsp=", crsp), paste0("ff=", ff),
           paste0("out=", file.path(scr, "geiB_bad.rds")), "draws=1-1", "do_real=0", "methods=scc")
  stopifnot(r$status != 0, any(grepl("unknown method", r$out)))
  cat("PASS Geithner: an unknown method name stops the run\n")
}

# Review: a valid method list that cannot produce any row (per-firm-only methods with no per-firm draws, or 做法 1
# switched off) must stop instead of saving a malformed rds that combine later rejects with an unrelated error
if (part %in% c("sim", "all")) {
  r <- run("11_jointsc_sim.R", "run", paste0("ff=", ff), paste0("out=", file.path(scr, "sim_empty.rds")),
           "scenarios=weak_blur", "reps=1-1", "pf_reps=0", "methods=jsc2")
  stopifnot(r$status != 0, any(grepl("no method selected", r$out)))
  cat("PASS sim: a method list that produces no rows stops the run\n")
}
if (part %in% c("geithner", "all")) {
  pa_rds <- file.path(outd, "geithner_panel.rds")
  r <- run("12_jointsc_geithner.R", "run", paste0("panel=", pa_rds), "panel_type=A", paste0("ff=", ff),
           paste0("out=", file.path(scr, "geiA_empty.rds")), "draws=1-1", "pf_draws=0", "do_real=0", "methods=csc_ff3_cand")
  r2 <- run("12_jointsc_geithner.R", "run", paste0("panel=", pa_rds), "panel_type=A", paste0("ff=", ff),
            paste0("out=", file.path(scr, "geiA_empty2.rds")), "draws=1-1", "do_real=0", "methods=jsc1_ff3", "jsc1=0")
  stopifnot(r$status != 0, any(grepl("no method selected", r$out)), r2$status != 0, any(grepl("no method selected", r2$out)))
  cat("PASS Geithner: a method list that produces no rows stops the run\n")
}
