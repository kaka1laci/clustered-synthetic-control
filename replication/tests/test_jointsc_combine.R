# Review C1. combine must not read float noise as a significant paired difference: when two methods give the same
# estimates up to rounding (做法 2 staying at its ClusterSC start), the paired row reports mean 0 and no t / p.
# Usage: Rscript --vanilla replication/tests/test_jointsc_combine.R <src_dir> <scratch_dir>
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(3)
noisy <- function(x) x * (1 + 1e-15 * sample(c(-1, 0, 1), length(x), replace = TRUE))   # equal up to rounding

# simulation: csc_ff3/pf and jsc2_ff3/pf equal up to rounding, jsc1_ff3/pf genuinely different
e_csc <- rnorm(10, 0.2, 0.1); a_km <- runif(10, 0.3, 0.6)
row <- function(method, version, err0, ari = NA) data.frame(scenario = "weak_clear", rep = 1:10, method = method,
  version = version, err0 = err0, post_bias = 0, w_g1 = 0.8, ari = ari, runtime = 1, mipgap = NA_real_, n_nonopt = 0L,
  converged = NA, n_iter = NA_integer_)
met <- rbind(row("csc_ff3", "pf", e_csc), row("jsc2_ff3", "pf", noisy(e_csc), noisy(a_km)),
             row("jsc1_ff3", "pf", e_csc + rnorm(10, -0.1, 0.02)), row("csc_ff3", "avg", rnorm(10), a_km))
sim <- list(metrics = met, choices = data.frame(scenario = "weak_clear", rep = 1, method = "jsc1_ff3", version = "pf",
                                                 target = 1, m = 50, c = 1),
            dgp = data.frame(scenario = "weak_clear", rep = 1:10, cor_mkt = 0.9, cor_smb = 0.5))
saveRDS(sim, file.path(out, "fake_sim.rds"))
stopifnot(system2("Rscript", c("--vanilla", file.path(src, "11_jointsc_sim.R"), "combine", file.path(out, "fake_sim"),
                               file.path(out, "fake_sim.rds")), stdout = FALSE) == 0)
pr <- read.csv(file.path(out, "fake_sim_paired.csv")); ar <- read.csv(file.path(out, "fake_sim_ari.csv"))
same <- pr[pr$comparison == "jsc2_ff3/pf vs csc_ff3/pf", ]; diff <- pr[pr$comparison == "jsc1_ff3/pf vs csc_ff3/pf", ]
print(rbind(same, diff), digits = 3)
stopifnot(nrow(same) == 1, same$mean_diff == 0, is.na(same$t), is.na(same$p),
          nrow(diff) == 1, is.finite(diff$p), diff$p < 0.05,
          ar$mean_diff[ar$method == "jsc2_ff3/pf"] == 0, is.na(ar$p[ar$method == "jsc2_ff3/pf"]))
cat("PASS simulation combine: identical-up-to-rounding methods get no t / p\n")

# Geithner: the same on decimal-return estimates, both paired measures
e <- rnorm(10, -0.003, 0.004)
prow <- function(method, estimate) data.frame(pool = "financial", draw = 1:10, method = method, version = "pf",
  estimate = estimate, pre_rmse = 0.02, fin = 0.4, runtime = 1, mipgap = NA_real_, n_nonopt = 0L, converged = NA,
  n_iter = NA_integer_)
gei <- list(placebo = rbind(prow("csc_ff3", e), prow("jsc2_ff3", noisy(e)), prow("jsc1_ff3", e + rnorm(10, 0.002, 0.0005))),
            real = NULL, choices = data.frame(pool = "financial", draw = 1, method = "jsc1_ff3", version = "pf",
                                              target = 1, m = 50, c = 1))
saveRDS(gei, file.path(out, "fake_gei.rds"))
stopifnot(system2("Rscript", c("--vanilla", file.path(src, "12_jointsc_geithner.R"), "combine", file.path(out, "fake_gei"),
                               file.path(out, "fake_gei.rds")), stdout = FALSE) == 0)
pg <- read.csv(file.path(out, "fake_gei_paired.csv"))
sg <- pg[pg$comparison == "jsc2_ff3/pf vs csc_ff3/pf", ]; dg <- pg[pg$comparison == "jsc1_ff3/pf vs csc_ff3/pf", ]
print(rbind(sg, dg), digits = 3)
stopifnot(nrow(sg) == 2, all(sg$mean_diff == 0), all(is.na(sg$t)), all(is.na(sg$p)), nrow(dg) == 2, all(is.finite(dg$p)))
cat("PASS Geithner combine: identical-up-to-rounding methods get no t / p\n")

# Review M4. combine must refuse duplicate rows (e.g. the same method from two runs passed together)
fails_with <- function(script, prefix, rds, pattern) {
  res <- suppressWarnings(system2("Rscript", c("--vanilla", file.path(src, script), "combine", file.path(out, prefix), rds),
                                  stdout = TRUE, stderr = TRUE))
  !is.null(attr(res, "status")) && attr(res, "status") != 0 && any(grepl(pattern, res))
}
dup_sim <- sim; dup_sim$metrics <- rbind(met, met[1, ])
saveRDS(dup_sim, file.path(out, "dup_sim.rds"))
dup_gei <- gei; dup_gei$placebo <- rbind(gei$placebo, gei$placebo[1, ])
saveRDS(dup_gei, file.path(out, "dup_gei.rds"))
dup_real <- gei; dup_real$real <- cbind(connection = "schedule", n_treated = 12, gei$placebo[c(1, 1), -(1:2)])
saveRDS(dup_real, file.path(out, "dup_real.rds"))
stopifnot(fails_with("11_jointsc_sim.R", "dup_sim", file.path(out, "dup_sim.rds"), "duplicate"),
          fails_with("12_jointsc_geithner.R", "dup_gei", file.path(out, "dup_gei.rds"), "duplicate"),
          fails_with("12_jointsc_geithner.R", "dup_real", file.path(out, "dup_real.rds"), "duplicate"))
cat("PASS combine stops on duplicate rows (simulation, Geithner placebo, Geithner real)\n")
