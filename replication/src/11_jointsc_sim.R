# 方案 B simulation: when does joint group selection beat SC and ClusterSC? (true effect known)
#
# DGP (Table 1's two-factor model plus a group-specific factor), percentage points:
#   r_it = rf_t + b_mkt_i (mkt_t - rf_t) + b_smb_i smb_t + delta_i g_{k(i),t} + e_it
#   400 units in 3 groups (134 / 133 / 133); loading centres (b_mkt, b_smb) = (0.8, 0.6) / (1.0, 1.2) / (1.3, 0.8),
#   within-group sd 0.1 ("clear") or 0.3 ("blur"); g_kt ~ N(0, (kappa * sd_smb)^2), independent across groups,
#   kappa in {0, 0.5, 1} ("none" / "weak" / "strong"), sd_smb = sd of daily SMB over 1926-07..2022-12;
#   delta_i ~ N(1, 0.3^2); e_it ~ N(0, 1.209^2) (Table 1 calibration). Factors: a random block of 250
#   consecutive FF days (no timing selection). T0 = 239, event day t = 240, days 241..250 after; effect 3 on the
#   event day for the 12 treated units, all from group 1. Donors = the other 388 units.
# Methods (simplex weights throughout), cohort-average target (avg) and per-firm targets (pf):
#   sc        SC on all donors                                   avg, pf
#   csc_ff3   ClusterSC: k-means (k = 3) on FF loadings, SC in the target's cluster      avg, pf
#   csc_svd   ClusterSC: k-means (k = 3) on SVD features (r = 10)                        avg
#   jsc1_ff3  做法 1 (jsc1_tune) with FF-loading distances, ClusterSC start              avg, pf
#   jsc1_svd  做法 1 with SVD distances                                                   avg
#   jsc2_ff3  做法 2 (K = 3) on the 12 treated units, FF loadings, k-means start          pf
#   jsc2_ff3_c<c>_gap<g>  optional extra 做法 2 arm with its own gamma factor and MIP gap (jsc2_alt=c:g)   pf
#   gsynth    Gsynth (gsynth 1.2.1, unit FE, r by CV over 0..5; gsynth_fit.R) on the 12 treated units and the 388
#             donors, every rep; counted as per-firm since it averages per-unit effects        pf
# Metrics per (scenario, rep, method): event-day error (estimate - 3), mean gap on days 241..250 (truth 0),
#   weight on group-1 donors, ARI of the donor partition (k-means for csc_*, final partition for jsc2), wall time.
#
# Usage:
#   run:     Rscript --vanilla replication/src/11_jointsc_sim.R run ff=<ff_csv> out=<rds> [scenarios=all|a,b]
#              [reps=1-30] [pf_reps=30] [grid=full|small] [tl1=60] [tl2=300] [gap=0.01] [jsc2_c=1]
#              [jsc2_alt=c:gap] [methods=all] [seed=2026]
#   combine: Rscript --vanilla replication/src/11_jointsc_sim.R combine <out_prefix> <rds>...
#   pf_reps: per-firm methods run on reps <= pf_reps only; grid=small: m = 50, c in {0.1, 1, 10}.
#   jsc2_alt: at the main c = 1 做法 2 stays at its k-means start (pilot, 10-26); this arm (e.g. 0.1:0.001) runs
#     next to the main arm and gets the main arm's comparisons in combine.
#   methods: run only the listed methods (sc, csc_ff3, csc_svd, jsc1_ff3, jsc1_svd, jsc2 = every 做法 2 arm, gsynth);
#     the data depend only on (scenario, rep), so methods=gsynth adds Gsynth to an earlier run (combine them).

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
kv <- function(a) { s <- strsplit(a, "=", fixed = TRUE); setNames(lapply(s, `[`, 2), sapply(s, `[`, 1)) }
mx <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
paired_stats <- function(dd) {
  se <- sd(dd) / sqrt(length(dd))
  data.frame(n = length(dd), mean_diff = mean(dd), se = se,
             t = if (se > 0) mean(dd) / se else NA_real_, p = if (se > 0) t.test(dd)$p.value else NA_real_)
}
# x - y with differences at rounding level set to 0, so identical methods (e.g. 做法 2 staying at its ClusterSC start)
# get mean 0 and no t / p instead of a test on float noise
pdiff <- function(x, y) { d <- x - y; d[abs(d) <= 1e-10 * pmax(abs(x), abs(y))] <- 0; d }
CMP <- list(c("jsc1_ff3/avg", "csc_ff3/avg"), c("jsc1_svd/avg", "csc_svd/avg"), c("jsc1_ff3/pf", "csc_ff3/pf"),
            c("jsc2_ff3/pf", "csc_ff3/pf"), c("jsc2_ff3/pf", "jsc1_ff3/pf"),
            c("csc_ff3/avg", "sc/avg"), c("csc_svd/avg", "sc/avg"), c("jsc1_ff3/avg", "sc/avg"),
            c("jsc1_svd/avg", "sc/avg"), c("csc_ff3/pf", "sc/pf"), c("jsc1_ff3/pf", "sc/pf"), c("jsc2_ff3/pf", "sc/pf"),
            c("gsynth/pf", "sc/avg"), c("gsynth/pf", "sc/pf"), c("gsynth/pf", "csc_ff3/avg"), c("gsynth/pf", "csc_ff3/pf"),
            c("jsc1_ff3/avg", "gsynth/pf"), c("jsc1_ff3/pf", "gsynth/pf"), c("jsc2_ff3/pf", "gsynth/pf"))

if (mode == "combine") {
  prefix <- args[2]
  parts <- lapply(args[-(1:2)], readRDS)
  met <- do.call(rbind, lapply(parts, `[[`, "metrics"))
  cho <- do.call(rbind, lapply(parts, `[[`, "choices"))
  dgp <- unique(do.call(rbind, lapply(parts, `[[`, "dgp")))
  gsr <- do.call(rbind, lapply(parts, `[[`, "gsynth"))               # NULL for runs without Gsynth
  met$key <- paste(met$method, met$version, sep = "/")
  dup <- duplicated(met[, c("scenario", "rep", "method", "version")])
  if (any(dup)) stop(sprintf("combine: %d duplicate (scenario, rep, method, version) rows, first: %s", sum(dup),
                             paste(met[which(dup)[1], c("scenario", "rep", "method", "version")], collapse = " ")))
  summ <- do.call(rbind, lapply(split(met, list(met$scenario, met$key), drop = TRUE), function(d) data.frame(
    scenario = d$scenario[1], method = d$key[1], reps = nrow(d), bias0 = mean(d$err0), rmse0 = sqrt(mean(d$err0^2)),
    post_bias = mean(d$post_bias), w_g1 = mean(d$w_g1), ari = mean(d$ari), runtime_median = median(d$runtime),
    mipgap_max = mx(d$mipgap), n_nonopt = sum(d$n_nonopt), converged = mean(d$converged))))
  summ <- summ[order(summ$scenario, summ$rmse0), ]
  # paired within (scenario, rep): squared event-day error of the first minus the second (negative = first better)
  alt_keys <- grep("^jsc2_ff3_c.*_gap.*/pf$", unique(met$key), value = TRUE)
  CMP <- c(CMP, do.call(c, lapply(alt_keys, function(k)
    lapply(c("csc_ff3/pf", "jsc2_ff3/pf", "jsc1_ff3/pf", "sc/pf", "gsynth/pf"), function(b) c(k, b)))))
  paired <- do.call(rbind, lapply(unique(met$scenario), function(s) do.call(rbind, lapply(CMP, function(c2) {
    m <- merge(met[met$scenario == s & met$key == c2[1], c("rep", "err0")],
               met[met$scenario == s & met$key == c2[2], c("rep", "err0")], by = "rep")
    if (nrow(m) < 2) return(NULL)
    cbind(scenario = s, comparison = paste(c2, collapse = " vs "), paired_stats(pdiff(m$err0.x^2, m$err0.y^2)))
  }))))
  # H3: donor-partition ARI of 做法 2 vs k-means (the csc_ff3 rows carry the k-means ARI)
  ari_cmp <- do.call(rbind, lapply(unique(met$scenario), function(s) do.call(rbind, lapply(c("jsc2_ff3/pf", alt_keys), function(k) {
    m <- merge(met[met$scenario == s & met$key == k, c("rep", "ari")],
               met[met$scenario == s & met$key == "csc_ff3/avg", c("rep", "ari")], by = "rep")
    if (nrow(m) < 2) return(NULL)
    cbind(scenario = s, method = k, ari_jsc2 = mean(m$ari.x), ari_kmeans = mean(m$ari.y), paired_stats(pdiff(m$ari.x, m$ari.y)))
  }))))
  write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)
  if (!is.null(paired)) write.csv(paired, paste0(prefix, "_paired.csv"), row.names = FALSE)
  if (!is.null(ari_cmp)) write.csv(ari_cmp, paste0(prefix, "_ari.csv"), row.names = FALSE)
  if (!is.null(cho)) {
    tab <- as.data.frame(table(method = paste(cho$method, cho$version, sep = "/"), m = cho$m, c = cho$c))
    write.csv(tab[tab$Freq > 0, ], paste0(prefix, "_choices.csv"), row.names = FALSE)
  }
  if (!is.null(gsr)) {
    gtab <- as.data.frame(table(scenario = gsr$scenario, r = gsr$r))
    write.csv(gtab[gtab$Freq > 0, ], paste0(prefix, "_gsynth_r.csv"), row.names = FALSE)
  }
  cat("Event-day error (truth 3), percentage points:\n"); print(summ, digits = 3, row.names = FALSE)
  cat("\nPaired differences in squared event-day error; negative = first is better:\n"); print(paired, digits = 3, row.names = FALSE)
  cat("\nDonor-partition ARI, 做法 2 vs k-means:\n"); print(ari_cmp, digits = 3, row.names = FALSE)
  if (!is.null(gsr)) { cat("\nGsynth: factors chosen by cross-validation (reps):\n"); print(table(gsr$scenario, gsr$r)) }
  cat("\nDGP check, correlation of estimated and true loadings:\n")
  print(aggregate(cbind(cor_mkt, cor_smb) ~ scenario, dgp, mean), digits = 3)
  saveRDS(list(metrics = met, choices = cho, dgp = dgp, summary = summ, paired = paired, ari = ari_cmp, gsynth = gsr),
          paste0(prefix, ".rds"))
  quit(save = "no")
}

stopifnot(mode == "run")
opt <- modifyList(list(scenarios = "all", reps = "1-30", pf_reps = "30", grid = "full", tl1 = "60", tl2 = "300", gap = "0.01",
                       jsc2_c = "1", jsc2_alt = "", methods = "all", seed = "2026"), kv(args[-1]))
stopifnot(!is.null(opt$ff), !is.null(opt$out), opt$grid %in% c("full", "small"))
source(file.path(script_dir, "sc_qp.R"))
source(file.path(script_dir, "clustersc.R"))
source(file.path(script_dir, "jointsc.R"))
source(file.path(script_dir, "gsynth_fit.R"))
METHOD_NAMES <- c("sc", "csc_ff3", "csc_svd", "jsc1_ff3", "jsc1_svd", "jsc2", "gsynth")
methods <- if (opt$methods == "all") METHOD_NAMES else strsplit(opt$methods, ",")[[1]]
if (!all(methods %in% METHOD_NAMES)) stop("unknown method(s): ", paste(setdiff(methods, METHOD_NAMES), collapse = ", "),
                                          " (choose from ", paste(METHOD_NAMES, collapse = ", "), ")")
want <- function(m) m %in% methods
if (want("gsynth")) use_gsynth121(script_dir)
t_start <- Sys.time()
log_step <- function(...) { cat(sprintf("[%6.0fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = ""); flush.console() }

SCEN <- expand.grid(sep = c("clear", "blur"), kappa = c(0, 0.5, 1), stringsAsFactors = FALSE)
SCEN$name <- paste(c("none", "weak", "strong")[match(SCEN$kappa, c(0, 0.5, 1))], SCEN$sep, sep = "_")
scen <- if (opt$scenarios == "all") SCEN$name else strsplit(opt$scenarios, ",")[[1]]
stopifnot(all(scen %in% SCEN$name))
rr <- as.integer(strsplit(opt$reps, "-")[[1]]); reps <- rr[1]:rr[length(rr)]
pf_reps <- as.integer(opt$pf_reps); seed <- as.integer(opt$seed)
jsc2_c <- as.numeric(strsplit(opt$jsc2_c, ",")[[1]])
m_grid <- if (opt$grid == "full") c(20, 50, 100) else 50
c_grid <- if (opt$grid == "full") c(0, 0.1, 1, 10) else c(0.1, 1, 10)
p1 <- list(TimeLimit = as.numeric(opt$tl1), MIPGap = as.numeric(opt$gap))
p2 <- list(TimeLimit = as.numeric(opt$tl2), MIPGap = as.numeric(opt$gap))
# 做法 2 arms: name, gamma factor, Gurobi parameters
arms2 <- lapply(jsc2_c, function(cc) list(name = if (length(jsc2_c) == 1) "jsc2_ff3" else paste0("jsc2_ff3_c", cc), c = cc, params = p2))
if (nzchar(opt$jsc2_alt)) {
  alt <- strsplit(opt$jsc2_alt, ":", fixed = TRUE)[[1]]
  stopifnot(length(alt) == 2)
  arms2[[length(arms2) + 1]] <- list(name = sprintf("jsc2_ff3_c%s_gap%s", alt[1], alt[2]), c = as.numeric(alt[1]),
                                     params = modifyList(p2, list(MIPGap = as.numeric(alt[2]))))
}

# ---- Fama-French daily factors, percent (as in 06_simulation_table1.R) -------------------------
ff_lines <- readLines(opt$ff)
ff <- read.csv(text = ff_lines[grep("^,Mkt-RF", ff_lines):length(ff_lines)], strip.white = TRUE)
names(ff) <- c("date", "mkt_rf", "smb", "hml", "rf")
ff <- ff[grepl("^[0-9]{8}$", ff$date), ]
ff$date <- as.Date(as.character(ff$date), "%Y%m%d")
ff[, -1] <- lapply(ff[, -1], as.numeric)
ff <- ff[ff$date >= as.Date("1926-07-01") & ff$date <= as.Date("2022-12-31"), ]
Sig <- cov(ff[, c("mkt_rf", "smb")])
SIGMA <- sqrt(sum(diag(Sig)) * (1 + 0.3^2) + 2 * Sig[1, 2])          # Table 1 calibration
stopifnot(abs(SIGMA - 1.209) < 0.001)
SD_SMB <- sd(ff$smb)

N <- 400; TT <- 250; T0 <- 239; T_EVENT <- 240; EFFECT <- 3; N_TREAT <- 12; K <- 3; R_RANK <- 10L
GROUP <- rep(1:3, c(134, 133, 133))
CENTRE <- rbind(c(0.8, 0.6), c(1.0, 1.2), c(1.3, 0.8))
pre <- seq_len(TT) <= T0; post <- seq_len(TT) >= T_EVENT; after <- seq_len(TT) > T_EVENT
log_step(sprintf("scenarios %s | reps %d-%d | pf_reps %d | grid %s | sigma_e %.3f | sd_smb %.3f",
                 paste(scen, collapse = ","), min(reps), max(reps), pf_reps, opt$grid, SIGMA, SD_SMB))

simulate <- function(si, rep) {
  set.seed(seed + 1000 * si + rep)                                  # data depend only on (scenario, rep)
  start <- sample.int(nrow(ff) - TT + 1, 1)
  Fac <- ff[start:(start + TT - 1), c("mkt_rf", "smb", "rf")]
  sd_w <- if (SCEN$sep[si] == "clear") 0.1 else 0.3
  B <- CENTRE[GROUP, ] + matrix(rnorm(2 * N, 0, sd_w), N)
  G <- matrix(rnorm(K * TT, 0, SCEN$kappa[si] * SD_SMB), K)
  delta <- rnorm(N, 1, 0.3)
  R <- outer(rep(1, N), Fac$rf) + outer(B[, 1], Fac$mkt_rf) + outer(B[, 2], Fac$smb) +
       delta * G[GROUP, ] + matrix(rnorm(N * TT, 0, SIGMA), N)
  treated <- sort(sample(which(GROUP == 1), N_TREAT))
  R[treated, T_EVENT] <- R[treated, T_EVENT] + EFFECT
  list(R = R, Fac = Fac, B = B, treated = treated, donors = setdiff(seq_len(N), treated))
}

timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; v$elapsed <- proc.time()[["elapsed"]] - t0; v }

run_methods <- function(dat, do_pf) {
  R <- dat$R; dn <- dat$donors; tr <- dat$treated
  Y0 <- R[dn, ]; g1 <- GROUP[dn] == 1
  Ff <- features_ff3(R[c(dn, tr), pre], as.matrix(dat$Fac[pre, c("mkt_rf", "smb")]), dat$Fac$rf[pre])
  st <- standardize_features(Ff[seq_along(dn), , drop = FALSE], Ff[-seq_along(dn), , drop = FALSE])
  km <- kmeans_features(st$donor, k = K, seed = 1)
  cl <- cluster_donors(Y0[, pre], r = R_RANK, k = K, seed = 1)
  Fs <- Y0[, pre] %*% cl$V

  sc_on <- function(members, y) {
    f <- simplex_sc(Y0[members, , drop = FALSE], y, pre, post)
    w <- numeric(length(dn)); w[members] <- f$weights
    list(gap = y - f$fit, w = w, mipgap = NA_real_, n_nonopt = 0L)
  }
  jsc1_on <- function(y, d, labels, g) {
    f <- jsc1_tune(Y0, y, pre, post, d, start_order_from_clusters(d, labels, g), m_grid, c_grid, params = p1)
    list(gap = f$gap, w = f$weights, mipgap = mx(c(f$tuning$mipgap, f$mipgap)),
         n_nonopt = sum(c(f$tuning$status, f$status) != "OPTIMAL"), m = f$m, c = f$c)
  }
  one_target <- function(y, f_ff, version) {
    g <- nearest_center(f_ff, km$centers)
    out <- list()
    if (want("sc")) out$sc <- timed(sc_on(seq_along(dn), y))
    if (want("csc_ff3")) out$csc_ff3 <- timed(sc_on(which(km$labels == g), y))
    if (want("jsc1_ff3")) out$jsc1_ff3 <- timed(jsc1_on(y, feature_dist(st$donor, f_ff), km$labels, g))
    if (version == "avg") {
      u <- crossprod(cl$V, y[pre]); gs <- nearest_center(u, cl$centers)
      if (want("csc_svd")) out$csc_svd <- timed(sc_on(which(cl$labels == gs), y))
      if (want("jsc1_svd")) out$jsc1_svd <- timed(jsc1_on(y, feature_dist(Fs, u), cl$labels, gs))
    }
    out
  }
  rows_of <- function(fits, version) do.call(rbind, lapply(names(fits[[1]]), function(mt) {
    f <- lapply(fits, `[[`, mt)
    data.frame(method = mt, version = version,
               err0 = mean(sapply(f, function(x) x$gap[T_EVENT])) - EFFECT,
               post_bias = mean(sapply(f, function(x) mean(x$gap[after]))),
               w_g1 = mean(sapply(f, function(x) sum(x$w[g1]))),
               ari = switch(mt, csc_ff3 = ari(km$labels, GROUP[dn]), csc_svd = ari(cl$labels, GROUP[dn]), NA_real_),
               runtime = sum(sapply(f, `[[`, "elapsed")), mipgap = mx(sapply(f, `[[`, "mipgap")),
               n_nonopt = sum(sapply(f, `[[`, "n_nonopt")), converged = NA, n_iter = NA_integer_)
  }))
  choice_rows <- function(fits, version) do.call(rbind, lapply(seq_along(fits), function(i)
    do.call(rbind, lapply(grep("^jsc1", names(fits[[i]]), value = TRUE), function(mt)
      data.frame(method = mt, version = version, target = i, m = fits[[i]][[mt]]$m, c = fits[[i]][[mt]]$c)))))

  avg <- list(one_target(colMeans(R[tr, ]), colMeans(st$target), "avg"))
  metrics <- rows_of(avg, "avg"); choices <- choice_rows(avg, "avg")
  if (do_pf) {
    pf <- lapply(seq_along(tr), function(i) one_target(R[tr[i], ], st$target[i, ], "pf"))
    metrics <- rbind(metrics, rows_of(pf, "pf")); choices <- rbind(choices, choice_rows(pf, "pf"))
    if (want("jsc2")) {
      init <- list(labels = km$labels, centers = km$centers,
                   target_groups = apply(st$target, 1, nearest_center, centers = km$centers))
      for (arm in arms2) {
        f2 <- timed(jsc2_fit(Y0, R[tr, ], pre, post, st$donor, K,
                             jsc2_gamma(Y0, R[tr, ], pre, st$donor, init, arm$c), init, params = arm$params))
        metrics <- rbind(metrics, data.frame(
          method = arm$name, version = "pf",
          err0 = mean(f2$gaps[, T_EVENT]) - EFFECT, post_bias = mean(f2$gaps[, after]),
          w_g1 = mean(rowSums(f2$weights[, g1, drop = FALSE])), ari = ari(f2$labels, GROUP[dn]),
          runtime = f2$elapsed, mipgap = f2$mipgap, n_nonopt = sum(f2$history$status != "OPTIMAL"),
          converged = f2$converged, n_iter = f2$n_iter))
      }
    }
  }
  gsr <- NULL
  if (want("gsynth")) {                                               # every rep, like the cohort-average methods
    g <- timed(gsynth_att(Y0, R[tr, , drop = FALSE], pre, post))
    metrics <- rbind(metrics, data.frame(method = "gsynth", version = "pf", err0 = mean(g$gaps[, T_EVENT]) - EFFECT,
      post_bias = mean(g$gaps[, after]), w_g1 = NA_real_, ari = NA_real_, runtime = g$elapsed, mipgap = NA_real_,
      n_nonopt = 0L, converged = NA, n_iter = NA_integer_))
    gsr <- data.frame(r = g$r)
  }
  dgp <- data.frame(cor_mkt = cor(Ff[seq_along(dn), 1], dat$B[dn, 1]), cor_smb = cor(Ff[seq_along(dn), 2], dat$B[dn, 2]))
  list(metrics = metrics, choices = choices, dgp = dgp, gsynth = gsr)
}

rows <- list(); choices <- list(); dgps <- list(); gsrs <- list()
for (sname in scen) {
  si <- match(sname, SCEN$name)
  for (rep in reps) {
    res <- run_methods(simulate(si, rep), rep <= pf_reps)
    if (!is.null(res$metrics)) rows[[length(rows) + 1]] <- cbind(scenario = sname, rep = rep, res$metrics)
    if (!is.null(res$choices)) choices[[length(choices) + 1]] <- cbind(scenario = sname, rep = rep, res$choices)
    dgps[[length(dgps) + 1]] <- cbind(scenario = sname, rep = rep, res$dgp)
    if (!is.null(res$gsynth)) gsrs[[length(gsrs) + 1]] <- cbind(scenario = sname, rep = rep, res$gsynth)
    log_step(sprintf("%s rep %d: %s", sname, rep,
                     paste(sprintf("%s/%s %.0fs", res$metrics$method, res$metrics$version, res$metrics$runtime), collapse = ", ")))
    saveRDS(list(metrics = do.call(rbind, rows), choices = do.call(rbind, choices), dgp = do.call(rbind, dgps),
                 gsynth = do.call(rbind, gsrs), opt = opt), opt$out)  # checkpoint after every rep
  }
}
log_step(paste("saved", opt$out))
