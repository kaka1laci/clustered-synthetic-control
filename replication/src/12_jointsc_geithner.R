# 方案 B on Geithner Panel B (or Panel A): placebo cohorts of N1 untreated firms (true effect 0) and
# the real treated cohorts. Estimate = average daily effect over days 0..10 (as in Table 2).
# Methods as in 11_jointsc_sim.R: sc, csc_ff3, csc_svd, jsc1_ff3, jsc1_svd (avg) and sc, csc_ff3, jsc1_ff3,
# jsc2_ff3 (pf). Features: FF3 loadings (Mkt-RF, SMB, HML on days -256..-31), standardized by the untreated
# firms; SVD features r = 10 (the clusters of 09 / 10). Clusters (k by silhouette, k_max 10) are computed once
# on all untreated firms; a placebo draw's own firms leave the donor pool. SC and ClusterSC use all donors;
# 做法 1 uses each target's n_keep most correlated donors, 做法 2 the n_keep with the highest average correlation.
# Controls for the pre-screen (so 做法 vs SC / ClusterSC is not a candidate-set effect): sc_cand = SC on the
# target's own candidates (avg, pf; 做法 1's set), csc_ff3_cand = SC on candidates in the target's start group (pf;
# 做法 2's start). jsc1=0 skips 做法 1 (e.g. to rerun the other methods and combine with an earlier 做法 1 run).
# Placebo draws use the seeds of 10_clustersc_perfirm.R, so the cohorts are the same as that baseline's.
# panel_type=A: Panel A (01_build_geithner_panel.R, 545 firms, all financial); donors = the firms without any
#   connection (497); one placebo pool "all"; no crsp= needed; the financial share is NA. n_keep >= donors keeps every
#   donor in its original order, so sc_cand = sc and csc_ff3_cand = csc_ff3 exactly.
# gsynth: Gsynth (gsynth 1.2.1, unit FE, r by CV over 0..5; gsynth_fit.R) on the cohort's firms and the same donors,
#   every draw; counted as per-firm since it averages per-unit effects.
# methods: run only the listed methods (sc, sc_cand, csc_ff3, csc_svd, csc_ff3_cand, jsc1_ff3, jsc1_svd, jsc2, gsynth);
#   the draws do not depend on it, so methods=gsynth adds Gsynth to earlier runs (combine them).
#
# Usage:
#   run:     Rscript --vanilla replication/src/12_jointsc_geithner.R run panel=<panel_rds> [panel_type=B|A] [crsp=<crsp_dir>]
#              ff=<ff_csv> out=<rds> [pools=all,financial] [draws=1-30] [N1=12] [seed=1] [do_real=1]
#              [pf_draws=30] [grid=full|small] [tl1=60] [tl2=300] [gap=0.01] [n_keep=1000] [jsc1=1] [methods=all]
#   combine: Rscript --vanilla replication/src/12_jointsc_geithner.R combine <out_prefix> <rds>...

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
            # pre-screen controls: the method on the same candidates, and the pre-screen effect itself
            c("jsc1_ff3/avg", "sc_cand/avg"), c("jsc1_svd/avg", "sc_cand/avg"), c("jsc1_ff3/pf", "sc_cand/pf"),
            c("jsc2_ff3/pf", "csc_ff3_cand/pf"), c("sc_cand/avg", "sc/avg"), c("sc_cand/pf", "sc/pf"),
            c("csc_ff3_cand/pf", "csc_ff3/pf"),
            c("gsynth/pf", "sc/avg"), c("gsynth/pf", "sc/pf"), c("gsynth/pf", "csc_ff3/avg"), c("gsynth/pf", "csc_ff3/pf"),
            c("jsc1_ff3/avg", "gsynth/pf"), c("jsc1_ff3/pf", "gsynth/pf"), c("jsc2_ff3/pf", "gsynth/pf"))

if (mode == "combine") {
  prefix <- args[2]
  parts <- lapply(args[-(1:2)], readRDS)
  plac <- do.call(rbind, lapply(parts, `[[`, "placebo"))
  real <- do.call(rbind, lapply(parts, `[[`, "real"))
  cho <- do.call(rbind, lapply(parts, `[[`, "choices"))
  gsr <- do.call(rbind, lapply(parts, `[[`, "gsynth"))                # NULL for runs without Gsynth
  plac$key <- paste(plac$method, plac$version, sep = "/")
  no_dup <- function(d, cols, what) {
    dup <- duplicated(d[, cols])
    if (any(dup)) stop(sprintf("combine: %d duplicate %s rows, first: %s", sum(dup), what,
                               paste(d[which(dup)[1], cols], collapse = " ")))
  }
  no_dup(plac, c("pool", "draw", "method", "version"), "placebo (pool, draw, method, version)")
  if (!is.null(real)) no_dup(real, c("connection", "method", "version"), "real (connection, method, version)")
  summ <- do.call(rbind, lapply(split(plac, list(plac$pool, plac$key), drop = TRUE), function(d) data.frame(
    pool = d$pool[1], method = d$key[1], draws = nrow(d), bias = 100 * mean(d$estimate),
    bias_se = 100 * sd(d$estimate) / sqrt(nrow(d)), rmse = 100 * sqrt(mean(d$estimate^2)),
    pre_rmse = 100 * mean(d$pre_rmse), fin_share = mean(d$fin), runtime_median = median(d$runtime),
    mipgap_max = mx(d$mipgap), n_nonopt = sum(d$n_nonopt), converged = mean(d$converged))))
  summ <- summ[order(summ$pool, summ$rmse), ]
  paired <- do.call(rbind, lapply(unique(plac$pool), function(pool) do.call(rbind, lapply(CMP, function(c2) {
    m <- merge(plac[plac$pool == pool & plac$key == c2[1], c("draw", "estimate")],
               plac[plac$pool == pool & plac$key == c2[2], c("draw", "estimate")], by = "draw")
    if (nrow(m) < 2) return(NULL)
    lab <- paste(c2, collapse = " vs ")
    rbind(cbind(pool = pool, comparison = lab, measure = "squared error (pp^2)", paired_stats(pdiff(1e4 * m$estimate.x^2, 1e4 * m$estimate.y^2))),
          cbind(pool = pool, comparison = lab, measure = "estimate (pp)", paired_stats(pdiff(100 * m$estimate.x, 100 * m$estimate.y))))
  }))))
  write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)
  if (!is.null(paired)) write.csv(paired, paste0(prefix, "_paired.csv"), row.names = FALSE)
  if (!is.null(real)) write.csv(transform(real, estimate = 100 * estimate, pre_rmse = 100 * pre_rmse),
                                paste0(prefix, "_real.csv"), row.names = FALSE)
  if (!is.null(cho)) {
    tab <- as.data.frame(table(method = paste(cho$method, cho$version, sep = "/"), m = cho$m, c = cho$c))
    write.csv(tab[tab$Freq > 0, ], paste0(prefix, "_choices.csv"), row.names = FALSE)
  }
  if (!is.null(gsr)) {
    gtab <- as.data.frame(table(pool = gsr$pool, r = gsr$r))
    write.csv(gtab[gtab$Freq > 0, ], paste0(prefix, "_gsynth_r.csv"), row.names = FALSE)
    cat("Gsynth: factors chosen by cross-validation (draws):\n"); print(table(gsr$pool, gsr$r))
  }
  cat("Placebo (true effect 0), percent per day:\n"); print(summ, digits = 3, row.names = FALSE)
  cat("\nPaired differences (negative = first has lower squared error / lower estimate):\n"); print(paired, digits = 3, row.names = FALSE)
  if (!is.null(real)) { cat("\nReal treated cohorts, percent per day:\n"); print(transform(real, estimate = round(100 * estimate, 3)), row.names = FALSE) }
  saveRDS(list(placebo = plac, real = real, choices = cho, summary = summ, paired = paired, gsynth = gsr),
          paste0(prefix, ".rds"))
  quit(save = "no")
}

stopifnot(mode == "run")
user_opt <- kv(args[-1])
opt <- modifyList(list(panel_type = "B", pools = "all,financial", draws = "1-30", N1 = "12", seed = "1", do_real = "1",
                       pf_draws = "30", grid = "full", tl1 = "60", tl2 = "300", gap = "0.01", n_keep = "1000", jsc1 = "1",
                       methods = "all"), user_opt)
if (opt$panel_type == "A" && is.null(user_opt$pools)) opt$pools <- "all"
stopifnot(opt$panel_type %in% c("A", "B"), !is.null(opt$panel), !is.null(opt$ff), !is.null(opt$out),
          opt$grid %in% c("full", "small"))
if (opt$panel_type == "B" && is.null(opt$crsp)) stop("panel_type=B needs crsp=<crsp_dir> (SIC codes for the financial pool)")
dr <- as.integer(strsplit(opt$draws, "-")[[1]]); draws <- dr[1]:dr[length(dr)]
N1 <- as.integer(opt$N1); seed <- as.integer(opt$seed); pf_draws <- as.integer(opt$pf_draws); N_KEEP <- as.integer(opt$n_keep)
pool_names <- strsplit(opt$pools, ",")[[1]]
if (opt$panel_type == "A" && !identical(pool_names, "all")) stop("Panel A has one pool (all its firms are financial): use pools=all")
do_jsc1 <- opt$jsc1 == "1"
METHOD_NAMES <- c("sc", "sc_cand", "csc_ff3", "csc_svd", "csc_ff3_cand", "jsc1_ff3", "jsc1_svd", "jsc2", "gsynth")
methods <- if (opt$methods == "all") METHOD_NAMES else strsplit(opt$methods, ",")[[1]]
if (!all(methods %in% METHOD_NAMES)) stop("unknown method(s): ", paste(setdiff(methods, METHOD_NAMES), collapse = ", "),
                                          " (choose from ", paste(METHOD_NAMES, collapse = ", "), ")")
want <- function(m) m %in% methods && (do_jsc1 || !grepl("^jsc1", m))
sel <- Filter(want, METHOD_NAMES)
if (!length(sel) || (all(sel %in% c("csc_ff3_cand", "jsc2")) && opt$do_real != "1" && !any(draws <= pf_draws)))
  stop("no method selected produces rows (methods=", opt$methods, ", jsc1=", opt$jsc1, ", pf_draws=", pf_draws,
       "): csc_ff3_cand and jsc2 run on per-firm draws only, and jsc1=0 drops jsc1_*")
m_grid <- if (opt$grid == "full") c(20, 50, 100) else 50
c_grid <- if (opt$grid == "full") c(0, 0.1, 1, 10) else c(0.1, 1, 10)
p1 <- list(TimeLimit = as.numeric(opt$tl1), MIPGap = as.numeric(opt$gap))
p2 <- list(TimeLimit = as.numeric(opt$tl2), MIPGap = as.numeric(opt$gap))
R_RANK <- 10L

source(file.path(script_dir, "sc_qp.R"))
source(file.path(script_dir, "clustersc.R"))
source(file.path(script_dir, "jointsc.R"))
source(file.path(script_dir, "gsynth_fit.R"))
if (want("gsynth")) use_gsynth121(script_dir)
t_start <- Sys.time()
log_step <- function(...) { cat(sprintf("[%6.0fs] ", as.numeric(difftime(Sys.time(), t_start, units = "secs"))), ..., "\n", sep = ""); flush.console() }

p <- readRDS(opt$panel)
pre <- p$dif < 0; post <- p$dif >= 0
Y <- p$Y
ctrl <- which(rowSums(p$treat) == 0)
if (opt$panel_type == "B") {
  nm <- readRDS(file.path(opt$crsp, "crsp_names.rds"))
  nm$namedt <- as.Date(nm$namedt); nm$nameendt <- as.Date(nm$nameendt)
  ev <- nm[nm$namedt <= as.Date("2008-11-21") & nm$nameendt >= as.Date("2008-11-21"), ]
  sic <- ev$siccd[match(as.integer(p$permno), ev$permno)]
  financial <- !is.na(sic) & sic >= 6000 & sic <= 6799
} else financial <- rep(NA, nrow(Y))                                  # Panel A: every firm is financial

ff_lines <- readLines(opt$ff)
ff <- read.csv(text = ff_lines[grep("^,Mkt-RF", ff_lines):length(ff_lines)], strip.white = TRUE)
names(ff) <- c("date", "mkt_rf", "smb", "hml", "rf")
ff <- ff[grepl("^[0-9]{8}$", ff$date), ]
ff$date <- as.Date(as.character(ff$date), "%Y%m%d")
fac <- sapply(ff[match(p$dates, ff$date), c("mkt_rf", "smb", "hml", "rf")], as.numeric) / 100   # decimals, like Y
stopifnot(!anyNA(fac))

Ff <- features_ff3(Y[, pre], fac[pre, c("mkt_rf", "smb", "hml")], fac[pre, "rf"])
Fall <- standardize_features(Ff[ctrl, ], Ff)$target                    # every firm, scaled by the untreated firms
km <- kmeans_features(Fall[ctrl, ], k = NULL, k_max = 10, seed = seed)
cl <- cluster_donors(Y[ctrl, pre], r = R_RANK, k_max = 10, seed = seed)  # the r = 10 clusters of 09 / 10
Fs <- Y[, pre] %*% cl$V
lab_ff <- lab_sv <- integer(nrow(Y)); lab_ff[ctrl] <- km$labels; lab_sv[ctrl] <- cl$labels
fin_note <- if (opt$panel_type == "B") sprintf(" (%d financial)", sum(financial[ctrl])) else " (Panel A, all financial)"
fin_share <- if (opt$panel_type == "B") paste0(", financial share ", paste(round(tapply(financial[ctrl], km$labels, mean), 2),
                                                                         collapse = "/")) else ""
log_step(sprintf("Panel %s: %d untreated firms%s | FF3 clusters k = %d, sizes %s%s | SVD clusters k = %d", opt$panel_type,
                 length(ctrl), fin_note, km$k, paste(table(km$labels), collapse = "/"), fin_share, cl$k))

timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- expr; v$elapsed <- proc.time()[["elapsed"]] - t0; v }

estimate_cohort <- function(firms, donors, do_pf) {
  cand_of <- function(target_pre) if (N_KEEP >= length(donors)) donors else donors[prescreen_corr(Y[donors, pre], target_pre, N_KEEP)]
  sc_on <- function(members, y) {
    f <- simplex_sc(Y[members, , drop = FALSE], y, pre, post)
    list(estimate = f$estimate, pre_rmse = f$pre_rmse, fin = sum(f$weights * financial[members]),
         mipgap = NA_real_, n_nonopt = 0L)
  }
  jsc1_on <- function(y, cand, d, labels, g) {
    f <- jsc1_tune(Y[cand, ], y, pre, post, d, start_order_from_clusters(d, labels, g), m_grid, c_grid, params = p1)
    list(estimate = f$estimate, pre_rmse = f$pre_rmse, fin = sum(f$weights * financial[cand]),
         mipgap = mx(c(f$tuning$mipgap, f$mipgap)), n_nonopt = sum(c(f$tuning$status, f$status) != "OPTIMAL"),
         m = f$m, c = f$c)
  }
  one_target <- function(y, f_ff, version) {
    need_cand <- want("sc_cand") || want("jsc1_ff3") || (version == "avg" && want("jsc1_svd"))
    cand <- if (need_cand) cand_of(y[pre]) else NULL
    g <- nearest_center(f_ff, km$centers, allowed = sort(unique(lab_ff[donors])))
    out <- list()
    if (want("sc")) out$sc <- timed(sc_on(donors, y))
    if (want("sc_cand")) out$sc_cand <- timed(sc_on(cand, y))
    if (want("csc_ff3")) out$csc_ff3 <- timed(sc_on(donors[lab_ff[donors] == g], y))
    if (want("jsc1_ff3")) out$jsc1_ff3 <- timed(jsc1_on(y, cand, feature_dist(Fall[cand, , drop = FALSE], f_ff), lab_ff[cand], g))
    if (version == "avg") {
      u <- crossprod(cl$V, y[pre])
      gs <- nearest_center(u, cl$centers, allowed = sort(unique(lab_sv[donors])))
      if (want("csc_svd")) out$csc_svd <- timed(sc_on(donors[lab_sv[donors] == gs], y))
      if (want("jsc1_svd")) out$jsc1_svd <- timed(jsc1_on(y, cand, feature_dist(Fs[cand, , drop = FALSE], u), lab_sv[cand], gs))
    }
    out
  }
  rows_of <- function(fits, version) do.call(rbind, lapply(names(fits[[1]]), function(mt) {
    f <- lapply(fits, `[[`, mt)
    data.frame(method = mt, version = version, estimate = mean(sapply(f, `[[`, "estimate")),
               pre_rmse = mean(sapply(f, `[[`, "pre_rmse")), fin = mean(sapply(f, `[[`, "fin")),
               runtime = sum(sapply(f, `[[`, "elapsed")), mipgap = mx(sapply(f, `[[`, "mipgap")),
               n_nonopt = sum(sapply(f, `[[`, "n_nonopt")), converged = NA, n_iter = NA_integer_)
  }))
  choice_rows <- function(fits, version) do.call(rbind, lapply(seq_along(fits), function(i)
    do.call(rbind, lapply(grep("^jsc1", names(fits[[i]]), value = TRUE), function(mt)
      data.frame(method = mt, version = version, target = i, m = fits[[i]][[mt]]$m, c = fits[[i]][[mt]]$c)))))

  avg <- list(one_target(colMeans(Y[firms, , drop = FALSE]), colMeans(Fall[firms, , drop = FALSE]), "avg"))
  rows <- rows_of(avg, "avg"); choices <- choice_rows(avg, "avg")
  if (do_pf) {
    pf <- lapply(firms, function(i) one_target(Y[i, ], Fall[i, ], "pf"))
    rows <- rbind(rows, rows_of(pf, "pf")); choices <- rbind(choices, choice_rows(pf, "pf"))
    if (want("csc_ff3_cand") || want("jsc2")) {
      cand <- cand_of(Y[firms, pre, drop = FALSE])
      allowed <- sort(unique(lab_ff[cand]))
      init <- list(labels = lab_ff[cand], centers = km$centers,
                   target_groups = sapply(firms, function(i) nearest_center(Fall[i, ], km$centers, allowed)))
      Yt <- Y[firms, , drop = FALSE]; Fc <- Fall[cand, , drop = FALSE]
      if (want("csc_ff3_cand")) {
        start <- lapply(seq_along(firms), function(t) list(csc_ff3_cand = timed(sc_on(cand[init$labels == init$target_groups[t]], Yt[t, ]))))
        rows <- rbind(rows, rows_of(start, "pf"))
      }
      if (want("jsc2")) {
        f2 <- timed(jsc2_fit(Y[cand, ], Yt, pre, post, Fc, km$k, jsc2_gamma(Y[cand, ], Yt, pre, Fc, init, 1), init, params = p2))
        rows <- rbind(rows, data.frame(method = "jsc2_ff3", version = "pf", estimate = mean(f2$estimates),
          pre_rmse = mean(sqrt(rowMeans(f2$gaps[, pre, drop = FALSE]^2))), fin = mean(f2$weights %*% financial[cand]),
          runtime = f2$elapsed, mipgap = f2$mipgap, n_nonopt = sum(f2$history$status != "OPTIMAL"),
          converged = f2$converged, n_iter = f2$n_iter))
      }
    }
  }
  gs_r <- NULL
  if (want("gsynth")) {                                               # every draw, like the cohort-average methods
    g <- timed(gsynth_att(Y[donors, , drop = FALSE], Y[firms, , drop = FALSE], pre, post))
    rows <- rbind(rows, data.frame(method = "gsynth", version = "pf", estimate = g$estimate, pre_rmse = g$pre_rmse,
      fin = NA_real_, runtime = g$elapsed, mipgap = NA_real_, n_nonopt = 0L, converged = NA, n_iter = NA_integer_))
    gs_r <- g$r
  }
  list(rows = rows, choices = choices, gsynth_r = gs_r)
}

save_all <- function() saveRDS(list(placebo = plac, real = real, choices = choices, gsynth = gsr, opt = opt,
                                     k_ff = km$k, k_svd = cl$k), opt$out)
real <- NULL; plac <- NULL; choices <- NULL; gsr <- NULL
if (opt$do_real == "1") {
  for (cn in names(p$treat)) {
    in_sample <- if (is.null(p$sample)) rep(TRUE, nrow(Y)) else p$sample[[cn]]   # Panel A has no per-connection sample
    tr <- which(p$treat[[cn]] == 1 & in_sample)
    r <- estimate_cohort(tr, ctrl, TRUE)
    if (!is.null(r$rows)) real <- rbind(real, cbind(connection = cn, n_treated = length(tr), r$rows))
    if (!is.null(r$choices)) choices <- rbind(choices, cbind(pool = paste0("real_", cn), draw = 0L, r$choices))
    if (!is.null(r$gsynth_r)) gsr <- rbind(gsr, data.frame(pool = paste0("real_", cn), draw = 0L, r = r$gsynth_r))
    log_step(sprintf("real %s (%d firms): %s", cn, length(tr), paste(sprintf("%s/%s %.0fs", r$rows$method, r$rows$version, r$rows$runtime), collapse = ", ")))
    save_all()
  }
}
pools <- if (opt$panel_type == "B") list(all = ctrl, financial = ctrl[financial[ctrl]]) else list(all = ctrl)
for (pool in pool_names) for (b in draws) {
  set.seed(seed + match(pool, names(pools)) * 100000 + b)            # same draws as 10_clustersc_perfirm.R
  pseudo <- sample(pools[[pool]], N1)
  r <- estimate_cohort(pseudo, setdiff(ctrl, pseudo), b <= pf_draws)
  if (!is.null(r$rows)) plac <- rbind(plac, cbind(pool = pool, draw = b, r$rows))   # per-firm-only methods skip later draws
  if (!is.null(r$choices)) choices <- rbind(choices, cbind(pool = pool, draw = b, r$choices))
  if (!is.null(r$gsynth_r)) gsr <- rbind(gsr, data.frame(pool = pool, draw = b, r = r$gsynth_r))
  log_step(sprintf("placebo %s draw %d: %s", pool, b, paste(sprintf("%s/%s %.0fs", r$rows$method, r$rows$version, r$rows$runtime), collapse = ", ")))
  save_all()
}
save_all()
log_step(paste("saved", opt$out))
