# Replication: Goldsmith-Pinkham & Lyu, *Causal Inference in Financial Event Studies*

Targets: **Table 2** (Geithner as Treasury Secretary) — Panel A (bank controls) and Panel B
(all NYSE/AMEX/NASDAQ common stocks as controls).
The authors have not released code; this rebuilds it from public data, CRSP, and packages.

## Layout

```
replication/
├── data/raw/kermani/   Acemoglu et al. (2016) replication package (data only; their .do/.m files are reference, never run)
├── data/raw/ff/        Fama-French daily 3 factors (Kenneth French data library)
├── data/raw/crsp/      CRSP extracts from WRDS (licensed: keep local, do not share)
├── rlib/gsynth121/     gsynth 1.2.1 (the version the paper cites) + lfe
├── src/                our scripts
├── tests/              checks for the SC solver
└── output/             panels (.rds), result tables, logs, diagnostics/
```

## Run (from the project root)

```bash
# Panel A
Rscript --vanilla replication/src/01_build_geithner_panel.R \
  "replication/data/raw/kermani/Stata Files" \
  replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv \
  replication/output/geithner_panel.rds
Rscript --vanilla replication/src/02_estimate_table2.R \
  replication/output/geithner_panel.rds replication/output/table2_panelA_schedule.csv schedule 100 1000
#   args: <panel> <out_csv> <connections: all | schedule,personal,newyork> <placebo reps> <gsynth bootstraps>
#   run one process per connection to parallelise

# Panel B (needs a WRDS account; password in ~/.pgpass, chmod 600)
Rscript --vanilla replication/src/04_pull_crsp.R replication/data/raw/crsp
Rscript --vanilla replication/src/05_build_panelB.R replication/data/raw/crsp \
  replication/output/geithner_panel.rds replication/output/geithner_panelB_complete.rds complete
Rscript --vanilla replication/src/02_estimate_table2.R \
  replication/output/geithner_panelB_complete.rds replication/output/table2_panelB_schedule.csv schedule 100 1000

# Comparison tables
Rscript --vanilla replication/src/03_make_table2_comparison.R \
  replication/output/table2_panelA_comparison.md replication/output/table2_panelA_{schedule,personal,newyork}.csv

# Tests
Rscript --vanilla replication/tests/test_sc_qp.R replication/src replication/output/geithner_panel.rds
Rscript --vanilla replication/tests/test_panelA_regression.R <dir_with_new_csvs> replication/output   # after an upgrade
```

R 4.6.1 (CRAN arm64). R packages: `haven`, `quadprog`, `DBI`, `RPostgres`, `synthdid` (GitHub: synth-inference/synthdid,
commit 70c1ce3), `gurobi` 13.0-3 + `slam` (方案 B), and **gsynth 1.2.1** + `lfe` in `replication/rlib/gsynth121`
(the scripts put it first on the library path). `rlib/gsynth121-R4.5` is the same build for R 4.5.1, kept as a fallback
(switch `/Library/Frameworks/R.framework/Versions/Current` back to `4.5-arm64` and rename the folder).
Paper values for both panels are in `src/paper_table2.csv`.

## What each estimator is (pinned down by matching the paper's point estimates)

| Column | Implementation | Note |
|---|---|---|
| Average | mean daily return, treated minus control, days 0..10 | SE: two-sample t-test on firm-days |
| DID, SDID | `synthdid::did_estimate`, `synthdid::synthdid_estimate` | placebo SE, 100 reps |
| Market | mean of (R − R_m) over **treated firms only** | Definition 3: no control group is subtracted, so Panel B = Panel A |
| CAPM | treated-only mean of R − α̂ − β̂·R_m, β from firm OLS on days **−280..−31** | the Acemoglu window (as in GP&L §5.1.2), not the −256..−31 stated for Table 2; −256 gives 0.0154/0.0151/0.0108 vs paper 0.016/0.016/0.011. S&P 500 = Datastream row 679; day 0 = 3–4pm return |
| FF3F | treated-only mean abnormal return on Mkt-RF, SMB, HML | published full-day factors on day 0 |
| SC | `sc_qp.R`: the `synthdid::sc_estimate` objective solved with `quadprog` | synthdid's Frank-Wolfe needs ~15 min per fit on Panel A; QP gives the same 0.0161 in 0.5 s. For > 600 donors a working-set method reaches the same optimum (tested) |
| Gsynth | `gsynth` **1.2.1**, `force = "unit"`, CV over r = 0..5, parametric bootstrap 1000 | gsynth 1.4 (fect wrapper) picks r = 0/1 by CV and gives 0.045 / 0.020 for schedule; 1.2.1 picks r = 3 and gives 0.0127 (paper 0.012) |

Building gsynth 1.2.1 on current R needed two edits to `src/Makevars`: `CXX_STD = CXX17`
(RcppArmadillo now requires C++14+) and `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS)` (no gfortran
toolchain installed; R's own LAPACK already links R's libgfortran). It also needs `lfe` in the same library.
Rebuilt the same way for R 4.6.1 (2026-10); Panel A point estimates are unchanged (`test_panelA_regression.R`).

## Data sources

- Acemoglu, Johnson, Kermani, Kwak & Mitton (2016) data and code:
  https://faculty.haas.berkeley.edu/amir/research/data%20and%20code/Data%20and%20Programs.zip
- Fama-French daily factors:
  https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip
- CRSP via WRDS: `crsp.dsenames`, `crsp.dsf`, `crsp.dsi`, 2007-09-01 to 2008-12-31 (`04_pull_crsp.R`)

## Panel A sample (matches the paper's N = 545 firms x 237 days)

1. Acemoglu et al. filters: >= 225 return obs in days -280..-31 and all 11 obs in days 0..10 (666 firms).
2. Day-0 return = 3pm-4pm return on Nov 21, 2008 (`Return Nov. 21 Data.dta`); drop firms without it (603).
3. Drop firms whose return correlation with Citigroup exceeds the 90th percentile; the
   percentile is computed over `fullsample == 1` firms, the cutoff applied to all firms (545).
   Requiring `fullsample == 1` as well gives 525, which does not match the paper.
4. GP&L window: pre days -256..-31 (226) + event days 0..10 (11); balanced for all 545 firms.

Treated firms: schedule 15, personal 9, New York 38. Controls for each = all other firms in the 545.

## Panel B sample (matches the paper's N = 4,078 firms x 237 days)

1. CRSP share code 10/11 and exchange code 1–3 on Nov 21, 2008 (4,366 securities).
2. A CRSP return on every one of the 237 window days (4,078). All returns are CRSP close-to-close,
   treated firms included: there is no 3–4pm day-0 return in Panel B.
3. Treated = the Acemoglu connected firms, matched to PERMNO by ticker on the event date and kept
   if they are common stock: AB, BX, IEP, OZM (shrcd 71) and LAZ (shrcd 12) drop out.
   With rule `complete` (main): schedule 12, personal 7, New York 32.
4. Rule `allow_missing` also keeps treated firms that miss a day (CARV and CMSB, both 2008-03-10;
   filled with 0): personal 8, New York 34. It reproduces the paper's personal Average/DID
   (0.0199 / 0.0206 vs 0.020 / 0.021) but not New York's (0.0084 vs 0.011); rule `complete` does the
   reverse (personal 0.0274, New York 0.0106). No single rule matches all three connections.
5. Not explained: the paper's Average column has 45,045 obs (4,095 firms x 11), not 4,078 x 11.

## Table 1 simulation (`06_simulation_table1.R`, comparison `07_make_table1_comparison.R`)

```bash
for p in A B C D; do Rscript --vanilla replication/src/06_simulation_table1.R \
  replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv replication/output/sim_table1_$p.rds $p 50 200 & done; wait
Rscript --vanilla replication/src/07_make_table1_comparison.R replication/output/table1_comparison.md \
  replication/output/sim_table1_{A,B,C,D}.rds
```

Definitions the paper does not state, reverse-engineered from its Table 1 numbers (details in the script header):
units are percentage points; "All" / "Untreated" E(Bias) and MAD are for the average over the 11 / last 10
periods; RMSE is for the cumulative (summed) error; "Coverage" is the rejection rate of a 5% test.
sigma_e = 1.209 pp makes the residual variance half of the total, as the paper states (its literal "0.1" does not).

Result (50 sims per panel, gsynth bootstrap 200; `output/table1_comparison.md`): rejection rates and the
pattern of biases match — every estimator is unbiased under random timing; under timing selection only CAPM
is biased; under both selections Simple Means is biased too; CAPM over-rejects on untreated days (0.33–0.43 vs
paper 0.40–0.46). **The size of the timing-selection bias is about half the paper's**: CAPM event-day bias
1.23 (C) / 1.00 (D) vs 2.71 / 2.26, Simple Means (D) −0.29 vs −0.52. Both scale with the event-day SMB return
(CAPM bias ≈ 0.9 × SMB, corr 0.88). Our selection rule picks roughly the largest SMB day of the 250-day block
(median 1.2 pp); the paper's numbers imply an event-day SMB near 3 pp, i.e. a more extreme selection than the
paper's description produces.

## Panel B results (`output/table2_panelB_complete_comparison.md`, `..._allow_missing_comparison.md`)

Placebo SE 100 reps, gsynth bootstrap 1000. The paper's qualitative finding holds at the 5% level: with
all-CRSP controls, SC, SDID and Gsynth are insignificant for every connection under both treated rules,
while Average and DID stay significant. The one exception is schedule-connection Gsynth (0.0132, SE 0.0078,
significant at 10% only; paper 0.001). Standard errors are close to the paper's (schedule: SC 0.0070 vs
0.007, SDID 0.0064 vs 0.006, Gsynth 0.0078 vs 0.008). Point estimates of the weighted estimators differ by
0.001–0.004 (SC, SDID) and up to 0.012 (Gsynth); the cause is not yet known.
Runtime note: synthdid placebo on 4,078 firms took ~125 s per refit with five processes in parallel.

## S&P 500 inclusions, 1980s (`08_sp500_1980s.R`; summary `output/sp500_1980s_comparison.md`)

```bash
Rscript --vanilla replication/src/04_pull_crsp.R replication/data/raw/crsp_1979_1990 1978-12-01 1990-03-31
Rscript --vanilla replication/src/08_sp500_1980s.R replication/data/raw/crsp_1979_1990 \
  replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv replication/output/sp500_1980s.rds 1
```

Only 1980 – Sep 1989 is done: GP&L infer those announcement dates from the Wednesday-announce /
Thursday-effective rule, which CRSP effective dates reproduce. Later decades need Siblis announcement dates
(effective dates do not locate the announcement jump after 1989; `check_sp500_effective_dates.R`).
Diff-in-Means 3.25% vs 3.27%; model columns 0.08–0.20 high; Gsynth removes ~75% of the pre-inclusion drift.
gsynth 1.2.1 note: with `se = FALSE`, `fit$att` is named by the panel's time index, not relative time.

## ClusterSC baseline on Panel B (`clustersc.R`, `09_clustersc_panelB.R`, `10_clustersc_perfirm.R`)

```bash
Rscript --vanilla replication/tests/test_clustersc.R replication/src     # reproduces ClusterSC's gain on its own DGP
Rscript --vanilla replication/src/09_clustersc_panelB.R replication/output/geithner_panelB_complete.rds \
  replication/data/raw/crsp replication/output/clustersc_panelB.rds 100 12 1   # ~45 min
```

Cohort-average target, 6 estimators ({simplex, robust} x {all donors, assigned cluster, random subset}),
placebo cohorts of 12 untreated firms (all / financial), 100 draws each. Daily returns are far from low
rank (95% of singular values needs 207 of 226 components), clusters do not separate industries, and
clustering never improves placebo error significantly; it is significantly worse in two comparisons.
Financial placebo cohorts show a common bias of about −0.25 pp/day under every estimator.
Per-firm version (`10_clustersc_perfirm.R`, as in the ClusterSC paper): clustering moves placebo error
slightly in the right direction (t between −1.5 and −1.8) but never significantly; the financial bias
persists; per-firm averaging has lower placebo RMSE than cohort-average SC (not a paired comparison).

## 方案 B: joint clustering + SC prototype (`jointsc.R`, `11_jointsc_sim.R`, `12_jointsc_geithner.R`; results `output/jointsc_sim_*.csv`, `output/jointsc_geithner_*.csv`)

```bash
Rscript --vanilla replication/tests/test_jointsc.R replication/src      # Gurobi smoke test + 做法 1 / 2 correctness
Rscript --vanilla replication/tests/test_jointsc_combine.R replication/src <scratch_dir>   # identical methods get no t-test
# simulation: one process per scenario pair (weak_clear,weak_blur -> _p2.rds; strong_clear,strong_blur -> _p3.rds)
Rscript --vanilla replication/src/11_jointsc_sim.R run ff=replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv \
  out=replication/output/jointsc_sim_p1.rds scenarios=none_clear,none_blur pf_reps=10 grid=small jsc2_alt=0.1:0.001
Rscript --vanilla replication/src/11_jointsc_sim.R combine replication/output/jointsc_sim replication/output/jointsc_sim_p*.rds
# Geithner: full runs (all pool + real cohorts; financial pool with do_real=0 -> jointsc_geithner_financial.rds) ...
Rscript --vanilla replication/src/12_jointsc_geithner.R run panel=replication/output/geithner_panelB_complete.rds \
  crsp=replication/data/raw/crsp ff=replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv \
  out=replication/output/jointsc_geithner_all.rds pools=all do_real=1 pf_draws=10 grid=small
# ... then both pools again with jsc1=0 after the review fixes (-> jointsc_geithner_{all,financial}_v2.rds);
# keep only the 做法 1 rows of the full runs (its code did not change; same seeded cohorts) and combine
Rscript --vanilla -e 'for (p in c("all", "financial")) { x <- readRDS(sprintf("replication/output/jointsc_geithner_%s.rds", p)); k <- function(d) if (is.null(d)) NULL else d[grepl("^jsc1_", d$method), ]; saveRDS(list(placebo = k(x$placebo), real = k(x$real), choices = x$choices, opt = x$opt), sprintf("replication/output/jointsc_geithner_%s_jsc1.rds", p)) }'
Rscript --vanilla replication/src/12_jointsc_geithner.R combine replication/output/jointsc_geithner \
  replication/output/jointsc_geithner_{all,financial}_jsc1.rds replication/output/jointsc_geithner_{all,financial}_v2.rds
```

做法 1 (`jsc1_fit`, `jsc1_tune`): per target, a convex MIQP chooses a group of m donors and simplex weights jointly,
with a penalty lambda on the group's feature distance to the target (m and lambda tuned on the last 20% of the
pre-period). 做法 2 (`jsc2_fit`): one partition of the donors into K groups shared by all targets, each target using
only its group's donors; alternates a joint MIQP (partition + all weights, given centres) with a centre update.
Features: FF3 loadings (main) or SVD (r = 10). Gurobi 13.0.3, MIPGap 1%, 60 s per 做法 1 fit, 300 s per 做法 2 iteration.
Compute level L2: per-firm versions on the first 10 reps / draws, 做法 1 grid m = 50, c in {0.1, 1, 10}; one extra
做法 2 arm in the simulation (c = 0.1, MIPGap 0.1%). Time limits are wall-clock, so keep the machine awake
(we ran each process under `caffeinate -i`, on AC power).
Result: neither approach beats SC or ClusterSC consistently, and hypotheses H1–H3 do not hold. 做法 1's MIQP often
stops at the 60 s limit: about a third of cohort-average fits end above a 1% gap, even with only 388 donors. Under
the current gamma normalisation, the cluster term is about c·n/N1 times the fit term at the start, so 做法 2 never leaves its
k-means start in the simulation and equals per-firm ClusterSC; on Geithner it moves and becomes the worst method
(financial-pool RMSE 1.09 pp/day). On Geithner the 1,000-donor correlation pre-screen alone lowers placebo
estimates by 0.06–0.37 pp/day; net of it, cohort-average 做法 1 matches SC on the same candidates, and per-firm
做法 1 shifts estimates up by about 0.23 in both pools.

## Phase 5: Panel A and the Gsynth comparison (`gsynth_fit.R`; results `output/jointsc_geithnerA_*.csv` and the Gsynth rows of the 方案 B CSVs)

```bash
Rscript --vanilla replication/tests/test_gsynth_fit.R replication/src replication/output     # Gsynth wrapper = Table 2 Gsynth
Rscript --vanilla replication/tests/test_phase5_modes.R replication/src replication/output \
  replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv replication/data/raw/crsp <scratch_dir>
# Gsynth added to the earlier runs (the data / draws depend only on the seeds), then combined with them
Rscript --vanilla replication/src/11_jointsc_sim.R run ff=replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv \
  out=replication/output/jointsc_sim_gsynth.rds methods=gsynth
Rscript --vanilla replication/src/12_jointsc_geithner.R run panel=replication/output/geithner_panelB_complete.rds \
  crsp=replication/data/raw/crsp ff=replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv \
  out=replication/output/jointsc_geithner_gsynth.rds methods=gsynth do_real=1
# Panel A, all methods (three processes: draws=1-5 do_real=1, draws=6-10, draws=11-30)
Rscript --vanilla replication/src/12_jointsc_geithner.R run panel=replication/output/geithner_panel.rds panel_type=A \
  ff=replication/data/raw/ff/F-F_Research_Data_Factors_daily.csv out=replication/output/jointsc_geithnerA_p1.rds \
  draws=1-5 pf_draws=10 grid=small do_real=1
Rscript --vanilla replication/src/12_jointsc_geithner.R combine replication/output/jointsc_geithnerA replication/output/jointsc_geithnerA_p1.rds replication/output/jointsc_geithnerA_p2.rds replication/output/jointsc_geithnerA_p3.rds
```

Gsynth (`gsynth_att`): gsynth 1.2.1 with the Table 2 settings (unit fixed effects, r by cross-validation over 0..5),
point estimates only; it reproduces the replicated Table 2 Gsynth exactly (Panel A 0.0127, Panel B 0.0132 for
schedule). It averages per-firm effects, so it is reported as `gsynth/pf` and run on every draw.
Panel A (`panel_type=A`): donors are the 497 firms without any connection; one placebo pool; methods and parameters as
in the 方案 B runs.
Result: on Panel A, clustering and joint group selection again do not lower placebo error (no squared-error
difference against SC is significant), although they shift the cohort-average estimates up by 0.10–0.26 pp/day.
In the simulation, whose data come from a linear factor model, Gsynth is never significantly worse than the methods
it was paired with and often significantly better, but it is indistinguishable from SC on both Geithner panels. The
negative placebo bias seen for financial firms in Panel B does not go away on Panel A, where placebo firms and donors
are all banks (SC: −0.12 cohort-average, not distinguishable from 0 at 30 draws; −0.35 per firm, significant at 10
draws), so a financial-industry shock is not supported as its explanation, though the evidence cannot rule it out. On Panel A every method gives positive effects for all three
connections (schedule 1.43–2.06, New York 0.74–1.19 % per day).
