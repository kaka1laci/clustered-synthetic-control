# S&P 500 inclusions, 1980 – Sep 1989 (GP&L §5.2): ours vs paper

224 inclusions of CRSP common stocks (248 additions in the window; Oct–Dec 1989 dropped). Day 0 = CRSP
effective date (229 of 241 are Thursdays); the literal Wednesday "announcement date" gives ~0.1% for every
estimator, so the paper's day 0 is the effective date. Loadings on days −250..−101.
Script: `replication/src/08_sp500_1980s.R`; results: `replication/output/sp500_1980s.rds`.

## Table 5, 1980–1989 row (announcement-day effect, %)

| | Diff-in-Means | Market | CAPM | FF3F | Gsynth |
|---|---|---|---|---|---|
| Ours | **3.25** | 3.36 (S&P) / 3.33 (CRSP VW) | 3.30 (S&P) / 3.28 (VW) / 3.28 (FF Mkt) | 3.21 | 3.26 |
| Paper | 3.27 | 3.25 | 3.15 | 3.05 | 3.06 |

Diff-in-Means matches; the model-based columns are 0.08–0.20 higher. The paper's point — model choice barely
matters for short-horizon effects with randomly timed events (Theorem 2.1) — holds: our estimators span
0.15 pp, the paper's 0.22 pp.

## Table 4 Panel A, betas of included firms (days −250..−100)

| | CAPM (S&P) | CAPM (CRSP VW) | FF3F Mkt | FF3F SMB | FF3F HML |
|---|---|---|---|---|---|
| Ours mean (sd) | 0.866 (0.457) | 1.020 (0.498) | 1.145 (0.489) | 0.536 (0.628) | −0.199 (0.954) |
| Paper mean (sd) | 0.961 (0.523) | | 1.108 (0.539) | 0.558 (0.604) | −0.148 (0.987) |

## Pre-inclusion drift (cumulative from day −100, %)

| Day | −80 | −50 | −15 | −1 | 0 | +20 |
|---|---|---|---|---|---|---|
| Market-adjusted (S&P 500) | 2.21 | 4.58 | 6.44 | 6.62 | 9.98 | 10.53 |
| Gsynth | 0.97 | 1.54 | 1.52 | 1.20 | 4.46 | 3.99 |

Gsynth removes about three quarters of the market-adjusted drift over days −100..−15 (1.52 vs 6.44), as the
paper reports. The paper's Figure 5 pools all decades, so these 1980s numbers are not directly comparable.
