# Clustered synthetic control for financial event studies

Code, tests and results for two pieces of work on synthetic control (SC) in financial event studies:

1. **A replication of Goldsmith-Pinkham & Lyu, *Causal Inference in Financial Event Studies*** ([arXiv 2511.15123](https://arxiv.org/abs/2511.15123)).
   The authors have not released code, so the estimators are rebuilt from public data, CRSP and R packages.
   It covers Table 2 (the Geithner nomination: Panel A with bank controls, Panel B with all CRSP common stocks), the Table 1 simulation, and S&P 500 inclusions in 1980 – Sep 1989.
2. **Donor selection by clustering for SC on daily stock returns**, tested on the same data:
   - a **ClusterSC** baseline (Rho et al., AISTATS 2025, [arXiv 2503.21629](https://arxiv.org/abs/2503.21629)): k-means on donor features, then SC within the target's cluster;
   - a **joint clustering + SC prototype** that chooses donor groups and SC weights in one mixed-integer quadratic program, solved with Gurobi. Approach 1 (做法 1) picks a group of m donors and its weights for each target. Approach 2 (做法 2) shares one partition of the donors across all targets and alternates the joint solve with a centre update.

This is research code. The results are prototype evidence, not a published method.

## Results

| Topic | Result | Details |
|---|---|---|
| Table 2, Panel A | 22 of 24 point estimates match the paper to three decimals | [`table2_panelA_comparison.md`](replication/output/table2_panelA_comparison.md) |
| Table 2, Panel B | The paper's main finding holds at 5%: SC, SDID and Gsynth are insignificant for every connection, Average and DID are significant | [`table2_panelB_complete_comparison.md`](replication/output/table2_panelB_complete_comparison.md) |
| Table 1 simulation | Same qualitative pattern; the timing-selection bias is about half the paper's | [`table1_comparison.md`](replication/output/table1_comparison.md) |
| S&P 500 inclusions, 1980s | Diff-in-Means 3.25% vs the paper's 3.27% | [`sp500_1980s_comparison.md`](replication/output/sp500_1980s_comparison.md) |
| ClusterSC on Panel B placebos | Clustering never lowers placebo error significantly and is worse in two comparisons; financial placebo cohorts share a bias of about −0.25 pp/day under every estimator | [`09_clustersc_panelB.R`](replication/src/09_clustersc_panelB.R), [`10_clustersc_perfirm.R`](replication/src/10_clustersc_perfirm.R) |
| Joint clustering + SC | Neither approach beats SC or ClusterSC consistently, in a simulation with group-specific factors or on Geithner placebos; pre-screening donors by correlation alone lowers placebo estimates by 0.06–0.37 pp/day | [`jointsc_sim_summary.csv`](replication/output/jointsc_sim_summary.csv), [`jointsc_geithner_summary.csv`](replication/output/jointsc_geithner_summary.csv) |
| Panel A and Gsynth (Phase 5) | With bank-only donors (Panel A), clustering and joint selection again do not lower placebo error; Gsynth is never significantly worse than the paired methods in the simulation but matches SC on both Geithner panels; the financial-firm placebo bias does not go away on Panel A | [`jointsc_geithnerA_summary.csv`](replication/output/jointsc_geithnerA_summary.csv), [`gsynth_fit.R`](replication/src/gsynth_fit.R) |

For the clustering work the result tables are the CSV files in `replication/output/`; written summaries and run logs are not published.

## Layout

```
replication/
├── README.md   how to run each step, implementation details, data sources
├── src/        R scripts
│   ├── 01–08   replication: Geithner panels, Table 2 estimators, Table 1 simulation, S&P 500
│   ├── 09–10   ClusterSC baseline (cohort-average and per-firm targets)
│   ├── 11–12   joint clustering + SC and Gsynth: simulation, Geithner Panels A and B
│   ├── sc_qp.R, clustersc.R, jointsc.R, gsynth_fit.R   SC solver, ClusterSC, joint clustering + SC, Gsynth
│   ├── check_*.R                         diagnostics
│   └── paper_table{1,2}.csv              the paper's published values
├── tests/      correctness tests
└── output/     result tables (.csv, .md), diagnostics/
```

## Data (not included)

- **CRSP** (via WRDS) is licensed and is not distributed here. `replication/src/04_pull_crsp.R` downloads it with your own WRDS account.
- The **Acemoglu et al. (2016)** replication data and the **Fama-French** daily factors are public; download links are in [`replication/README.md`](replication/README.md#data-sources).
- Intermediate panels and per-run outputs (`.rds`) are not tracked; the scripts rebuild them.

## Requirements

- R 4.6.1 with `haven`, `quadprog`, `DBI`, `RPostgres`, `synthdid` (GitHub, commit 70c1ce3), and **gsynth 1.2.1** + `lfe` (the version the paper cites; build notes in [`replication/README.md`](replication/README.md))
- For the joint clustering + SC prototype only: Gurobi 13.0.3 with a licence, and its R package `gurobi` 13.0-3 + `slam`

## Running

Every step's command is in [`replication/README.md`](replication/README.md); run them from the repository root. Tests:

```bash
Rscript --vanilla replication/tests/test_sc_qp.R replication/src replication/output/geithner_panel.rds   # needs the Panel A panel (01)
Rscript --vanilla replication/tests/test_clustersc.R replication/src
Rscript --vanilla replication/tests/test_jointsc.R replication/src                                        # needs Gurobi
Rscript --vanilla replication/tests/test_jointsc_combine.R replication/src <scratch_dir>
```

## License

The code is released under the [MIT License](LICENSE). CRSP and other third-party data are not included and keep their own terms.
