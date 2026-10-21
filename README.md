# JAK-STAT × MDD imaging-transcriptomics: a multi-modal evaluation with positive controls

Code, derived results, and figures for the manuscript:

> **The cortical principal gradient accounts for reported transcriptomic–imaging convergence of JAK-STAT signaling in major depressive disorder: a multi-modal evaluation with positive controls**

This repository accompanies a *methods-first* re-evaluation of imaging-transcriptomic associations between JAK-STAT signaling and major depressive disorder (MDD) cortical morphometry. The headline result is **not** a positive association: the observed spatial correspondence (ρ = 0.315) is **fully absorbed by the cortical principal gradient** (PC1–PC5 partial Spearman ρ = 0.005, BF₀₁ = 16.4), and is not specific to the JAK-STAT gene set under matched random-set nulls (empirical *P* = 0.985). Genetic (cis-MR at *IL6R* rs2228145) and postmortem (GSE53987) layers are reported alongside **positive controls** (schizophrenia, same cohort) and explicit **minimum detectable effect (MDE)** bounds so that null results are interpretable.

## What is in here

```
repo_jakstat_mdd/
├── README.md                 ← this file
├── LICENSE                   ← MIT
├── requirements.txt          ← Python 3.10+ dependencies
├── .gitignore
├── scripts/                  ← all analysis scripts (run in the order below)
│   ├── gradient_confound_check.py       # PC1–PC5 partial Spearman (the key analysis)
│   ├── sign_check.py                    # per-gene / per-parcel sign diagnostics
│   ├── test4_null_distributions.py      # 3× matched random gene-set nulls (10,000 draws)
│   ├── test4b_spatial_geneset_null.py   # spatial gene-set null (mean|rho|, grid top-18)
│   ├── test2_celltype_spatial.py        # cell-type spatial panel (brainrnaseq.org)
│   ├── run_il6r_mr.py                   # cis-MR: rs2228145 Asp358Ala → global CT
│   ├── run_il6r_panel.py                # regional CT + total surface area panel
│   ├── run_de_gse53987.py               # postmortem DE (MDD / SCZ / BD), OLS + covariates
│   ├── gse53987_addendum.py             # multi-probe rule, best-probe flag, FDR13/FDR39
│   ├── compute_bayes_factors.py         # JZS / BIC / Wakefield ABF battery
│   └── make_figures.py                  # Fig 1–3 (300 DPI, PNG/PDF/TIFF)
├── results/                  ← derived tables (all numbers quoted in the manuscript)
│   ├── bayes_results.json
│   ├── de_bf_all.csv
│   ├── gse53987_de_results.tsv / .json
│   ├── il6r_mr_results.json / .tsv
│   ├── il6r_mr_panel_results.json
│   ├── test2_celltype_spatial.csv / .md
│   ├── sign_check/           (partial_control_PC1.csv, sign_check_13genes.csv, sign_check_per_parcel.csv)
│   └── test4/                (matched-null CSVs, per-gene grid, Table-1 reconciliation)
├── figures/                  ← Fig 1–3 in PNG + PDF + TIFF (300 DPI) + fig1C underlying data
└── docs/                     ← full manuscript draft, submission strategy, reference list
```

## Key numbers (for cross-checking)

| Layer | Quantity | Result |
|---|---|---|
| Spatial (raw) | ρ(d_MDD, expression_JAK-STAT) | 0.315 (BF₁₀ > 10⁶) |
| Spatial (PC1 only) | partial ρ | 0.179 (BF₁₀ ≈ 14) |
| **Spatial (PC1–PC5)** | **partial ρ** | **0.005 (BF₀₁ = 16.4; 98.5% attenuation)** |
| Matched gene-set null | observed count = 3 | empirical *P* = 0.985 |
| cis-MR (*IL6R* rs2228145) | Wald ratio θ | −0.0115 SD/allele (95% CI −0.0267, +0.0036), *P* = 0.135 |
| MR MDE | 80% power floor | 0.0224 SD/allele (observed effect ≈ 53% of MDE) |
| Postmortem MDD (39 tests) | FDR hits | 0/39; BF₀₁ median 4.8 (0.55–5.61) |
| Postmortem SCZ positive control (39 tests) | FDR hits | 3/39 (IL6ST/PFC, STAT3/PFC, JAK1/HIP); BF₁₀ up to 33.4 |

The full numeric trace is in `results/` and the manuscript (`docs/路径1_全文初稿.md`, English IMRaD).

## Reproduction

```bash
# 1. Environment
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 2. Point the scripts at the raw inputs (see "Data sources" below).
#    Each script reads its inputs from a path constant at the top of the file;
#    edit those constants to match your local layout.

# 3. Run in this order (dependencies are light; each script is standalone):
python scripts/gradient_confound_check.py       # → PC1/PC1–5 partial correlations
python scripts/sign_check.py                    # → sign_check/*.csv
python scripts/test4_null_distributions.py      # → matched null distributions (.npy + .csv)
python scripts/test4b_spatial_geneset_null.py   # → spatial gene-set null
python scripts/test2_celltype_spatial.py        # → cell-type panel
python scripts/run_il6r_mr.py                   # → MR Wald ratio
python scripts/run_il6r_panel.py                # → regional MR panel
python scripts/run_de_gse53987.py               # → postmortem DE (414 rows)
python scripts/gse53987_addendum.py             # → multi-probe / FDR addendum
python scripts/compute_bayes_factors.py         # → bayes_results.json, de_bf_all.csv
python scripts/make_figures.py                  # → figures/Fig{1,2,3}_*.{png,pdf,tiff}
```

Random seeds are fixed inside the null-distribution scripts (`numpy.random.default_rng(20240517)` where applicable). The matched random gene-set procedure matches on mean expression, expression variance, and gene length; the exact matching rule is documented in the manuscript Methods and in `docs/路径1_全文初稿.md`.

## Data sources (not redistributed here)

| Input | Source | Used for |
|---|---|---|
| AHBA regional expression (Schaefer-400) | `expr_schaefer400_allgenes.csv` (derived from Allen Human Brain Atlas) | spatial layer |
| ENIGMA MDD Cohen's *d* (cortical thickness) | Schmaal et al. *Mol Psychiatry* 2017;22(6):900–909 (GCST90271790–14) | spatial layer |
| Global mean cortical thickness GWAS | Makowski et al. (GCST90271815) | cis-MR outcome |
| rs2228145 MDD GWAS | PGC (Wray et al. *Nat Genet* 2018) | cis-MR outcome |
| GSE53987 series matrix | NCBI GEO | postmortem DE (MDD/SCZ/BD) |
| Cell-type marker panels | brainrnaseq.org (Zhang et al. *Neuron* 2016;89(1):37–53; Zhang et al. *J Neurosci* 2014;34(36):11929–11947) | cell-source layer |

AHBA-derived expression data and ENIGMA summary statistics are third-party; please obtain and cite them from the original sources before running the scripts.

## Statistical conventions used in the scripts

- **Partial correlation**: rank-transform each variable (x, y, and each column of Z) → regress out Z → Pearson correlation of residuals. This is *not* `spearmanr(residuals)`. The implementation is in `gradient_confound_check.py:partial_spearman()` and is reproduced verbatim in `make_figures.py`.
- **Multiple comparisons**: BH-FDR applied separately within the 13-gene JAK-STAT panel (FDR13) and across the full 39 tests (FDR39 = 13 genes × 3 regions).
- **Multi-probe rule** (microarray): pre-registered as "highest mean expression across all samples"; flagged in `is_best_probe`.
- **Bayes factors**: JZS *t*-test (r = 0.707, numerical integration), BIC approximation (Wagenmakers 2007), and Wakefield ABF for MR (prior SD W swept 0.01–0.10). MR BFs are reported as "inconclusive" unless BF crosses 3; see manuscript for the reporting discipline.
- **MDE** is always reported alongside null MR/DE results (two-sided α = 0.05, 80% power).

## Citation

If you use this code or these derived results, please cite the manuscript (see `docs/`) and the primary data sources listed above.

```bibtex
@software{jakstat_mdd_gradient_2026,
  author    = {[作者姓名]},
  title     = {JAK-STAT x MDD imaging-transcriptomics: analysis code and derived results},
  year      = {2026},
  version   = {v1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.[回填]},
  url       = {https://github.com/[用户名]/jakstat-mdd-gradient}
}
```

Release/DOI procedure: see `RELEASE.md`.

## License

MIT — see `LICENSE`.
