# pyCaCTS

A fast, dependency-light **Python reimplementation of CaCTS** (Cancer Core Transcription factor
Specificity), applied to the **DepMap/CCLE** cell-line panel to nominate master transcription factors
(MTFs) at every level of the disease hierarchy, lineage, primary disease, subtype, model type, and each
individual cell line.

> **Credit / original method.** CaCTS is the work of the Lawrenson lab,
> **Reddy J, Fonseca MAS, Corona RI, *et al.*, "Predicting master transcription factors from pan-cancer
> expression data," *Sci. Adv.* 2021;7(48):eabf6123** (PMID 34818047; DOI 10.1126/sciadv.abf6123),
> code at <https://github.com/lawrenson-lab/CaCTS>. This repository is an **independent Python port** of
> their Jensen-Shannon-divergence specificity score, validated numerically against the original R, and
> extended to cell-line panels. All credit for the method is theirs; any errors in this port are mine.
> Their original R is GPL and is **not** redistributed here.

## Why
CaCTS identifies MTFs from **expression alone**: how *specifically* a transcription factor is expressed in
one group relative to the diversity of all groups (a Jensen-Shannon-divergence measure), combined with an
absolute-expression filter. It needs **no ChIP-seq / epigenomic data**. The original was run on TCGA
tumor types; pyCaCTS runs the same method on the DepMap/CCLE cell-line panel (~1,450 lines with
expression), so you can read out the candidate master regulators of any cancer group, down to a single
cell line.

## What it does
- Computes the CaCTS specificity score for every TF × group.
- Groupings supported (one engine): **OncotreeLineage**, **OncotreePrimaryDisease**, **OncotreeSubtype**,
  **DepmapModelType**, and **individual cell line**.
- Applies the CaCTS specific / non-specific MTF definitions (top-5% specificity ∩ top-5% expression;
  and the high-expression / low-specificity "non-specific" category).
- Adds an **empirical-null FDR** per TF per group as a non-arbitrary alternative to the top-5% cutoff.

## Validation & performance
pyCaCTS's specificity score is **numerically identical to the original CaCTS R** and dramatically faster.
Validated by scoring the same representative matrices with both `pycacts` (vectorized) and the original
`run_CaCTS_score` (looped), on CCLE-derived inputs:

| input (TFs × groups) | pyCaCTS | original R | speed-up | max &#124;pyCaCTS − R&#124; |
| :-- | --: | --: | --: | --: |
| 1,651 × 29 (≈ TCGA scale) | **1.4 ms** | 6,025 ms | **~4,300×** | 8e-16 (identical) |
| 1,651 × 400 (per-line) | **32 ms** | 103,433 ms | **~3,200×** | 1e-15 (identical) |

Reproduce: `python scripts/benchmark_vs_r.py` (needs R on PATH and the original `reference_R/JSD.R`, kept
locally and **not** redistributed). The speed-up comes from an exact O(TFs × groups) vectorization of the
per-query JSD loop.

## Explore the results
A static, no-backend **dashboard** lives in `dashboard/`, a browsable MTF atlas across the disease
hierarchy, a sortable per-group **TF-scores** table with empirical-null FDR, and a TF finder. Run locally
with `cd dashboard && python3 -m http.server 8000`, or deploy via the included GitHub Pages workflow
(Settings → Pages → Source → GitHub Actions). Regenerate its data with `scripts/stage_dashboard_data.py`
and `scripts/stage_line_data.py`.

## Layout
```
pycacts/        the package (score / grouping / io / filter)
scripts/        runners (run_divisions.py, stage_dashboard_data.py, stage_line_data.py,
                benchmark_vs_r.py, cacts_reference.R)
dashboard/      static results explorer (index.html + css/ + js/ + data/); GitHub-Pages-ready
data/           bundled TF list + data pointers (large inputs are NOT committed, see data/README.md)
results/        outputs (score matrices, MTF lists): not committed
```

## Data (not committed)
Point pyCaCTS at a downloaded DepMap release via the `PYCACTS_DEPMAP` environment variable (default
`data/depmap/`). See `data/README.md` for the exact files and where to get them. The CaCTS 1,671-TF
catalogue is bundled at `data/CaCTS_merged_1671_TFs.txt`.

## License
pyCaCTS is released under the **GNU General Public License v3.0-or-later** (`LICENSE`),
© 2026 Harlan Barker. It is an independent reimplementation written from the published CaCTS method, the
original CaCTS R (also GPL) is credited above and is **not** redistributed here. GPL is a deliberate choice
to keep pyCaCTS, and anything built from it, open, mirroring how the original was released.
