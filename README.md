# pyCaCTS

A fast, dependency-light **Python reimplementation of CaCTS** (Cancer Core Transcription factor
Specificity), applied to the **DepMap/CCLE** cell-line panel to nominate master transcription factors
(MTFs) at every level of the disease hierarchy, lineage, primary disease, subtype, model type, and each
individual cell line.

![pyCaCTS dashboard](docs/overview.png)

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

## Usage
pyCaCTS is a small importable package (numpy + pandas). Point it at a DepMap release (or any group-mean
expression matrix) and score every TF × group in one call.

```python
from pycacts import io, grouping, score, filter as cfilter

# load a DepMap release + the CaCTS TF universe
expr  = io.load_expression("depmap/OmicsExpressionProteinCodingGenesTPMLogp1.csv")  # genes x lines
model = io.load_model("depmap/Model.csv")                                           # indexed by ModelID
tfs   = io.load_tf_universe("data/CaCTS_merged_1671_TFs.txt")                       # 1,671 TFs

# per-group mean expression at the lineage level, then the CaCTS specificity score
rep, sizes = grouping.build_rep_matrix(expr, model, "lineage", tf_universe=tfs)     # TFs x lineages
scores = score.cacts_score_matrix(rep)          # TFs x lineages, lower = more specific (a few ms)

# master TFs of the Skin lineage: specific = top-5% score AND top-5% expression
mtfs = cfilter.mtf_categories(scores, rep, "Skin")
print(mtfs["specific"])                          # ['SOX10', 'LEF1', 'MITF']
print(score.rank_specific(scores, "Skin").head())
```

Any column in `Model.csv` is a valid grouping, `"OncotreeSubtype"`, `"DepmapModelType"`, or `"line"` (each
cell line its own group, ~1,450 groups scored in ~0.2 s):

```python
rep, _ = grouping.build_rep_matrix(expr, model, "line", tf_universe=tfs)
S = score.cacts_score_matrix(rep)
```

It is not tied to DepMap, score any **TFs × groups** mean-expression matrix (TCGA tumor types, GTEx
tissues, single-cell pseudobulk, …):

```python
import pandas as pd
from pycacts import score
rep = pd.read_csv("my_group_means.tsv", sep="\t", index_col=0)   # rows = TFs/genes, cols = groups
S = score.cacts_score_matrix(rep)                                # CaCTS JSD; lower = more group-specific
```

Runnable end-to-end example: `python examples/quickstart.py [LINEAGE]` (set `PYCACTS_DEPMAP`; see
`data/README.md`).

## Validation & performance
pyCaCTS's specificity score is **numerically identical to the original CaCTS R** and dramatically faster.
Scoring the same representative matrices with both `pycacts` (vectorized) and the original `run_CaCTS_score`
(looped), across the **whole DepMap/CCLE panel at all five hierarchy levels**:

| input (TFs × groups) | pyCaCTS | original R | speed-up | max &#124;pyCaCTS − R&#124; |
| :-- | --: | --: | --: | --: |
| 1,651 × 29 (lineage) | **1.3 ms** | 5.5 s | **~4,180×** | 8e-16 |
| 1,651 × 79 (primary disease) | **4.7 ms** | 15.1 s | **~3,240×** | 9e-16 |
| 1,651 × 191 (subtype) | **12.1 ms** | 40.2 s | **~3,320×** | 1e-15 |
| 1,651 × 192 (model type) | **11.9 ms** | 40.3 s | **~3,390×** | 1e-15 |
| 1,651 × 1,450 (per cell line) | **161 ms** | 10.8 min | **~4,020×** | 3e-15 |
| **full panel (all 5 levels)** | **0.19 s** | **12.5 min** | **~3,900×** | identical |

All scores match to floating-point precision (max |Δ| ≤ 3e-15). Reproduce: `python scripts/benchmark_vs_r.py`
(needs R on PATH and the original `reference_R/JSD.R`, kept locally and **not** redistributed).

### Where the speed-up comes from
`scripts/profile_speedup.py` times three Python implementations of the *identical* score on the same input,
isolating each source of the gain (measured at the lineage level, 1,651 TFs × 29 groups):

| step | time | gain | what it removes |
| :-- | --: | --: | :-- |
| original R (`run_CaCTS_score`) | 5,506 ms | 1× | baseline: O(T·G²), + a scores file written per query + `rbind`-in-loop |
| ↳ literal Python transcription | 462 ms | **~12×** | R's per-query file I/O and O(T²) `rbind` growth |
| ↳ vectorized over TFs (numpy) | 9.6 ms | **~48×** | the interpreted per-TF / per-query loops |
| ↳ closed-form (pyCaCTS) | 1.3 ms | **~7.5×** | re-summing over all groups per query |

The three gains are independent and multiply (≈ 12 × 48 × 7.5 ≈ 4,300×, matching the observed ~4,180×). The
first two are ~constant factors; the **closed-form** gain grows with the number of groups (~7× at 29 groups,
~38× at 191, **~286× at 1,450**), because CaCTS's JSD to a one-hot ideal has a closed form in a single value,
so the whole TF × group matrix is one O(T·G) elementwise pass instead of an O(T·G²) per-query loop. That is
what makes per-cell-line scoring tractable (44 s → 0.15 s).

## Explore the results
A static, no-backend **dashboard** lives in `dashboard/`, a browsable MTF atlas across the disease
hierarchy, a sortable per-group **TF-scores** table with empirical-null FDR, and a TF finder. Run locally
with `cd dashboard && python3 -m http.server 8000`, or deploy via the included GitHub Pages workflow
(Settings → Pages → Source → GitHub Actions). Regenerate its data with `scripts/stage_dashboard_data.py`
and `scripts/stage_line_data.py`.

![TF-scores table](docs/tf-scores.png)

## Layout
```
pycacts/        the package (score / grouping / io / filter)
scripts/        runners: run_divisions, stage_dashboard_data, stage_line_data, build_gene_info,
                stage_essentiality, benchmark_vs_r, profile_speedup, cacts_reference.R
examples/       quickstart.py (runnable end-to-end example)
dashboard/      static results explorer (index.html + css/ + js/ + data/); GitHub-Pages-ready
data/           bundled TF list + data pointers (large inputs are NOT committed, see data/README.md)
results/        outputs (score matrices, MTF lists, benchmark / profile CSVs): not committed
```

## Data (not committed)
Point pyCaCTS at a downloaded DepMap release via the `PYCACTS_DEPMAP` environment variable (default
`data/depmap/`). See `data/README.md` for the exact files and where to get them. The CaCTS 1,671-TF
catalogue is bundled at `data/CaCTS_merged_1671_TFs.txt`.

## Data sources & citations
pyCaCTS reimplements a published method and its dashboard ships aggregate statistics derived from public
resources. If you use it, please cite:

- **CaCTS (the method reimplemented here):** Reddy J, Fonseca MAS, Corona RI, *et al.* "Predicting master
  transcription factors from pan-cancer expression data." *Science Advances* 2021;7(48):eabf6123.
  DOI [10.1126/sciadv.abf6123](https://doi.org/10.1126/sciadv.abf6123) · code:
  <https://github.com/lawrenson-lab/CaCTS>.
- **DepMap / CCLE** (expression, CRISPR gene-effect, cell-line annotation): Broad Institute DepMap Public
  release, <https://depmap.org>. DepMap data are released under **CC BY 4.0**; the dashboard redistributes
  only per-group aggregates (means, ranks) derived from them; please cite the DepMap release and the
  CCLE/Achilles papers listed at depmap.org. Raw DepMap inputs are not committed (see `data/README.md`).
- **TF universe & families:** Lambert SA, *et al.* "The Human Transcription Factors." *Cell*
  2018;172(4):650–665. DOI 10.1016/j.cell.2018.01.029 (humantfs.ccbr.utoronto.ca); the CaCTS 1,671-TF list
  also draws on Saint-André V, *et al.*, *Genome Research* 2016.

A machine-readable `CITATION.cff` is included, so GitHub shows a "Cite this repository" button.

## License
pyCaCTS is released under the **GNU General Public License v3.0-or-later** (`LICENSE`),
© 2026 Harlan Barker. It is an independent reimplementation written from the published CaCTS method, the
original CaCTS R (also GPL) is credited above and is **not** redistributed here. GPL is a deliberate choice
to keep pyCaCTS, and anything built from it, open, mirroring how the original was released.
