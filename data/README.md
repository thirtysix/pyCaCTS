# Data pointers (large inputs are NOT committed)

pyCaCTS reads a DepMap/CCLE release from a directory you supply. Point at it with the `PYCACTS_DEPMAP`
environment variable (default `data/depmap/`):

```bash
export PYCACTS_DEPMAP=/path/to/depmap_release
```

## Required inputs (download from the DepMap portal, https://depmap.org/portal/data_page/)

| File | What | Used by |
| :-- | :-- | :-- |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | protein-coding expression, log2(TPM+1), ~1,450 lines × ~19k genes | all runners |
| `Model.csv` | cell-line annotation (`OncotreeLineage`, `OncotreePrimaryDisease`, `OncotreeSubtype`, `DepmapModelType`, names) | all runners |

Place both under `$PYCACTS_DEPMAP`. The DepMap portal versions its releases; any recent release with the
above two files works (column names are stable).

## Bundled in the repo

| File | What | Source |
| :-- | :-- | :-- |
| `CaCTS_merged_1671_TFs.txt` | the CaCTS TF universe (Lambert 2018 + Saint-André 2016) | github.com/lawrenson-lab/CaCTS `files/` |

## Optional (validation / benchmark only)

| Input | What | Source |
| :-- | :-- | :-- |
| `reference_R/JSD.R` | the original CaCTS R scoring code, kept **locally only**, GPL, not redistributed | github.com/lawrenson-lab/CaCTS |
| TCGA pan-cancer RNA-seq | to reproduce the published TCGA scores | github.com/lawrenson-lab/CaCTS `RNA.TCGA.pancancer.url.txt` |

The benchmark (`scripts/benchmark_vs_r.py`) needs R on `PATH` and a local `reference_R/JSD.R`; it is the
only path that touches the original GPL source, which is never committed to this repo.
