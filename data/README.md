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
| `CRISPRGeneEffect.csv` | CRISPR (Chronos) gene-effect, ~1,150 lines × ~18k genes | `stage_essentiality.py` (dashboard essentiality column) |

Place these under `$PYCACTS_DEPMAP`. The DepMap portal versions its releases; any recent release with these
files works (column names are stable). `CRISPRGeneEffect.csv` is only needed for the dashboard's per-group
essentiality column.

## TF annotation (dashboard gene-info)

`scripts/build_gene_info.py` builds the dashboard's per-TF dictionary (gene name, DNA-binding-domain family,
Entrez / Ensembl IDs) from the **Lambert et al. 2018** human-TF table. Point at it with `PYCACTS_TF_ANNOT`
(default `data/DatabaseExtract_v_1.01.csv`):

| File | What | Source |
| :-- | :-- | :-- |
| `DatabaseExtract_v_1.01.csv` | Lambert 2018 human-TF table (HGNC symbol, DBD family, EntrezGene description + ID, Ensembl ID) | humantfs.ccbr.utoronto.ca / Lambert et al., *Cell* 2018 |

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

## TCGA panel (dashboard TCGA data)

`scripts/stage_tcga.py` builds the dashboard's `data/tcga/` bundle (tumor type / molecular subtype / sample
type). It needs three downloads (paths via env vars; see the README worked example):

| Input | env var | What | Source |
| :-- | :-- | :-- | :-- |
| TCGA expression (TPM) | `PYCACTS_TCGA_EXPR` | Xena **Toil** `tcga_RSEM_gene_tpm` re-encoded to log2(TPM+1) by `convert_toil_tcga_tpm.py` (same units as DepMap; a literal 1 TPM floor) | toil-xena-hub S3 |
| sample → tumor type | `PYCACTS_TCGA_TYPES` | CaCTS `SuppTable1-34-TCGAID.txt` (33 types) | github.com/lawrenson-lab/CaCTS `files/` |
| sample → molecular subtype | `PYCACTS_TCGA_SUBTYPE` | UCSC Xena `TCGASubtype.20170308.tsv` (`Subtype_Selected`) | pancanatlas.xenahubs.net |
| patient → AJCC stage | `PYCACTS_TCGA_STAGE` | `TCGA_stage.tsv` (patient, major stage I/II/III/IV) from `build_tcga_stage.py`, for the within-cancer Stage axis | GDC clinical API |

`convert_toil_tcga_tpm.py` prepares `PYCACTS_TCGA_EXPR` from two Toil-hub downloads: the RSEM gene-TPM
matrix (`PYCACTS_TOIL_TPM`, log2(TPM+0.001), Ensembl-keyed) and the gencode.v23 gene probemap
(`PYCACTS_TOIL_MAP`, Ensembl → symbol). It maps to symbols, keeps the CaCTS TF universe, and re-encodes to
log2(TPM+1). The earlier batch-corrected `EB++AdjustPANCAN` matrix was replaced because it is not TPM (so a
1 TPM abundance floor could not be applied literally); Toil is uniformly re-quantified TPM.

The staged `data/tcga/*.tsv` are aggregate per-group statistics derived from the above; the raw downloads
are not committed. TCGA/Xena data are freely available; see the citations in the top-level README.

Both panels call an MTF **specific** when its empirical-null **FDR &lt; 0.10** (data-driven specificity gate,
replacing CaCTS's fixed top-5%-by-score cutoff) **and** its mean expression is **&ge; 1 TPM** (a light
abundance floor that keeps genuinely expressed lineage TFs, e.g. ovarian SOX17/WT1, while dropping
near-silent JSD artifacts). `pycacts/stats.py` computes the FDR; `filter.mtf_categories_fdr` does the call.
