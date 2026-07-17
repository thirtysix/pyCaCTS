# pyCaCTS: design notes & roadmap

## What this is
A fast, dependency-light Python reimplementation of **CaCTS** (Reddy *et al.*, *Sci. Adv.* 2021), applied
to the DepMap/CCLE cell-line panel. CaCTS nominates master transcription factors from **expression alone**
(no ChIP-seq): for each TF it scores how *specifically* the factor is expressed in a query group versus the
diversity of all groups, then combines that with an absolute-expression filter.

## Method (what we reimplement)
For each TF, using a representative (mean) expression profile per group:
1. **score (JSD)**: per (TF, group): a Jensen-Shannon-divergence specificity of the TF's normalized
   cross-group expression profile vs a one-hot "ideal" at that group (mirrors the original `JSD.R`).
   Lower score = more group-specific.
2. **empirical-null FDR** (the significance gate, `pycacts/stats.py`): model the group's non-specific
   (high-JSD) bulk as the null, take a left-tail p-value per TF from its score, Benjamini-Hochberg correct.
3. **call** (`filter.mtf_categories_fdr`): **specific MTFs** = empirical-null FDR < 0.10 ∩ mean expression
   ≥ 1 TPM; **non-specific MTFs** = high expression (top-5%) but not FDR-significant (candidate ubiquitous /
   multi-cancer regulators). This is the tool-wide call: it replaces CaCTS's original fixed top-5%-by-score
   ∩ top-5%-by-expression cutoffs with a data-driven significance threshold plus a light 1 TPM abundance
   floor, which keeps genuinely expressed lineage TFs the aggressive top-5%-expression gate dropped (e.g.
   ovarian SOX17 / WT1 / MECOM) while still excluding near-silent JSD artifacts. The original fixed-cutoff
   call is kept available as `filter.mtf_categories`.

Scores are validated to be numerically identical to the original R (max |Δ| < 1e-15 on identical input).

## Groupings (one engine, five resolutions)
- **OncotreeLineage** (29 groups): most robust.
- **OncotreePrimaryDisease** (79) and **OncotreeSubtype** (191).
- **DepmapModelType** (192): the finest Oncotree code.
- **Individual cell line** (~1,450 singleton groups): each line scored on its own profile vs the whole
  panel; noisier (single-sample), presented with that caveat.

## TF universe
- **Bundled:** CaCTS merged 1,671 (Lambert 2018 + Saint-André 2016), `data/CaCTS_merged_1671_TFs.txt`,
  direct method fidelity. 1,651 of these are expressed in the panel and scored.

## Status
- [x] **Engine**: `pycacts` package (io, grouping, vectorized JSD score, filter, stats).
- [x] **Validate**: numerically identical to the original CaCTS R on shared input; unit tests (`tests/`).
- [x] **Run**: all five resolutions scored; MTF lists + score matrices emitted.
- [x] **Benchmark**: pyCaCTS vs the original R on the same input (~3,000–4,300× faster).
- [x] **Dashboard**: static results explorer (MTF atlas, TF scores + FDR, Tumor-vs-model, Within-cancer, TF finder).
- [x] **Public release**: GPL-3.0, public GitHub + live GitHub Pages dashboard; two verified worked examples.
- [x] **TCGA arm**: full TCGA panel beside DepMap (Xena Toil RSEM-TPM), tumor-type / subtype / sample-type,
  a tumor-vs-cell-line-model comparison, and a within-cancer view (subtype / tumor-vs-normal / AJCC-stage).
- [x] **FDR + 1 TPM call**: the tool-wide MTF definition (see Method), replacing the fixed top-5% cutoffs.

## Possible future work (none committed; nice-to-have)
- **Permutation null** as a cross-check on the parametric empirical-null FDR: permute group labels, rebuild
  the score distribution, and compare the resulting FDR to the current closed-form null. Would validate the
  significance gate rather than assume its parametric form.
- **Per-cell-line CRISPR essentiality**: essentiality is staged per *group* (mean Chronos), not for the
  single-cell-line level (that column is blank there). Stage per-line Chronos so the cell-line resolution
  also carries a dependency column.
- **Bump the deprecated GitHub Actions**: `dashboard/.github/workflows/pages.yml` uses Node-20-deprecated
  action versions (checkout@v4 etc.); non-blocking warnings, bump when convenient.
- **AJCC-stage axis is broad-signal**: stage reflects progression/spread more than lineage identity, so the
  within-cancer Stage axis is weaker than subtype; early stages often have few stage-specific TFs. Could add
  a grade axis (for cancers that grade rather than stage) from the same GDC clinical source.

## Design principles
- Modular + importable (`from pycacts import score`), dependency-light (numpy/pandas), vectorized JSD.
- Deterministic; large inputs never committed (see `.gitignore`, `data/README.md`).
- Honest crediting of the original CaCTS method throughout; original GPL R never redistributed.

## Licensing
**Chosen: GPL-3.0-or-later** (`LICENSE`, © 2026 Harlan Barker). pyCaCTS is an independent reimplementation
from the method description / paper, it does not copy the original source (their `R/*.R` is kept locally,
gitignored, for reference/benchmarking only), so it carried no obligation to be GPL; GPL was chosen
deliberately to keep it and its derivatives open, mirroring the original CaCTS release. Source files carry
`SPDX-License-Identifier: GPL-3.0-or-later`. The original CaCTS is credited prominently throughout.

## Open decisions
- Per-line noise handling (rank + robustness flag; optional min-expression gate).
- Packaging surface: released GitHub-only so far; whether to publish to PyPI (and add a small CLI) is open.
