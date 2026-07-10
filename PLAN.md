# pyCaCTS — design notes & roadmap

## What this is
A fast, dependency-light Python reimplementation of **CaCTS** (Reddy *et al.*, *Sci. Adv.* 2021), applied
to the DepMap/CCLE cell-line panel. CaCTS nominates master transcription factors from **expression alone**
(no ChIP-seq): for each TF it scores how *specifically* the factor is expressed in a query group versus the
diversity of all groups, then combines that with an absolute-expression filter.

## Method (what we reimplement)
For each TF, using a representative (mean) expression profile per group:
1. **score (JSD)** — per (TF, group): a Jensen-Shannon-divergence specificity of the TF's normalized
   cross-group expression profile vs a one-hot "ideal" at that group (mirrors the original `JSD.R`).
   Lower score = more group-specific.
2. **filter** — top-5% by CaCTS score ∩ top-5% by mean expression = **specific MTFs**; high-expression but
   low-specificity (outside the top-5% score) = **non-specific MTFs** (candidate ubiquitous / multi-cancer
   regulators).
3. **empirical-null FDR** (pyCaCTS extension) — model the group's non-specific bulk as the null, take a
   left-tail p-value per TF, Benjamini-Hochberg correct. A data-driven alternative to the top-5% cutoff.

Scores are validated to be numerically identical to the original R (max |Δ| < 1e-15 on identical input).

## Groupings (one engine, five resolutions)
- **OncotreeLineage** (29 groups) — most robust.
- **OncotreePrimaryDisease** (79) and **OncotreeSubtype** (191).
- **DepmapModelType** (192) — the finest Oncotree code.
- **Individual cell line** (~1,450 singleton groups) — each line scored on its own profile vs the whole
  panel; noisier (single-sample), presented with that caveat.

## TF universe
- **Bundled:** CaCTS merged 1,671 (Lambert 2018 + Saint-André 2016), `data/CaCTS_merged_1671_TFs.txt` —
  direct method fidelity. 1,651 of these are expressed in the panel and scored.

## Status
- [x] **Engine** — `pycacts` package (io, grouping, vectorized JSD score, filter).
- [x] **Validate** — numerically identical to the original CaCTS R on shared input.
- [x] **Run** — all five resolutions scored; MTF lists + score matrices emitted.
- [x] **Benchmark** — pyCaCTS vs the original R on the same input (~3,000–4,300× faster).
- [x] **Dashboard** — static results explorer (MTF atlas, TF scores + FDR, TF finder).
- [ ] **Package + release** — packaging metadata, a small CLI, example, and public GitHub push.
- [ ] **Optional TCGA arm** — reproduce the published TCGA scores end-to-end as an extra validation.

## Design principles
- Modular + importable (`from pycacts import score`), dependency-light (numpy/pandas), vectorized JSD.
- Deterministic; large inputs never committed (see `.gitignore`, `data/README.md`).
- Honest crediting of the original CaCTS method throughout; original GPL R never redistributed.

## Licensing
**Chosen: GPL-3.0-or-later** (`LICENSE`, © 2026 Harlan Barker). pyCaCTS is an independent reimplementation
from the method description / paper — it does not copy the original source (their `R/*.R` is kept locally,
gitignored, for reference/benchmarking only), so it carried no obligation to be GPL; GPL was chosen
deliberately to keep it and its derivatives open, mirroring the original CaCTS release. Source files carry
`SPDX-License-Identifier: GPL-3.0-or-later`. The original CaCTS is credited prominently throughout.

## Open decisions
- Per-line noise handling (rank + robustness flag; optional min-expression gate).
- Whether to add the full TCGA validation/benchmark arm now or later.
- Packaging surface (CLI shape, PyPI or GitHub-only).
