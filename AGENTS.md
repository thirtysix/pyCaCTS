# AGENTS.md

Orientation for AI coding agents (and humans) arriving at this repo: what it is, how to run it, and the conventions worth knowing before making changes.

## What this is
pyCaCTS: a fast, dependency-light Python reimplementation of CaCTS (Cancer Core Transcription factor Specificity), the Jensen-Shannon-divergence specificity score of Reddy, Fonseca, Corona et al. (Sci. Adv. 2021). It scores every transcription factor against every group at five nested levels of the Oncotree disease hierarchy on the DepMap/CCLE cell-line panel, and ships a static results dashboard. It is an independent port validated numerically against the original R; the original CaCTS R is credited but not redistributed here.

## Stack & layout
- **Package**: `pycacts/` (`score`, `grouping`, `io`, `filter`, `stats`); numpy + pandas only.
- **Scripts**: `scripts/` (division runners, dashboard/line/TCGA staging, `benchmark_vs_r.py`, `profile_speedup.py`).
- **Examples**: `examples/quickstart.py` (DepMap) and `examples/tcga.py` (TCGA), runnable end to end.
- **Dashboard**: `dashboard/` (static `index.html` + `css/` + `js/` + `data/`), GitHub-Pages-ready.
- **Tests**: `tests/` (`test_score.py`, `test_stats.py`).

## Run, test, lint
```bash
pytest                                     # tests import pycacts from the repo root
python examples/quickstart.py [LINEAGE]    # needs PYCACTS_DEPMAP set; see data/README.md
cd dashboard && python3 -m http.server 8000
```
No build step and no `pyproject.toml`; run from the repo root so `import pycacts` resolves.

## Conventions
- The specificity score is lower = more group-specific (JSD to a one-hot ideal, computed in closed form as one O(T*G) pass).
- A "specific" master TF is called by empirical-null FDR < 0.10 AND mean >= 1 TPM (`filter.mtf_categories_fdr`); the high-expression / low-specificity "non-specific" category is retained.
- Licensed GPL-3.0-or-later, a deliberate choice to keep the port and derivatives open.

## Gotchas
- Large inputs (DepMap release, TCGA matrices) are NOT committed. Point at a downloaded DepMap release via the `PYCACTS_DEPMAP` environment variable (default `data/depmap/`); see `data/README.md`. The bundled TF universe is `data/CaCTS_merged_1671_TFs.txt`; `results/` is not committed.
- `scripts/benchmark_vs_r.py` needs R on PATH and the original `reference_R/JSD.R`, which is kept locally and not redistributed.
