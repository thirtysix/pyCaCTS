# pyCaCTS dashboard

A static, no-backend results explorer for the pyCaCTS cell-line MTF atlas. Loads precomputed TSV/JSON
from `data/`; no build step.

## Tabs
- **Overview** — validation + benchmark + what to explore.
- **MTF atlas** — specific / non-specific MTFs for every group, at five resolutions of the DepMap disease
  hierarchy (lineage, primary disease, subtype, model type, and each individual cell line).
- **TF scores** — every TF's CaCTS specificity score for a chosen group, sortable, with MTF class and an
  empirical-null FDR you can filter by.
- **TF finder** — where a given TF is a specific master regulator.
- **About & methods** — method, MTF-calling rule, and credit to the original CaCTS.

## Run locally
A server is required (`fetch()` of local TSVs is blocked over `file://`):
- macOS / Linux: `python3 -m http.server 8000` then open http://localhost:8000
- Node: `npx serve -l 8000`

## Data
Regenerate `data/` from the analysis with:
- `scripts/stage_dashboard_data.py` — per-division score + expression matrices, MTF tables, manifest.
- `scripts/stage_line_data.py` — the per-cell-line JSONs, the shared `tf_names.json`, and the line index.

Divisions shipped: lineage (29 groups), primary disease (79), subtype (191), model type (192), plus the
per-cell-line level (~1,450 lines). All from DepMap/CCLE expression scored against the CaCTS 1,671-TF
catalogue.

## Deploy (GitHub Pages)
`.github/workflows/pages.yml` publishes this `dashboard/` directory on push to `main`. One manual step:
**Settings → Pages → Source → GitHub Actions**.
