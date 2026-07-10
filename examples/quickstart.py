#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""pyCaCTS quickstart: score master TFs across DepMap and read out one group's MTFs.

Set PYCACTS_DEPMAP to a DepMap release directory holding
OmicsExpressionProteinCodingGenesTPMLogp1.csv and Model.csv (see data/README.md).

    python examples/quickstart.py            # master TFs of the Skin lineage
    python examples/quickstart.py Lung       # any DepMap OncotreeLineage
"""
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pycacts import io, grouping, score, filter as cfilter

DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", ROOT / "data" / "depmap"))
TF = ROOT / "data" / "CaCTS_merged_1671_TFs.txt"
GROUP = sys.argv[1] if len(sys.argv) > 1 else "Skin"

# 1. load a DepMap release + the CaCTS TF universe
expr = io.load_expression(DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv")   # genes x lines
model = io.load_model(DEPMAP / "Model.csv")                                           # indexed by ModelID
tfs = io.load_tf_universe(TF, "cacts")                                                # 1,671 TFs

# 2. per-group mean expression at the lineage level, then the CaCTS specificity score
rep, sizes = grouping.build_rep_matrix(expr, model, "lineage", tf_universe=tfs)       # TFs x lineages
scores = score.cacts_score_matrix(rep)                       # TFs x lineages; lower = more specific (a few ms)

# 3. read out one group
print(f"{GROUP}: {int(sizes.get(GROUP, 0))} cell lines\n")
print("most group-specific TFs (by CaCTS score):")
print(score.rank_specific(scores, GROUP).head(10).to_string(index=False))

mtfs = cfilter.mtf_categories(scores, rep, GROUP)            # specific = top-5% score AND top-5% expression
print(f"\nspecific MTFs ({len(mtfs['specific'])}): {', '.join(mtfs['specific'])}")
print(f"non-specific MTFs ({len(mtfs['non_specific'])}): {', '.join(mtfs['non_specific'][:8])} ...")
