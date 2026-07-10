#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Stage the pyCaCTS dashboard's data/ bundle: per-level CaCTS score & expression matrices,
MTF tables, plus manifest / tf-index / benchmark / meta JSON. Run with the shared analysis venv.

Organizational levels: OncotreeLineage (29 groups), OncotreePrimaryDisease (79),
OncotreeSubtype (191), DepmapModelType (192). Depends on results/divisions/mtfs_by_*.tsv (run_divisions.py).
"""
import os, sys, json, shutil
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, grouping, score

DATA = HERE.parent / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))    # DepMap release dir (download; see README)
EXPR = DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL = DEPMAP / "Model.csv"; TF = DATA / "CaCTS_merged_1671_TFs.txt"
D = HERE.parent / "dashboard" / "data"; D.mkdir(parents=True, exist_ok=True)
RDIV = HERE.parent / "results" / "divisions"
DIVS = {"lineage": "OncotreeLineage", "disease": "OncotreePrimaryDisease",
        "subtype": "OncotreeSubtype", "modeltype": "DepmapModelType"}


def main():
    expr = io.load_expression(EXPR); model = io.load_model(MODEL); tfs = io.load_tf_universe(TF, "cacts")
    manifest = {"divisions": {}}; tf_index = {}
    for short, col in DIVS.items():
        rep, gs = grouping.build_rep_matrix(expr, model, col, tf_universe=tfs)
        S = score.cacts_score_matrix(rep)
        rep.round(2).to_csv(D / f"expr_{short}.tsv", sep="\t")
        S.round(6).to_csv(D / f"scores_{short}.tsv", sep="\t")
        manifest["divisions"][short] = {"col": col, "n_groups": rep.shape[1],
                                        "groups": {g: int(gs.get(g, 0)) for g in rep.columns}}
        t = pd.read_csv(RDIV / f"mtfs_by_{col}.tsv", sep="\t"); t.to_csv(D / f"mtfs_{short}.tsv", sep="\t", index=False)
        for r in t[t.category == "specific"].itertuples():
            tf_index.setdefault(r.tf, []).append([short, r.group, int(r.jsd_rank)])
        print(f"{short}: {rep.shape}")
    (D / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
    (D / "tf_index.json").write_text(json.dumps(tf_index, separators=(",", ":")))
    (D / "meta.json").write_text(json.dumps({"n_tfs": 1671, "n_lines": "~1,450",
        "validated": "identical to the original R (max |Δ| < 1e-15)", "speedup": "~3,000–4,300×"}, separators=(",", ":")))
    (D / "benchmark.json").write_text(json.dumps({"rows": [
        {"input": "1,651 TFs × 29 groups (TCGA scale)", "py_ms": 1.4, "r_ms": 6025, "speedup": "~4,300×"},
        {"input": "1,651 TFs × 400 groups (per-line)", "py_ms": 32, "r_ms": 103433, "speedup": "~3,200×"}]}, separators=(",", ":")))
    print("staged", len(list(D.iterdir())), "files to", D)


if __name__ == "__main__":
    main()
