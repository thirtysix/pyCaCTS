#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Stage per-cell-line CaCTS data for the dashboard's cell-line level, plus the model-type
code -> description map. One small JSON per line (lazy-fetched by the frontend) holding CaCTS scores +
expression in a fixed global TF order (names live once in tf_names.json), an index for the searchable
picker, and modeltype_desc.json. Data dir: PYCACTS_DEPMAP env var (default data/depmap/); see README."""
import os, sys, json
from pathlib import Path
import pandas as pd
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, grouping, score

DATA = HERE.parent / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))
EXPR = DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL = DEPMAP / "Model.csv"; TF = DATA / "CaCTS_merged_1671_TFs.txt"
D = HERE.parent / "dashboard" / "data"; LD = D / "lines"; LD.mkdir(parents=True, exist_ok=True)


def main():
    expr = io.load_expression(EXPR); model = io.load_model(MODEL); tfs = io.load_tf_universe(TF, "cacts")
    rep, _ = grouping.build_rep_matrix(expr, model, "line", tf_universe=tfs)   # TF x ACH (per-line log2TPM)
    S = score.cacts_score_matrix(rep)                                          # TF x ACH (per-line CaCTS)
    m = pd.read_csv(MODEL)
    name = {r.ModelID: (r.StrippedCellLineName if pd.notna(r.StrippedCellLineName) else r.ModelID) for r in m.itertuples()}
    subt = {r.ModelID: (r.OncotreeSubtype if pd.notna(r.OncotreeSubtype) else "") for r in m.itertuples()}
    order = sorted(rep.index)                                # fixed global TF order shared by every line
    (D / "tf_names.json").write_text(json.dumps(order, separators=(",", ":")))
    Sr, repr_ = S.reindex(order), rep.reindex(order)
    idx = []
    for ach in rep.columns:
        (LD / f"{ach}.json").write_text(json.dumps({
            "a": ach, "n": name.get(ach, ach), "s": subt.get(ach, ""),
            "c": [round(float(x), 3) for x in Sr[ach].values],   # CaCTS scores in global TF order (names in tf_names.json)
            "e": [round(float(x), 2) for x in repr_[ach].values],
        }, separators=(",", ":")))
        idx.append({"a": ach, "n": name.get(ach, ach), "s": subt.get(ach, "")})
    idx.sort(key=lambda r: r["n"].upper())
    (D / "lines_index.json").write_text(json.dumps(idx, separators=(",", ":")))
    # model-type code -> human-readable description (e.g. GB -> Glioblastoma)
    dd = {}
    for r in m.dropna(subset=["DepmapModelType"]).drop_duplicates("DepmapModelType").itertuples():
        dd[r.DepmapModelType] = (r.OncotreeSubtype if pd.notna(r.OncotreeSubtype)
                                 else (r.OncotreePrimaryDisease if pd.notna(getattr(r, "OncotreePrimaryDisease", None)) else r.DepmapModelType))
    (D / "modeltype_desc.json").write_text(json.dumps(dd, separators=(",", ":")))
    print(f"wrote {len(idx)} per-line JSON + lines_index.json + {len(dd)} modeltype descriptions")


if __name__ == "__main__":
    main()
