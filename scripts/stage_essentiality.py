#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Stage per-group mean CRISPR (Chronos) essentiality for the dashboard's TF-scores table.

For each of the four group divisions (lineage, primary disease, subtype, model type) writes
ess_<div>.tsv = TFs × groups, each cell the mean Chronos across that group's cell lines with CRISPR data
(lower / more negative = stronger dependency; ~ -1 = median common-essential, 0 = neutral). The
single-cell-line level is not staged here. Reads DepMap CRISPRGeneEffect.csv from PYCACTS_DEPMAP.
"""
import os, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
DASH = HERE.parent / "dashboard" / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))
CRISPR = DEPMAP / "CRISPRGeneEffect.csv"
MODEL = DEPMAP / "Model.csv"
DIVS = {"lineage": "OncotreeLineage", "disease": "OncotreePrimaryDisease",
        "subtype": "OncotreeSubtype", "modeltype": "DepmapModelType"}


def main():
    universe = set(json.loads((DASH / "tf_names.json").read_text()))
    header = pd.read_csv(CRISPR, nrows=0).columns                     # "SYMBOL (ENTREZ)" columns + index
    sym = {c: c.split(" (")[0] for c in header}                       # colname -> gene symbol
    keep = [c for c in header if sym[c] in universe]                  # only our TF genes
    print(f"CRISPR TF columns matched: {len(keep)} of {len(universe)} TFs")
    eff = pd.read_csv(CRISPR, index_col=0, usecols=[header[0]] + keep)  # lines × TF genes
    eff.columns = [sym[c] for c in eff.columns]
    model = pd.read_csv(MODEL).set_index("ModelID")
    print(f"CRISPR lines: {eff.shape[0]}")

    for div, col in DIVS.items():
        grp = model.loc[[i for i in eff.index if i in model.index], col].dropna()
        E = eff.loc[grp.index]
        m = E.groupby(grp).mean().T                                   # TFs × groups (mean Chronos)
        m.index.name = "tf"
        m.round(3).to_csv(DASH / f"ess_{div}.tsv", sep="\t")
        print(f"  ess_{div}.tsv: {m.shape[0]} TFs × {m.shape[1]} groups")


if __name__ == "__main__":
    main()
