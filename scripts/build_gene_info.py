#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Build dashboard/data/gene_info.json: a per-TF metadata dictionary the tables join against.

Fields per TF symbol: full gene name, TF family (DNA-binding domain), Entrez + Ensembl IDs (for
external links), and cross-group breadth (# of group-level calls where it is a *specific* MTF).

Name / family / IDs come from the Lambert et al. 2018 human-TF table (set its path via the
PYCACTS_TF_ANNOT env var, default data/DatabaseExtract_v_1.01.csv; see data/README.md). Breadth is
derived from the dashboard's own mtfs_<div>.tsv, so it needs stage_dashboard_data.py to have run first.
"""
import os, re, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
DASH = HERE.parent / "dashboard" / "data"
LAMBERT = Path(os.environ.get("PYCACTS_TF_ANNOT", DATA / "DatabaseExtract_v_1.01.csv"))
GROUP_DIVS = ["lineage", "disease", "subtype", "modeltype"]     # per-line breadth is computed client-side

_SRC = re.compile(r"\s*\[Source:.*$")                           # strip the "[Source:HGNC …]" suffix


def main():
    universe = json.loads((DASH / "tf_names.json").read_text())
    lam = pd.read_csv(LAMBERT).set_index("HGNC symbol")
    # cross-group breadth: number of (division, group) where the TF is a specific MTF
    breadth = {}
    for div in GROUP_DIVS:
        f = DASH / f"mtfs_{div}.tsv"
        if not f.exists():
            continue
        t = pd.read_csv(f, sep="\t")
        for tf in t.loc[t.category == "specific", "tf"]:
            breadth[tf] = breadth.get(tf, 0) + 1

    info = {}
    for tf in universe:
        rec = {"name": "", "family": "", "entrez": "", "ensembl": "", "breadth": int(breadth.get(tf, 0))}
        if tf in lam.index:
            r = lam.loc[tf]
            r = r.iloc[0] if isinstance(r, pd.DataFrame) else r
            desc = str(r.get("EntrezGene Description", "")).strip()
            rec["name"] = _SRC.sub("", desc) if desc and desc != "nan" else ""
            fam = str(r.get("DBD", "")).strip()
            rec["family"] = "" if fam in ("", "nan") else fam
            eid = r.get("EntrezGene ID", "")
            rec["entrez"] = "" if pd.isna(eid) else str(int(eid)) if str(eid).replace(".0", "").isdigit() else str(eid)
            ens = str(r.get("Ensembl ID", "")).strip()
            rec["ensembl"] = "" if ens in ("", "nan") else ens
        info[tf] = rec

    (DASH / "gene_info.json").write_text(json.dumps(info, separators=(",", ":")))
    n_named = sum(1 for v in info.values() if v["name"])
    n_fam = sum(1 for v in info.values() if v["family"])
    print(f"wrote gene_info.json: {len(info)} TFs, {n_named} with name, {n_fam} with family, "
          f"{sum(1 for v in info.values() if v['breadth'])} appear as a specific MTF somewhere")


if __name__ == "__main__":
    main()
