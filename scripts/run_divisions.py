#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Run pyCaCTS across several DepMap division schemes and emit MTFs for EVERY group.

Divisions: the Oncotree disease hierarchy (Lineage -> PrimaryDisease -> Subtype -> DepmapModelType) plus
cross-cutting biology (Sex, PrimaryOrMetastasis, GrowthPattern). For each division we score the full
TF x group matrix once, then write specific/non-specific MTFs for all groups (annotated with group size).

Data: set the DepMap release directory via the PYCACTS_DEPMAP env var (default: data/depmap/). See README.
"""
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, grouping, score, filter as cfilter

DATA = HERE.parent / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))    # DepMap release dir (download; see README)
EXPR = DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL = DEPMAP / "Model.csv"
TF_CACTS = DATA / "CaCTS_merged_1671_TFs.txt"                       # CaCTS TF universe (Lambert 2018 + Saint-Andre 2016)
OUT = HERE.parent / "results" / "divisions"; OUT.mkdir(parents=True, exist_ok=True)
DIVISIONS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "DepmapModelType",
             "Sex", "PrimaryOrMetastasis", "GrowthPattern"]


def mtfs_all_groups(S, rep, gsize, division):
    rows = []
    for g in rep.columns:
        cats = cfilter.mtf_categories(S, rep, g)
        rk = S[g].rank(method="first")
        for cat in ("specific", "non_specific"):
            for tf in cats[cat]:
                rows.append(dict(division=division, group=g, group_size=int(gsize.get(g, 0)),
                                 tf=tf, category=cat, jsd_rank=int(rk[tf]),
                                 cacts_score=round(float(S.loc[tf, g]), 6),
                                 group_expr_log2tpm=round(float(rep.loc[tf, g]), 2)))
    return pd.DataFrame(rows)


def main():
    t0 = time.time()
    print("loading expression / model / TFs ...")
    expr = io.load_expression(EXPR); model = io.load_model(MODEL); tfs = io.load_tf_universe(TF_CACTS, "cacts")
    print(f"  expr {expr.shape}; CaCTS TFs {len(tfs)}\n")
    manifest = []
    for div in DIVISIONS:
        rep, gsize = grouping.build_rep_matrix(expr, model, div, tf_universe=tfs, min_group_n=1)
        t = time.time(); S = score.cacts_score_matrix(rep); dt = time.time() - t
        tab = mtfs_all_groups(S, rep, gsize, div)
        tab.to_csv(OUT / f"mtfs_by_{div}.tsv", sep="\t", index=False)
        if rep.shape[1] <= 200:
            S.round(6).to_csv(OUT / f"scores_{div}.tsv", sep="\t")
        print(f"[{div:22s}] {rep.shape[0]} TFs x {rep.shape[1]:3d} groups ({dt*1000:5.0f} ms) "
              f"-> {tab.group.nunique()} groups w/ MTFs, {len(tab)} rows")
        manifest.append(dict(division=div, n_groups=rep.shape[1], n_mtf_rows=len(tab), score_ms=round(dt*1000, 1)))
    pd.DataFrame(manifest).to_csv(OUT / "divisions_manifest.tsv", sep="\t", index=False)
    print(f"\ntotal {time.time()-t0:.1f}s; wrote mtfs_by_*.tsv + score matrices + manifest to results/divisions/")


if __name__ == "__main__":
    main()
