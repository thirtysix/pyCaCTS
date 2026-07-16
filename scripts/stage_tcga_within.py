#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Stage WITHIN-cancer CaCTS for the dashboard's 'Within-cancer' tab: for each cancer, score CaCTS across
its OWN groups (reference set = that cancer's samples only), on two axes:

  subtype      : the cancer's molecular subtypes (UCSC Xena Subtype_Selected)
  tumornormal  : the cancer's [Tumor, Normal-tissue] samples (from the TCGA barcode)

(A stage axis would go here too, but the AJCC-stage clinical table isn't freely fetchable from Xena.)

Writes dashboard/data/tcga/within_<axis>_mtfs.tsv (cols: cancer, group, group_size, tf, rank, category,
cacts_score, group_expr_log2tpm) + within_manifest.json ({axis: {cancer: {group: n}}}). Each cancer's
samples are found via the patient barcode (first 12 chars) joined to the CaCTS primary sample map, so
metastatic / normal samples inherit their patient's cancer. Same inputs as stage_tcga.py."""
import os, sys, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, score, filter as cfilter

DATA = HERE.parent / "data"
EXPR = Path(os.environ.get("PYCACTS_TCGA_EXPR", DATA / "tcga" / "TCGA_pancan.geneExp.gz"))
SMAP = Path(os.environ.get("PYCACTS_TCGA_TYPES", DATA / "tcga" / "TCGA_sample_types.txt"))
SUBTYPE = Path(os.environ.get("PYCACTS_TCGA_SUBTYPE", DATA / "tcga" / "TCGASubtype.tsv.gz"))
TF = DATA / "CaCTS_merged_1671_TFs.txt"
OUT = HERE.parent / "dashboard" / "data" / "tcga"; OUT.mkdir(parents=True, exist_ok=True)
MIN_N = 10                                                      # min samples per within-cancer group


TOPN = 20                                                       # top TFs by specificity score kept per group


def score_groups(expr, labels):
    """labels: Series barcode -> within-group. Returns (rows, {group: size}) or (None, None) if <2 groups.
    Within one cancer the strict specific-MTF set is tiny (shared lineage), so keep each group's top-TOPN
    TFs by CaCTS score (most subtype/state-specific), annotated with the strict class."""
    lab = labels.dropna()
    lab = lab[[c for c in lab.index if c in expr.columns]]
    counts = lab.value_counts()
    keep = counts[counts >= MIN_N].index
    lab = lab[lab.isin(keep)]
    if lab.nunique() < 2:
        return None, None
    gsize = counts[keep].astype(int)
    rep = expr[list(lab.index)].T.groupby(lab).mean().T.dropna(how="all").fillna(0.0).clip(lower=0.0)
    S = score.cacts_score_matrix(rep)
    rows = []
    for g in rep.columns:
        cats = cfilter.mtf_categories(S, rep, g)
        spec = set(cats["specific"])
        expressed = list(spec | set(cats["non_specific"]))       # the top-5%-expressed TFs in this group
        ranked = S[g][expressed].sort_values()                   # most specific among the expressed (ascending)
        for rank, (tf, sc) in enumerate(ranked.head(TOPN).items(), 1):
            rows.append(dict(group=g, group_size=int(gsize.get(g, 0)), tf=tf, rank=rank,
                             category="specific" if tf in spec else "non_specific",
                             cacts_score=round(float(sc), 6), group_expr_log2tpm=round(float(rep.loc[tf, g]), 2)))
    return rows, gsize.to_dict()


def main():
    tfs = set(io.load_tf_universe(TF, "cacts"))
    print("loading TCGA expression ...")
    expr = pd.read_csv(EXPR, sep="\t", index_col=0)
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.loc[expr.index.intersection(tfs)]
    if expr.index.duplicated().any():
        expr = expr.groupby(level=0).mean()
    expr.columns = [str(c) for c in expr.columns]
    print(f"  {expr.shape[0]} TFs x {expr.shape[1]} samples")

    smap = pd.read_csv(SMAP, sep="\t")
    pat2cancer = dict(zip(smap["SampleId"].str[:12], smap["Cancer"]))     # patient -> cancer (from primaries)
    bc_cancer = pd.Series({c: pat2cancer.get(c[:12]) for c in expr.columns}).dropna()   # every sample -> cancer
    sub = pd.read_csv(SUBTYPE, sep="\t")
    bc_sub = dict(zip(sub["sampleID"].astype(str).str[:15], sub["Subtype_Selected"]))

    def clean_sub(c):                                          # subtype call, dropping the "NA"/".NA" non-calls
        v = bc_sub.get(c[:15])
        return None if (v is None or pd.isna(v) or v == "NA" or str(v).endswith(".NA")) else v

    axes = {"subtype": {}, "tumornormal": {}}
    manifest = {"subtype": {}, "tumornormal": {}}
    for cancer, samples in bc_cancer.groupby(bc_cancer).groups.items():
        cols = list(samples)
        # subtype axis: this cancer's molecular subtypes
        sublab = pd.Series({c: clean_sub(c) for c in cols})
        rows, gs = score_groups(expr, sublab)
        if rows:
            for r in rows: r["cancer"] = cancer
            axes["subtype"][cancer] = rows; manifest["subtype"][cancer] = gs
        # tumor vs normal axis
        tnlab = pd.Series({c: ("Normal" if c[13:15] == "11" else "Tumor") for c in cols})
        rows, gs = score_groups(expr, tnlab)
        if rows and set(gs) == {"Tumor", "Normal"}:
            for r in rows: r["cancer"] = cancer
            axes["tumornormal"][cancer] = rows; manifest["tumornormal"][cancer] = gs

    for axis in ("subtype", "tumornormal"):
        allrows = [r for rows in axes[axis].values() for r in rows]
        cols = ["cancer", "group", "group_size", "tf", "rank", "category", "cacts_score", "group_expr_log2tpm"]
        pd.DataFrame(allrows)[cols].to_csv(OUT / f"within_{axis}_mtfs.tsv", sep="\t", index=False)
        print(f"  {axis}: {len(manifest[axis])} cancers")
    (OUT / "within_manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
    print("wrote within_{subtype,tumornormal}_mtfs.tsv + within_manifest.json")


if __name__ == "__main__":
    main()
