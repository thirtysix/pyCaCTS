#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Stage TCGA data for the dashboard, mirroring the DepMap output format under dashboard/data/tcga/.

Scores CaCTS at three levels: tumor type (33, from the CaCTS sample map), molecular subtype (from the
UCSC Xena TCGASubtype `Subtype_Selected` calls), and sample type (primary / metastatic / normal / …,
derived from the TCGA barcode). No essentiality (no CRISPR for tumors) and no per-sample level.

Inputs (download; see the README 'Worked example'), paths via env vars:
  PYCACTS_TCGA_EXPR     TCGA expression as log2(TPM+1), gene-symbol-keyed (genes x samples): the Xena Toil
                        RSEM-TPM matrix re-encoded by convert_toil_tcga_tpm.py (same units as DepMap)
  PYCACTS_TCGA_TYPES    CaCTS SuppTable1-34-TCGAID.txt  (columns Cancer, Category, SampleId)
  PYCACTS_TCGA_SUBTYPE  Xena TCGASubtype.20170308.tsv   (columns sampleID, …, Subtype_Selected)
"""
import os, sys, json
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, score, filter as cfilter, stats

DATA = HERE.parent / "data"
EXPR = Path(os.environ.get("PYCACTS_TCGA_EXPR", DATA / "tcga" / "TCGA_toil_tpm_log2p1.tsv.gz"))
SMAP = Path(os.environ.get("PYCACTS_TCGA_TYPES", DATA / "tcga" / "TCGA_sample_types.txt"))
SUBTYPE = Path(os.environ.get("PYCACTS_TCGA_SUBTYPE", DATA / "tcga" / "TCGASubtype.tsv.gz"))
TF = DATA / "CaCTS_merged_1671_TFs.txt"
OUT = HERE.parent / "dashboard" / "data" / "tcga"; OUT.mkdir(parents=True, exist_ok=True)
MIN_N = 5                                                        # drop groups with fewer samples
FDR_MAX = 0.10                                                   # specific = empirical-null FDR <= this ...
EXPR_FLOOR = 1.0                                                 # ... AND mean expr >= 1 TPM (log2(TPM+1) units)

# TCGA barcode sample-type codes (char 14-15) -> label
SAMPLE_TYPE = {"01": "Primary Solid Tumor", "02": "Recurrent Solid Tumor", "03": "Primary Blood-Derived Cancer",
               "05": "Additional New Primary", "06": "Metastatic", "07": "Additional Metastatic",
               "11": "Solid Tissue Normal"}

# TCGA tumor-type code -> full study name (for the dashboard's expanded label, like DepMap's model type)
TCGA_TYPE_NAMES = {
    "ACC": "Adrenocortical Carcinoma", "BLCA": "Bladder Urothelial Carcinoma", "BRCA": "Breast Invasive Carcinoma",
    "CESC": "Cervical & Endocervical Carcinoma", "CHOL": "Cholangiocarcinoma", "COAD": "Colon Adenocarcinoma",
    "DLBC": "Diffuse Large B-cell Lymphoma", "ESCA": "Esophageal Carcinoma", "GBM": "Glioblastoma Multiforme",
    "HNSC": "Head & Neck Squamous Cell Carcinoma", "KICH": "Kidney Chromophobe",
    "KIRC": "Kidney Renal Clear Cell Carcinoma", "KIRP": "Kidney Renal Papillary Cell Carcinoma",
    "LAML": "Acute Myeloid Leukemia", "LGG": "Brain Lower-Grade Glioma", "LIHC": "Liver Hepatocellular Carcinoma",
    "LUAD": "Lung Adenocarcinoma", "LUSC": "Lung Squamous Cell Carcinoma", "MESO": "Mesothelioma",
    "OV": "Ovarian Serous Cystadenocarcinoma", "PAAD": "Pancreatic Adenocarcinoma",
    "PCPG": "Pheochromocytoma & Paraganglioma", "PRAD": "Prostate Adenocarcinoma", "READ": "Rectum Adenocarcinoma",
    "SARC": "Sarcoma", "SKCM": "Skin Cutaneous Melanoma", "STAD": "Stomach Adenocarcinoma",
    "TGCT": "Testicular Germ Cell Tumors", "THCA": "Thyroid Carcinoma", "THYM": "Thymoma",
    "UCEC": "Uterine Corpus Endometrial Carcinoma", "UCS": "Uterine Carcinosarcoma", "UVM": "Uveal Melanoma"}


def build_rep(expr, labels):
    """labels: Series sample-barcode -> group label (NaN = drop). Returns (TFs x groups mean, group_size)."""
    lab = labels.dropna()
    lab = lab[[c for c in lab.index if c in expr.columns]]
    counts = lab.value_counts()
    keep = counts[counts >= MIN_N].index
    lab = lab[lab.isin(keep)]
    rep = expr[list(lab.index)].T.groupby(lab).mean().T          # TFs x groups
    return rep, counts[keep].astype(int)


def score_and_write(expr, labels, div, col, manifest):
    rep, gsize = build_rep(expr, labels)
    # Toil expression is log2(TPM+1) >= 0; keep the fill for any all-NaN row from an empty group intersection.
    rep = rep.dropna(how="all").fillna(0.0)
    S = score.cacts_score_matrix(rep)
    S.round(4).to_csv(OUT / f"scores_{div}.tsv", sep="\t")
    rep.round(2).to_csv(OUT / f"expr_{div}.tsv", sep="\t")
    rows = []
    for g in rep.columns:
        f = stats.empirical_fdr_log10(S[g])                     # log10(empirical-null FDR) per TF
        cats = cfilter.mtf_categories_fdr(S, rep, g, fdr_log10=f, fdr_max=FDR_MAX, expr_floor=EXPR_FLOOR)
        rk = S[g].rank(method="first")
        for cat in ("specific", "non_specific"):
            for tf in cats[cat]:
                rows.append(dict(division=col, group=g, group_size=int(gsize.get(g, 0)), tf=tf, category=cat,
                                 jsd_rank=int(rk[tf]), cacts_score=round(float(S.loc[tf, g]), 6),
                                 group_expr_log2tpm=round(float(rep.loc[tf, g]), 2),
                                 fdr_log10=round(float(f[tf]), 3)))
    pd.DataFrame(rows).to_csv(OUT / f"mtfs_{div}.tsv", sep="\t", index=False)
    manifest["divisions"][div] = {"col": col, "n_groups": rep.shape[1], "groups": gsize.to_dict()}
    print(f"  {div:10s}: {rep.shape[1]} groups")


def main():
    tfs = set(io.load_tf_universe(TF, "cacts"))
    print("loading TCGA expression ...")
    expr = pd.read_csv(EXPR, sep="\t", index_col=0)              # genes x samples (symbol-keyed, log2)
    expr.index = expr.index.astype(str).str.upper()
    expr = expr.loc[expr.index.intersection(tfs)]
    if expr.index.duplicated().any():                           # collapse any duplicate gene symbols
        expr = expr.groupby(level=0).mean()
    expr.columns = [str(c) for c in expr.columns]
    print(f"  {expr.shape[0]} TFs x {expr.shape[1]} samples")

    manifest = {"divisions": {}}
    # 1. tumor type: CaCTS 33-type map, matched on the 15-char sample barcode
    smap = pd.read_csv(SMAP, sep="\t")
    bc2type = dict(zip(smap["SampleId"].str[:15], smap["Cancer"]))
    types = pd.Series({c: bc2type.get(c[:15]) for c in expr.columns})
    score_and_write(expr, types, "type", "TCGA_type", manifest)
    # 2. molecular subtype: Xena Subtype_Selected
    sub = pd.read_csv(SUBTYPE, sep="\t")
    bc2sub = dict(zip(sub["sampleID"].astype(str).str[:15], sub["Subtype_Selected"]))
    subt = pd.Series({c: (bc2sub.get(c[:15]) if pd.notna(bc2sub.get(c[:15])) else None) for c in expr.columns})
    score_and_write(expr, subt, "subtype", "TCGASubtype", manifest)
    # 3. sample type: from the barcode (primary / metastatic / normal / …)
    st = pd.Series({c: SAMPLE_TYPE.get(c[13:15]) for c in expr.columns})
    score_and_write(expr, st, "sampletype", "SampleType", manifest)

    (OUT / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))
    # tumor-type code -> full name, for the dashboard's expanded label (like DepMap modeltype_desc)
    types_present = manifest["divisions"]["type"]["groups"].keys()
    (OUT / "type_desc.json").write_text(json.dumps(
        {c: TCGA_TYPE_NAMES[c] for c in types_present if c in TCGA_TYPE_NAMES}, separators=(",", ":")))
    # cross-group breadth (TCGA-specific): # of type+subtype groups where the TF is a specific MTF
    breadth = {}
    for div in ("type", "subtype"):
        t = pd.read_csv(OUT / f"mtfs_{div}.tsv", sep="\t")
        for tf in t.loc[t.category == "specific", "tf"]:
            breadth[tf] = breadth.get(tf, 0) + 1
    (OUT / "breadth.json").write_text(json.dumps(breadth, separators=(",", ":")))
    print(f"wrote tcga/: {list(manifest['divisions'])}; breadth for {len(breadth)} TFs")


if __name__ == "__main__":
    main()
