#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""pyCaCTS on TCGA: the original CaCTS analysis (master TFs per tumor type) from downloaded files.

Download the two inputs first (see the README 'Worked example: master TFs of the TCGA tumor types'):

  1. TCGA pan-cancer expression, log2, from UCSC Xena PanCanAtlas (~1.6 GB):
     curl -L -o TCGA_pancan.geneExp.gz \\
       "https://pancanatlas.xenahubs.net/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"

  2. Sample -> tumor-type map (the 34-type list CaCTS used), from the original authors' repo:
     curl -L -o TCGA_sample_types.txt \\
       "https://raw.githubusercontent.com/lawrenson-lab/CaCTS/master/files/SuppTable1-34-TCGAID.txt"

Then:  python examples/tcga.py TCGA_pancan.geneExp.gz TCGA_sample_types.txt [TUMOR_TYPE]   # default SKCM
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pycacts import io, score, filter as cfilter

expr_path = sys.argv[1] if len(sys.argv) > 1 else "TCGA_pancan.geneExp.gz"
map_path = sys.argv[2] if len(sys.argv) > 2 else "TCGA_sample_types.txt"
TYPE = sys.argv[3] if len(sys.argv) > 3 else "SKCM"          # cutaneous melanoma; any of the 34 codes
TF = ROOT / "data" / "CaCTS_merged_1671_TFs.txt"

print("loading TCGA expression ...")
expr = pd.read_csv(expr_path, sep="\t", index_col=0)         # genes x samples, log2
expr.index = expr.index.astype(str).str.upper()

# map each sample to its tumor type; match on the 15-char sample barcode (TCGA-XX-XXXX-01)
smap = pd.read_csv(map_path, sep="\t")                       # columns: Cancer, Category, SampleId
sample2type = dict(zip(smap["SampleId"].str[:15], smap["Cancer"]))
cols = [c for c in expr.columns if c[:15] in sample2type]
types = pd.Series({c: sample2type[c[:15]] for c in cols})
print(f"matched {len(cols)} samples across {types.nunique()} tumor types")

# subset to the CaCTS TF universe, per-type mean = the CaCTS representative matrix, then score
tfs = io.load_tf_universe(TF, "cacts")
rep = expr.loc[expr.index.intersection(tfs), cols].T.groupby(types).mean().T
scores = score.cacts_score_matrix(rep)                      # TFs x tumor types

print(f"\n{TYPE}: most group-specific TFs (by CaCTS score):")
print(score.rank_specific(scores, TYPE).head(10).to_string(index=False))
mtfs = cfilter.mtf_categories(scores, rep, TYPE)
print(f"\nspecific MTFs: {', '.join(mtfs['specific'])}")
