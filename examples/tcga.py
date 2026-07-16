#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""pyCaCTS on TCGA: master TFs per tumor type, from downloaded files, using the same method as the dashboard
(empirical-null FDR < 0.10 AND mean >= 1 TPM).

Download three inputs first (see the README 'Worked example: master TFs of the TCGA tumor types'):

  1. TCGA RSEM gene-TPM matrix (log2(TPM+0.001), Ensembl-keyed), from the UCSC Xena Toil hub (~740 MB):
     curl -L -o tcga_RSEM_gene_tpm.gz \\
       "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_tpm.gz"
  2. gencode.v23 gene probemap (Ensembl id -> gene symbol):
     curl -L -o gencode.v23.gene.probemap \\
       "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap/gencode.v23.annotation.gene.probemap"
  3. Sample -> tumor-type map (the 33-code tumor-type list CaCTS used), from the original authors' repo:
     curl -L -o TCGA_sample_types.txt \\
       "https://raw.githubusercontent.com/lawrenson-lab/CaCTS/master/files/SuppTable1-34-TCGAID.txt"

Then:  python examples/tcga.py tcga_RSEM_gene_tpm.gz gencode.v23.gene.probemap TCGA_sample_types.txt [TYPE]
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pycacts import io, score, filter as cfilter, stats

tpm_path = sys.argv[1] if len(sys.argv) > 1 else "tcga_RSEM_gene_tpm.gz"
map_path = sys.argv[2] if len(sys.argv) > 2 else "gencode.v23.gene.probemap"
type_path = sys.argv[3] if len(sys.argv) > 3 else "TCGA_sample_types.txt"
TYPE = sys.argv[4] if len(sys.argv) > 4 else "SKCM"           # cutaneous melanoma; any of the 33 codes
TF = ROOT / "data" / "CaCTS_merged_1671_TFs.txt"

tfs = set(io.load_tf_universe(TF, "cacts"))
ens2sym = dict(pd.read_csv(map_path, sep="\t", usecols=["id", "gene"]).itertuples(index=False))
keep = {e for e, s in ens2sym.items() if str(s).upper() in tfs}

print("loading TCGA Toil TPM (TF rows only) ...")
chunks = [ch[ch.index.isin(keep)] for ch in
          pd.read_csv(tpm_path, sep="\t", index_col=0, chunksize=4000)]
log2p001 = pd.concat(chunks)                                   # log2(TPM+0.001), Ensembl x samples
tpm = np.clip(np.power(2.0, log2p001.to_numpy(float)) - 0.001, 0.0, None)
expr = pd.DataFrame(np.log2(tpm + 1.0), index=[str(ens2sym[e]).upper() for e in log2p001.index],
                    columns=log2p001.columns).groupby(level=0).mean()   # log2(TPM+1), same units as DepMap

smap = pd.read_csv(type_path, sep="\t")                        # columns: Cancer, Category, SampleId
sample2type = dict(zip(smap["SampleId"].str[:15], smap["Cancer"]))
cols = [c for c in expr.columns if c[:15] in sample2type]
types = pd.Series({c: sample2type[c[:15]] for c in cols})
print(f"matched {len(cols)} samples across {types.nunique()} tumor types")

rep = expr[cols].T.groupby(types).mean().T                    # TFs x tumor types (mean log2(TPM+1))
scores = score.cacts_score_matrix(rep)

print(f"\n{TYPE}: most group-specific TFs (by CaCTS score):")
print(score.rank_specific(scores, TYPE).head(10).to_string(index=False))
# specific MTF = empirical-null FDR < 0.10 AND mean >= 1 TPM (log2(TPM+1) >= 1.0), same as the dashboard
mtfs = cfilter.mtf_categories_fdr(scores, rep, TYPE, fdr_max=0.10, expr_floor=1.0)
print(f"\nspecific MTFs (FDR<0.10 & >=1 TPM): {', '.join(mtfs['specific'])}")
