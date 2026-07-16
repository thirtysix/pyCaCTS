#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Convert the UCSC Xena Toil TCGA RSEM gene-TPM matrix into a gene-symbol-keyed log2(TPM+1) matrix,
so TCGA expression is in the SAME units as DepMap's OmicsExpressionProteinCodingGenesTPMLogp1 and a literal
1 TPM abundance floor applies identically to both.

Toil `tcga_RSEM_gene_tpm` is Ensembl-gene-keyed and stored as log2(TPM+0.001); this maps Ensembl->symbol
via the gencode.v23 probemap, keeps only the CaCTS TF universe, un-logs to linear TPM, and re-encodes as
log2(TPM+1). Output: a TFs x samples TSV (gzip) usable as PYCACTS_TCGA_EXPR for stage_tcga(_within).py.

Env / paths:
  PYCACTS_TOIL_TPM    tcga_RSEM_gene_tpm.gz            (Xena Toil hub)
  PYCACTS_TOIL_MAP    gencode.v23 gene probemap        (id -> gene symbol)
  PYCACTS_TCGA_OUT    output path (default data/tcga/TCGA_toil_tpm_log2p1.tsv.gz)
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io

DATA = HERE.parent / "data"
TPM = Path(os.environ["PYCACTS_TOIL_TPM"])
PMAP = Path(os.environ["PYCACTS_TOIL_MAP"])
OUT = Path(os.environ.get("PYCACTS_TCGA_OUT", DATA / "tcga" / "TCGA_toil_tpm_log2p1.tsv.gz"))
TF = DATA / "CaCTS_merged_1671_TFs.txt"
OFFSET = 0.001                                                   # Toil encoding: log2(TPM + 0.001)


def main():
    tfs = set(io.load_tf_universe(TF, "cacts"))
    pm = pd.read_csv(PMAP, sep="\t", usecols=["id", "gene"])
    ens2sym = dict(zip(pm["id"], pm["gene"]))
    keep_ens = {e for e, s in ens2sym.items() if str(s).upper() in tfs}   # Ensembl IDs that are CaCTS TFs
    print(f"CaCTS TFs {len(tfs)}; probemap Ensembl IDs mapping to a TF: {len(keep_ens)}")

    chunks = []
    reader = pd.read_csv(TPM, sep="\t", index_col=0, chunksize=4000, dtype={0: str})
    seen = 0
    for ch in reader:
        seen += len(ch)
        hit = ch[ch.index.isin(keep_ens)]
        if len(hit):
            chunks.append(hit.astype("float32"))
    log2p001 = pd.concat(chunks)
    print(f"scanned {seen} genes; kept {log2p001.shape[0]} TF rows x {log2p001.shape[1]} samples")

    # log2(TPM+0.001) -> linear TPM -> log2(TPM+1), matching DepMap units
    tpm = np.power(2.0, log2p001.to_numpy(dtype="float64")) - OFFSET
    np.clip(tpm, 0.0, None, out=tpm)
    log2p1 = pd.DataFrame(np.log2(tpm + 1.0), index=log2p001.index, columns=log2p001.columns)

    log2p1.index = [str(ens2sym[e]).upper() for e in log2p1.index]        # Ensembl -> symbol
    log2p1 = log2p1.groupby(level=0).mean()                               # collapse duplicate symbols
    log2p1.index.name = "gene"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    log2p1.round(4).to_csv(OUT, sep="\t")
    print(f"wrote {log2p1.shape[0]} TF symbols x {log2p1.shape[1]} samples -> {OUT}")
    print("sanity (should be >=0, ~log2(TPM+1)):", f"min={log2p1.values.min():.3f} max={log2p1.values.max():.3f}")


if __name__ == "__main__":
    main()
