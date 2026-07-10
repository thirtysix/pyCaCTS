# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Loaders for pyCaCTS inputs (DepMap/CCLE expression, cell-line annotation, TF universes)."""
from __future__ import annotations
import re
from pathlib import Path
import pandas as pd

_ENTREZ = re.compile(r"\s*\(\d+\)$")


def load_expression(path: str | Path) -> pd.DataFrame:
    """DepMap OmicsExpression (cell lines x genes, headers 'SYMBOL (ENTREZ)').
    Returns genes x cell-lines (log2 TPM+1), gene symbols de-suffixed & upper-cased, dup symbols meaned."""
    df = pd.read_csv(path, index_col=0)                     # lines x genes
    df.columns = [_ENTREZ.sub("", c).strip().upper() for c in df.columns]
    df = df.T                                               # genes x lines
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def load_model(path: str | Path) -> pd.DataFrame:
    """DepMap Model.csv indexed by ModelID (all columns kept, so any column can be a grouping)."""
    return pd.read_csv(path).set_index("ModelID")


def load_tf_universe(path: str | Path, kind: str = "cacts") -> list[str]:
    """kind='cacts' -> merged 1,671 TF list (one symbol per line, header 'NameTF').
       kind='lambert' -> Lambert 2018 DatabaseExtract, 'Is TF? == Yes' (1,639)."""
    if kind == "cacts":
        tfs = [l.strip().upper() for l in open(path)][1:]
        return sorted({t for t in tfs if t})
    if kind == "lambert":
        d = pd.read_csv(path)
        col = [c for c in d.columns if c.strip().lower() == "is tf?"][0]
        sym = [c for c in d.columns if "symbol" in c.lower() or c.strip().lower() == "hgnc symbol"][0]
        keep = d[d[col].astype(str).str.strip().str.lower() == "yes"]
        return sorted({str(x).strip().upper() for x in keep[sym].dropna()})
    raise ValueError(f"unknown TF universe kind: {kind}")


def load_essentiality(path: str | Path) -> pd.DataFrame:
    """DepMap CRISPRGeneEffect.csv -> cell lines x genes (Chronos), gene symbols de-suffixed."""
    df = pd.read_csv(path, index_col=0)
    df.columns = [_ENTREZ.sub("", c).strip().upper() for c in df.columns]
    return df
