# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Build the CaCTS 'representative-sample' matrix (per-group MEAN expression) at three resolutions.

CaCTS `prepare.R` collapses each group to one representative column = the mean expression of its samples;
that per-group mean is also how it balances imbalanced group sizes. We support grouping by individual
cell line, OncotreeSubtype, or OncotreeLineage from the same expression matrix.
"""
from __future__ import annotations
import pandas as pd

_GROUPCOL = {"lineage": "OncotreeLineage", "subtype": "OncotreeSubtype"}   # convenience aliases


def build_rep_matrix(expr: pd.DataFrame, model: pd.DataFrame, grouping: str,
                     tf_universe: list[str] | None = None, min_group_n: int = 1) -> tuple[pd.DataFrame, pd.Series]:
    """expr: genes x cell-lines. model: indexed by ModelID.
    grouping: 'line' (each cell line its own group), an alias in _GROUPCOL, or ANY Model.csv column name.
    min_group_n: drop groups with fewer than this many lines (1 = keep all).
    Returns (rep_matrix TFs-or-genes x groups [per-group mean], group_size Series).
    """
    lines = [c for c in expr.columns if c in model.index]      # keep annotated lines
    X = expr[lines]
    if grouping == "line":
        rep = X.copy()
        gsize = pd.Series(1, index=rep.columns)
    else:
        col = _GROUPCOL.get(grouping, grouping)
        if col not in model.columns:
            raise ValueError(f"grouping column {col!r} not in Model.csv")
        grp = model.loc[lines, col]
        keep = grp.dropna()
        X = X[keep.index]
        rep = X.T.groupby(keep).mean().T                        # genes x groups (per-group mean)
        gsize = keep.value_counts()
        big = gsize[gsize >= min_group_n].index
        rep = rep[[g for g in rep.columns if g in set(big)]]
        gsize = gsize[rep.columns]
    if tf_universe is not None:
        rep = rep.loc[[g for g in rep.index if g in set(tf_universe)]]
    return rep, gsize
