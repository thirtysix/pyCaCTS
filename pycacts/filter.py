# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""CaCTS MTF categories.

Two ways to call a master TF for a query group:

`mtf_categories` (port of `Filter.R`, the original fixed cutoff): `specific` = top-pn by CaCTS score
(lowest JSD) INTERSECT top-pnE by expression; `non_specific` = top-pnE by expression MINUS the specific set.

`mtf_categories_fdr` (the significance-based call): `specific` = empirical-null FDR <= `fdr_max` INTERSECT a
light abundance floor (mean expression >= `expr_floor`, in the matrix's own units; 1.0 = 1 TPM on a
log2(TPM+1) matrix). This replaces the arbitrary top-5%-by-score cutoff with a data-driven significance
threshold and the aggressive top-5%-expression gate with a floor that only excludes near-silent genes, so
genuinely expressed lineage TFs (e.g. SOX17 / WT1) are kept while low-abundance JSD artifacts are dropped.
`non_specific` stays the CaCTS "candidate ubiquitous master regulator": high expression (top-pnE), not
group-specific.
"""
from __future__ import annotations
import math
import pandas as pd

from .stats import empirical_fdr_log10


def mtf_categories(score_matrix: pd.DataFrame, rep_matrix: pd.DataFrame, group: str,
                   pn: float = 0.05, pnE: float = 0.05) -> dict:
    n = score_matrix.shape[0]
    kJ, kE = round(pn * n), round(pnE * n)
    top_jsd = list(score_matrix[group].nsmallest(kJ).index)      # most specific
    top_expr = list(rep_matrix[group].nlargest(kE).index)        # highest expression in group
    tj, te = set(top_jsd), set(top_expr)
    specific = [t for t in top_jsd if t in te]                   # JSD-ordered
    non_specific = [t for t in top_expr if t not in tj]          # expr-ordered
    return {"group": group, "n_universe": n, "k_jsd": kJ, "k_expr": kE,
            "specific": specific, "non_specific": non_specific,
            "top_jsd": top_jsd, "top_expr": top_expr}


def mtf_categories_fdr(score_matrix: pd.DataFrame, rep_matrix: pd.DataFrame, group: str,
                       fdr_log10: pd.Series | None = None, fdr_max: float = 0.10,
                       expr_floor: float = 1.0, pnE: float = 0.05) -> dict:
    """FDR + light-floor MTF call. `fdr_log10` = precomputed log10(FDR) for this group (else computed here).

    specific     = FDR <= fdr_max AND mean expression >= expr_floor, ordered by CaCTS score (most specific first)
    non_specific = top-pnE by expression, NOT specific (high expression, not group-specific = ubiquitous)
    Returns the category lists plus `fdr_log10` (per TF) so callers can store the significance value.
    """
    if fdr_log10 is None:
        fdr_log10 = empirical_fdr_log10(score_matrix[group])
    n = score_matrix.shape[0]
    kE = round(pnE * n)
    thr = math.log10(fdr_max)
    expr = rep_matrix[group]
    sig = set(fdr_log10.index[fdr_log10 <= thr])
    expressed = set(expr.index[expr >= expr_floor])
    spec = sig & expressed
    specific = list(score_matrix.loc[list(spec), group].sort_values().index)   # by CaCTS score (asc = specific)
    top_expr = list(expr.nlargest(kE).index)
    non_specific = [t for t in top_expr if t not in spec]                       # high expr, not specific
    return {"group": group, "n_universe": n, "fdr_max": fdr_max, "expr_floor": expr_floor,
            "specific": specific, "non_specific": non_specific, "fdr_log10": fdr_log10}
