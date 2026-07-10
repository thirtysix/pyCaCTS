# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""CaCTS MTF categories (port of `Filter.R`).

For a query group: `specific` = top-pn by CaCTS score (lowest JSD) INTERSECT top-pnE by expression in the
group; `non_specific` = top-pnE by expression MINUS the specific set (high expression, low specificity =
their Table S7 category). N = size of the TF universe (rep_matrix rows).
"""
from __future__ import annotations
import pandas as pd


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
