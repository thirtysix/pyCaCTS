# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Tests for the empirical-null FDR (pycacts.stats) and the FDR + abundance-floor MTF call
(pycacts.filter.mtf_categories_fdr). Run with `pytest`."""
import numpy as np
import pandas as pd
from pycacts import score, stats, filter as cfilter


def _panel():
    """A deterministic 5-group panel: 20 group-specific TFs (S*, peak ~8 in one group, ~0.3 elsewhere) and
    20 ubiquitous TFs (U*, flat 2.0-3.9). All values are non-negative log2(TPM+1)-style expression."""
    g = [f"g{j}" for j in range(5)]
    rows = {}
    for i in range(20):
        v = [0.3] * 5; v[i % 5] = 8.0; rows[f"S{i}"] = v          # specific to group i % 5
    for i in range(20):
        rows[f"U{i}"] = [2.0 + i * 0.1] * 5                       # flat, not group-specific
    return pd.DataFrame(rows).T.set_axis(g, axis=1)


def test_empirical_fdr_flags_specific():
    """A TF that peaks in its group is significant there (FDR < 0.10); a flat TF is not."""
    S = score.cacts_score_matrix(_panel())
    fdr = stats.empirical_fdr(S["g0"])
    assert fdr["S0"] < 0.10                                        # S0 peaks in g0
    assert fdr["U0"] > 0.10                                        # ubiquitous, not group-specific


def test_fdr_monotonic_in_score():
    """FDR is a left-tail p of the CaCTS score, so it is non-decreasing as the score increases."""
    S = score.cacts_score_matrix(_panel())
    lf = stats.empirical_fdr_log10(S["g0"])
    lf_by_score = lf.reindex(S["g0"].sort_values().index).to_numpy()
    assert np.all(np.diff(lf_by_score) >= -1e-9)                  # more specific (lower score) => lower FDR


def test_empirical_fdr_log10_consistency():
    """empirical_fdr == 10**empirical_fdr_log10, and log10(FDR) <= 0."""
    S = score.cacts_score_matrix(_panel())
    lf = stats.empirical_fdr_log10(S["g0"])
    assert np.allclose(stats.empirical_fdr(S["g0"]).to_numpy(), np.power(10.0, lf.to_numpy()))
    assert (lf <= 0).all()


def test_mtf_categories_fdr_gates():
    """specific = FDR <= max AND mean >= floor; a significant-but-silent TF is excluded, and a highly
    expressed non-specific TF lands in non_specific (the ubiquitous category)."""
    R = _panel()
    R.loc["SILENT"] = [0.3, 0.0, 0.0, 0.0, 0.0]                   # peaks in g0 but below 1 TPM
    R.loc["BROAD"] = [12.0] * 5                                   # highest expression, not group-specific
    S = score.cacts_score_matrix(R)
    cats = cfilter.mtf_categories_fdr(S, R, "g0", fdr_max=0.10, expr_floor=1.0)
    assert "S0" in cats["specific"]                               # significant AND expressed
    assert "SILENT" not in cats["specific"]                      # significant but < 1 TPM -> excluded by floor
    assert "BROAD" in cats["non_specific"]                        # high expression, not specific -> ubiquitous
    assert "BROAD" not in cats["specific"]
