# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Validate the fast closed-form CaCTS score against the literal per-query transcription, and check that
lower scores mean more group-specific. Run with `pytest`."""
import numpy as np
import pandas as pd
from pycacts.score import cacts_score_matrix, _score_matrix_naive


def _toy(seed=0, T=40, G=8):
    rng = np.random.default_rng(seed)
    m = rng.random((T, G)) * rng.integers(0, 3, (T, G))          # sparse-ish non-negative expression
    return pd.DataFrame(m, index=[f"TF{i}" for i in range(T)], columns=[f"G{j}" for j in range(G)])


def test_fast_equals_naive():
    """The O(T·G) closed form must equal the literal O(T·G²) R transcription to float precision."""
    for seed in range(5):
        R = _toy(seed)
        fast = cacts_score_matrix(R).values
        naive = _score_matrix_naive(R).values
        assert np.nanmax(np.abs(fast - naive)) < 1e-9


def test_lower_score_is_more_specific():
    """A group-restricted TF scores lower (more specific) than a ubiquitous one, and lowest in its group."""
    R = pd.DataFrame(0.0, index=["specific", "ubiquitous"], columns=["g0", "g1", "g2"])
    R.loc["specific"] = [10.0, 0.0, 0.0]                          # only expressed in g0
    R.loc["ubiquitous"] = [3.0, 3.0, 3.0]                         # expressed everywhere
    S = cacts_score_matrix(R)
    assert S.loc["specific", "g0"] < S.loc["ubiquitous", "g0"]   # more specific to g0
    assert S.loc["specific", "g0"] < S.loc["specific", "g1"]     # most specific to its own group
