# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""CaCTS specificity score (Jensen-Shannon divergence), Python port of lawrenson-lab/CaCTS `JSD.R`.

For a TF, let `obs` be its representative-sample (per-group mean) expression profile normalized to sum 1
across the G groups, and `ideal_q` be the one-hot distribution for query group q (1 at q, EPS elsewhere,
EPS=1e-17). The CaCTS score of the TF for group q is the Jensen-Shannon divergence

    JSD(obs, ideal_q) = 0.5*KL(obs || m) + 0.5*KL(ideal_q || m),   m = 0.5*(obs + ideal_q)

with natural logs. **Lower JSD = more group-specific** (rank 1 = smallest JSD), matching the R.

`cacts_score_matrix` returns the full TF x group JSD matrix in O(T*G) via an exact decomposition
(so per-cell-line grouping with ~1,450 groups is tractable); `_score_matrix_naive` is the literal,
per-query R transcription used only to validate the fast path.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

EPS = 1e-17  # `rep.NA` / ideal background value used by CaCTS


def _normalize_rows(mat: np.ndarray) -> np.ndarray:
    """Row-normalize a (TF x group) matrix to sum 1 (CaCTS `obs`). NaNs treated as 0."""
    m = np.nan_to_num(mat, nan=0.0).astype(float)
    s = m.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return m / s


def _xlogxy(x: np.ndarray, denom: np.ndarray) -> np.ndarray:
    """x * log(x / denom) with the convention 0*log(0/·)=0 (elementwise)."""
    out = np.zeros_like(x)
    nz = x > 0
    out[nz] = x[nz] * np.log(x[nz] / denom[nz])
    return out


def cacts_score_matrix(rep_matrix: pd.DataFrame) -> pd.DataFrame:
    """CaCTS JSD score for every TF (rows) x group (cols). Lower = more specific.

    rep_matrix: representative-sample matrix (per-group MEAN expression), TFs x groups.
    Returns a DataFrame of JSD scores with the same index/columns.
    """
    obs = _normalize_rows(rep_matrix.values)                      # T x G
    # obs term pieces:  A = obs*log(obs/(0.5*obs+0.5*EPS))  (non-query form, all j)
    A = _xlogxy(obs, 0.5 * obs + 0.5 * EPS)
    B = _xlogxy(obs, 0.5 * obs + 0.5)                             # query form (ideal=1 at q)
    SA = A.sum(axis=1, keepdims=True)
    # ideal term pieces: logterm = 1*log(1/(0.5*obs+0.5))  at the query column
    logterm = np.log(1.0 / (0.5 * obs + 0.5))
    c = EPS * np.log(EPS / (0.5 * obs + 0.5 * EPS))               # background j-term (subtracted at q)
    C = c.sum(axis=1, keepdims=True)
    jsd = 0.5 * (SA - A + B) + 0.5 * (logterm + C - c)            # T x G, exact
    return pd.DataFrame(jsd, index=rep_matrix.index, columns=rep_matrix.columns)


def _score_matrix_naive(rep_matrix: pd.DataFrame) -> pd.DataFrame:
    """Literal transcription of the R loop (slow); for validating cacts_score_matrix only."""
    obs_all = _normalize_rows(rep_matrix.values)
    T, G = obs_all.shape
    out = np.zeros((T, G))
    with np.errstate(divide="ignore", invalid="ignore"):         # log(0) is masked by the where below
        for i in range(T):
            obs = obs_all[i]
            for q in range(G):
                ideal = np.full(G, EPS); ideal[q] = 1.0
                m = 0.5 * (obs + ideal)
                ot = np.where(obs > 0, obs * np.log(obs / m), 0.0).sum()
                it = np.where(ideal > 0, ideal * np.log(ideal / m), 0.0).sum()
                out[i, q] = 0.5 * ot + 0.5 * it
    return pd.DataFrame(out, index=rep_matrix.index, columns=rep_matrix.columns)


def rank_specific(score_matrix: pd.DataFrame, group: str) -> pd.DataFrame:
    """Per-TF table for one group: CaCTS score (JSD), -log10 score, and rank (1 = most specific)."""
    s = score_matrix[group].sort_values(ascending=True)          # smallest JSD first
    return pd.DataFrame({"tf": s.index, "cacts_score": s.values,
                         "neg_log10": -np.log10(s.values),
                         "rank": np.arange(1, len(s) + 1)}).reset_index(drop=True)


if __name__ == "__main__":  # self-test: fast path must equal the naive R transcription
    rng = np.random.default_rng(0)
    R = pd.DataFrame(rng.random((40, 8)) * rng.integers(0, 3, (40, 8)),
                     index=[f"TF{i}" for i in range(40)], columns=[f"G{j}" for j in range(8)])
    fast, naive = cacts_score_matrix(R).values, _score_matrix_naive(R).values
    md = np.nanmax(np.abs(fast - naive))
    print(f"max |fast - naive| = {md:.2e}  ->  {'OK' if md < 1e-9 else 'MISMATCH'}")
