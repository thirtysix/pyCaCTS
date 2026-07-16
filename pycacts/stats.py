# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Empirical-null significance for CaCTS scores.

For one group, model a null from the *non-specific* (high-JSD) side of that group's own score distribution
(robust median + one-sided MAD), take a left-tail p-value per TF, and Benjamini-Hochberg correct to an FDR.
It answers "is this TF significantly more group-specific than the background?", a data-driven alternative to
the fixed top-5%-by-score cutoff. Computed in natural-log space so the most specific TFs stay distinct
(a normal CDF underflows to 0 for z below about -8); returns log10(FDR) so callers threshold at, e.g.,
log10(0.10) = -1. Numerically matches the dashboard's util.js `empiricalFDR`.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd

_erf = np.vectorize(math.erf)
_LN2PI = math.log(2 * math.pi)


def _ln_left_tail(z: np.ndarray) -> np.ndarray:
    """ln of the standard-normal left tail Phi(z), stable for very negative z."""
    tiny = np.finfo(float).tiny
    central = np.log(np.clip(0.5 * (1.0 + _erf(z / math.sqrt(2))), tiny, None))
    asymp = -0.5 * z * z - np.log(-np.minimum(z, -tiny)) - 0.5 * _LN2PI   # asymptotic ln Phi(z), z << 0
    return np.where(z > -6, central, asymp)


def empirical_fdr_log10(scores: pd.Series) -> pd.Series:
    """log10(empirical-null FDR) per TF for one group's CaCTS scores (lower score = more specific).

    Returns a Series (index = TF) of values <= 0; log10(FDR) <= -1 means FDR <= 0.10.
    """
    x = np.asarray(scores, dtype=float)
    mu = float(np.median(x))
    dev = np.abs(x[x >= mu] - mu)
    sig = 1.4826 * float(np.median(dev)) if dev.size else 0.0
    if not sig:
        sig = 1e-9
    lp = _ln_left_tail((x - mu) / sig)                    # ln(left-tail p) per TF
    n = lp.size
    order = np.argsort(lp)                                # most significant (smallest ln p) first
    out = np.zeros(n)
    prev = 0.0                                            # ln of running BH minimum (ln 1 = 0)
    for k in range(n, 0, -1):                             # BH step-up: ln(fdr) = ln p + ln n - ln rank
        i = order[k - 1]
        prev = min(prev, lp[i] + math.log(n) - math.log(k))
        out[i] = min(0.0, prev)
    return pd.Series(out / math.log(10), index=scores.index)   # natural-log -> log10


def empirical_fdr(scores: pd.Series) -> pd.Series:
    """Empirical-null FDR per TF (linear, in (0, 1]). Convenience wrapper over `empirical_fdr_log10`."""
    return np.power(10.0, empirical_fdr_log10(scores))
