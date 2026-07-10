#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Decompose where pyCaCTS's speed-up over the original CaCTS R comes from.

Times three Python implementations of the identical score on the same matrices, a "ladder":
  1. naive       : the literal R algorithm: per-query, per-TF Python loops (O(T·G²), interpreted).
  2. per-query    : vectorize the naive loop over TFs, still one pass per query group (O(T·G²), compiled).
  3. closed-form  : pyCaCTS: the whole TF×group matrix at once via the exact decomposition (O(T·G)).

naive→per-query isolates the *vectorization* gain (constant factor); per-query→closed-form isolates the
*algorithmic* gain (a factor that grows with the number of groups). The original R sits above naive: it runs
the same O(T·G²) loop but adds per-query file I/O and rbind-in-loop growth (pass R times from
benchmark_vs_r.py for the full picture). Reads DepMap from PYCACTS_DEPMAP.
"""
import os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, grouping, score
from pycacts.score import _normalize_rows, _score_matrix_naive, cacts_score_matrix, EPS

DATA = HERE.parent / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))
EXPR = DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL = DEPMAP / "Model.csv"
TF = DATA / "CaCTS_merged_1671_TFs.txt"
LEVELS = [("lineage", "lineage"), ("subtype", "subtype"), ("line", "line")]


def score_perquery(rep):
    """Vectorized over TFs, one pass per query group (O(T·G²), compiled): the middle rung."""
    obs = _normalize_rows(rep.values)                            # T x G
    T, G = obs.shape
    out = np.empty((T, G))
    with np.errstate(divide="ignore", invalid="ignore"):         # log(0) is masked by the where below
        for q in range(G):
            ideal = np.full(G, EPS); ideal[q] = 1.0              # 1 x G
            m = 0.5 * (obs + ideal)                              # T x G (broadcast)
            ot = np.where(obs > 0, obs * np.log(obs / m), 0.0).sum(axis=1)
            it = (ideal * np.log(ideal / m)).sum(axis=1)         # ideal = EPS > 0 everywhere
            out[:, q] = 0.5 * ot + 0.5 * it
    return pd.DataFrame(out, index=rep.index, columns=rep.columns)


def timeit(fn, rep, n):
    ts = []
    for _ in range(n):
        t = time.perf_counter(); fn(rep); ts.append(time.perf_counter() - t)
    return float(np.median(ts))


def main():
    print("loading DepMap + building representative matrices ...")
    expr = io.load_expression(EXPR); model = io.load_model(MODEL); tfs = io.load_tf_universe(TF, "cacts")
    rows = []
    for key, gcol in LEVELS:
        rep, _ = grouping.build_rep_matrix(expr, model, gcol, tf_universe=tfs, min_group_n=1)
        rep = rep.dropna(how="all")
        T, G = rep.shape
        # correctness: all three agree
        fast = cacts_score_matrix(rep).values
        pq = score_perquery(rep).values
        assert np.nanmax(np.abs(fast - pq)) < 1e-9, "per-query != closed-form"
        # naive is O(T·G²) interpreted; only time it where it finishes quickly
        t_closed = timeit(cacts_score_matrix, rep, 5)
        t_pq = timeit(score_perquery, rep, 3)
        t_naive = timeit(_score_matrix_naive, rep, 1) if G <= 200 else None
        row = dict(level=key, T=T, G=G,
                   naive_s=t_naive, perquery_s=t_pq, closed_s=t_closed,
                   vec_gain=(t_naive / t_pq) if t_naive else None,          # naive -> per-query
                   algo_gain=(t_pq / t_closed))                             # per-query -> closed-form
        rows.append(row)
        n = f"{t_naive*1000:9.1f}" if t_naive else "     n/a"
        print(f"[{key:8s}] T={T} G={G:4d} | naive {n} ms | per-query {t_pq*1000:8.2f} ms | "
              f"closed {t_closed*1000:7.3f} ms | vec×{(t_naive/t_pq):.0f} algo×{t_pq/t_closed:.0f}" if t_naive
              else f"[{key:8s}] T={T} G={G:4d} | naive      n/a | per-query {t_pq*1000:8.2f} ms | "
                   f"closed {t_closed*1000:7.3f} ms | algo×{t_pq/t_closed:.0f}")
    tab = pd.DataFrame(rows)
    (HERE.parent / "results" / "profile_speedup.csv").write_text(tab.to_csv(index=False))
    print("\n=== ladder (each step is an independent multiplicative gain) ===")
    print(tab.to_string(index=False))
    print("\nnaive→per-query = vectorization (≈constant); per-query→closed-form = algorithm (grows with G).")


if __name__ == "__main__":
    main()
