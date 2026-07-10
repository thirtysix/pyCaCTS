#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Validate pyCaCTS against the original CaCTS R and benchmark speed on identical input.

Builds real representative matrices from CCLE at two scales (lineage ~ TCGA scale: ~1.6k TFs x ~29
groups; and per-line: ~1.6k x ~1.45k), scores each with (a) pyCaCTS (vectorized Python) and (b) the
original run_CaCTS_score R (looped), and reports max score difference + wall-clock speedup.
Run with the shared analysis venv (R must be on PATH; reference_R/JSD.R present locally).
"""
import os, sys, time, subprocess, tempfile
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from pycacts import io, grouping, score

DATA = HERE.parent / "data"
DEPMAP = Path(os.environ.get("PYCACTS_DEPMAP", DATA / "depmap"))    # DepMap release dir (download; see README)
EXPR = DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
MODEL = DEPMAP / "Model.csv"; TF_CACTS = DATA / "CaCTS_merged_1671_TFs.txt"
JSD_R = HERE.parent / "reference_R" / "JSD.R"                       # original CaCTS R (GPL; not redistributed, see README)
BENCH = HERE.parent / "results" / "bench"; BENCH.mkdir(parents=True, exist_ok=True)


def safe_cols(df):  # their R writes 'query-CaCTS-scores.txt', so no '/' or spaces in group names
    return df.rename(columns={c: c.replace("/", "_").replace(" ", "_") for c in df.columns})


def run_r(rep, tag):
    rp, op = BENCH / f"rep_{tag}.tsv", BENCH / f"rscore_{tag}.tsv"
    safe_cols(rep).to_csv(rp, sep="\t")
    r = subprocess.run(["Rscript", str(HERE / "cacts_reference.R"), str(rp), str(op), str(JSD_R)],
                       cwd=BENCH, capture_output=True, text=True)
    if r.returncode != 0:
        print("  R FAILED:\n" + r.stdout + r.stderr); return None, None
    rt = float(next(l.split()[1] for l in r.stdout.splitlines() if l.startswith("R_ELAPSED_SECONDS")))
    rs = pd.read_csv(op, sep="\t", index_col=0)
    return rs, rt


def bench(rep, tag):
    rep = rep.dropna(how="all")
    print(f"\n=== {tag}: {rep.shape[0]} TFs x {rep.shape[1]} groups ===")
    # pyCaCTS (median of 3 timed runs)
    py_t = []
    for _ in range(3):
        t = time.perf_counter(); Spy = score.cacts_score_matrix(rep); py_t.append(time.perf_counter() - t)
    pyt = float(np.median(py_t))
    Spy = safe_cols(Spy)
    # original R
    Sr, rt = run_r(rep, tag)
    if Sr is None: return
    common_c = [c for c in Spy.columns if c in Sr.columns]
    A = Spy.loc[Sr.index, common_c].values; B = Sr.loc[Sr.index, common_c].values
    maxdiff = float(np.nanmax(np.abs(A - B)))
    print(f"  pyCaCTS: {pyt*1000:8.2f} ms   |   original R: {rt*1000:9.1f} ms   |   speedup ~{rt/pyt:,.0f}x")
    print(f"  max |pyCaCTS - R| over {A.shape[0]}x{A.shape[1]} scores = {maxdiff:.2e}  "
          f"-> {'IDENTICAL' if maxdiff < 1e-9 else 'MISMATCH'}")
    return dict(tag=tag, shape=rep.shape, py_ms=pyt*1000, r_ms=rt*1000, speedup=rt/pyt, maxdiff=maxdiff)


def main():
    print("loading CCLE expression + building representative matrices ...")
    expr = io.load_expression(EXPR); model = io.load_model(MODEL); tfs = io.load_tf_universe(TF_CACTS, "cacts")
    rep_lin, _ = grouping.build_rep_matrix(expr, model, "lineage", tf_universe=tfs)   # ~TCGA scale
    rep_line, _ = grouping.build_rep_matrix(expr, model, "line", tf_universe=tfs)     # per-line (pyCaCTS extension)
    rows = []
    rows.append(bench(rep_lin, "lineage_TCGA-scale"))
    rows.append(bench(rep_line.iloc[:, :400], "per-line_400"))   # subset: full 1450 in R is very slow
    tab = pd.DataFrame([r for r in rows if r]).to_string(index=False)
    print("\n=== summary ===\n" + tab)
    (HERE.parent / "results" / "benchmark_summary.txt").write_text(tab + "\n")


if __name__ == "__main__":
    main()
