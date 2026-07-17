#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Build a TCGA patient -> AJCC major-stage map for the within-cancer Stage axis.

Queries the GDC clinical API for every TCGA case's `diagnoses.ajcc_pathologic_stage` and collapses it to a
major stage (I / II / III / IV); cases with no AJCC stage (e.g. GBM, LGG, LAML, PRAD, which stage/grade
differently) are dropped. Writes patient<TAB>stage to PYCACTS_TCGA_STAGE (default data/tcga/TCGA_stage.tsv).
A pre-downloaded GDC cases TSV can be supplied via PYCACTS_GDC_TSV to avoid the live query.
"""
import os, re, io, json, sys
import urllib.request, urllib.parse
from pathlib import Path
import pandas as pd

OUT = Path(os.environ.get("PYCACTS_TCGA_STAGE", "data/tcga/TCGA_stage.tsv"))
GDC_TSV = os.environ.get("PYCACTS_GDC_TSV")                     # optional pre-downloaded GDC cases TSV
API = "https://api.gdc.cancer.gov/cases"


def major_stage(v):
    """'Stage IIIA' / 'stage iv' -> 'III' / 'IV'; anything without a I-IV roman numeral -> None."""
    if not isinstance(v, str):
        return None
    m = re.match(r"stage\s+(iv|iii|ii|i)", v.strip().lower())   # longest-first so IVa->IV, IIIa->III, Ib->I
    return m.group(1).upper() if m else None


def main():
    if GDC_TSV:
        txt = Path(GDC_TSV).read_text()
    else:
        params = {"filters": json.dumps({"op": "in",
                    "content": {"field": "project.program.name", "value": ["TCGA"]}}),
                  "fields": "submitter_id,diagnoses.ajcc_pathologic_stage", "size": "15000", "format": "tsv"}
        print("querying GDC clinical API ...")
        txt = urllib.request.urlopen(API + "?" + urllib.parse.urlencode(params), timeout=180).read().decode()
    d = pd.read_csv(io.StringIO(txt), sep="\t", dtype=str)
    scol = [c for c in d.columns if c.endswith("ajcc_pathologic_stage")]   # diagnoses.0/1/2.ajcc_pathologic_stage
    d["stage"] = d[scol].apply(lambda r: next((s for s in (major_stage(v) for v in r) if s), None), axis=1)
    out = (d.dropna(subset=["stage"])[["submitter_id", "stage"]]
           .rename(columns={"submitter_id": "patient"}).drop_duplicates("patient"))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {len(out)} patients with a major stage -> {OUT}")
    print("  stage counts:", out["stage"].value_counts().to_dict())


if __name__ == "__main__":
    main()
