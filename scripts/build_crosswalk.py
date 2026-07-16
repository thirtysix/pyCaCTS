#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Harlan Barker
"""Build dashboard/data/crosswalk.json: a curated map from each TCGA tumor type to its best-matching
DepMap group, for the 'Tumor vs model' comparison. Each TCGA type points at a DepMap subtype where a clean
histology match exists (adequate sample count), otherwise the DepMap lineage. Validated against both
manifests (a mapping is dropped, with a warning, if either group is missing). Needs stage_dashboard_data.py
and stage_tcga.py to have run."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DD = HERE.parent / "dashboard" / "data"

# TCGA type code -> (DepMap division, DepMap group). Subtype where a clean histology match exists, else lineage.
CROSSWALK = {
    "SKCM": ("subtype", "Cutaneous Melanoma"), "OV": ("subtype", "High-Grade Serous Ovarian Cancer"),
    "GBM": ("subtype", "Glioblastoma"), "LUAD": ("subtype", "Lung Adenocarcinoma"),
    "LUSC": ("subtype", "Lung Squamous Cell Carcinoma"), "COAD": ("subtype", "Colon Adenocarcinoma"),
    "PAAD": ("subtype", "Pancreatic Adenocarcinoma"), "STAD": ("subtype", "Stomach Adenocarcinoma"),
    "BLCA": ("subtype", "Bladder Urothelial Carcinoma"), "LIHC": ("subtype", "Hepatocellular Carcinoma"),
    "KIRC": ("subtype", "Renal Clear Cell Carcinoma"), "LAML": ("subtype", "Acute Myeloid Leukemia"),
    "DLBC": ("subtype", "Diffuse Large B-Cell Lymphoma, NOS"), "UVM": ("subtype", "Uveal Melanoma"),
    "CHOL": ("subtype", "Intrahepatic Cholangiocarcinoma"), "UCEC": ("subtype", "Endometrial Carcinoma"),
    "ESCA": ("subtype", "Esophageal Squamous Cell Carcinoma"), "CESC": ("subtype", "Cervical Squamous Cell Carcinoma"),
    "PRAD": ("subtype", "Prostate Adenocarcinoma"),
    "BRCA": ("lineage", "Breast"), "HNSC": ("lineage", "Head and Neck"), "SARC": ("lineage", "Soft Tissue"),
    "THCA": ("lineage", "Thyroid"), "KIRP": ("lineage", "Kidney"), "MESO": ("lineage", "Pleura"),
    "TGCT": ("lineage", "Testis"),
}


def main():
    dep = json.loads((DD / "manifest.json").read_text())["divisions"]
    tcga = json.loads((DD / "tcga" / "manifest.json").read_text())["divisions"]
    tnames = json.loads((DD / "tcga" / "type_desc.json").read_text())
    tgroups = tcga["type"]["groups"]

    out = []
    for code, (div, group) in CROSSWALK.items():
        if code not in tgroups:
            print(f"  skip {code}: not a TCGA type here"); continue
        if group not in dep[div]["groups"]:
            print(f"  skip {code}: DepMap {div} group {group!r} missing"); continue
        out.append({"tcga": code, "tcga_label": f"{code} · {tnames.get(code, code)}",
                    "div": div, "group": group, "n_tcga": tgroups[code], "n_dep": dep[div]["groups"][group]})
    out.sort(key=lambda c: c["tcga_label"])
    (DD / "crosswalk.json").write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote crosswalk.json: {len(out)} cancers with a DepMap match")


if __name__ == "__main__":
    main()
