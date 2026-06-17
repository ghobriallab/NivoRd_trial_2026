"""Supp Fig S4A — M1 macrophage module score (DR vs P/NR × pre/post).

Purpose:      Render Supp Fig S4A: M1 macrophage module score (Hallmark classical activation).
Inputs:       CD138neg_immune_cells.h5ad (macrophage subset).
Outputs:      figures/SuppFigS4A.pdf / .svg.
Dependencies: scanpy, matplotlib, pandas, _00_common, _01_data_loaders, _02_boxplot.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import scanpy as sc

import _00_common as cm
from _01_data_loaders import load_scrna_immune
from _02_boxplot import boxplot_panel

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

M1_GENES = ["CD86", "HLA-DRA", "IL1B", "TNF", "NOS2",
            "IDO1", "CXCL10", "IL6", "IL12B"]

a = load_scrna_immune(nivord_only=True, exclude_doublets=True)
macs = a[a.obs["celltype_refined"].astype(str) == "Macrophages"].copy()
macs = macs[macs.obs["Pt"].astype(str).isin(cm.SCRNA_COHORT)].copy()

avail = [g for g in M1_GENES if g in macs.var_names]
sc.tl.score_genes(macs, gene_list=avail, score_name="M1_score",
                  random_state=0, n_bins=25)

rows = []
for (pt, tp), idx in macs.obs.groupby(["Pt", "Tp"]).groups.items():
    rows.append({"Pt": pt, "Tp": tp,
                 "value": float(macs.obs.loc[idx, "M1_score"].mean())})
df = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(2.6, 3.4))
boxplot_panel(ax, df, value_col="value",
              ylabel="M1 module score",
              title="M1 module (classical / IFN-γ-induced)",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "SuppFigS4A")
plt.close(fig)
print("wrote SuppFigS4A — M1 macrophage module")
