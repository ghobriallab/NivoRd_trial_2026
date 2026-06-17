"""Supp Fig S4B — M2 macrophage module score (DR vs P/NR × pre/post).

Per-patient mean "alternative activation" (M2-like) module score on the
macrophage compartment of the NivoRd CD138⁻ scRNA-seq cohort
(celltype_refined == "Macrophages"; 5 paired patients).

M2 gene list:
    CD163, MRC1, ARG1, IL10, MSR1, MERTK, TGFB1, CCL22

Plotted via the canonical _02_boxplot.boxplot_panel helper so the style
matches Fig 1B-C and Supp Figs S3A-C.
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

M2_GENES = ["CD163", "MRC1", "ARG1", "IL10", "MSR1",
            "MERTK", "TGFB1", "CCL22"]

a = load_scrna_immune(nivord_only=True, exclude_doublets=True)
macs = a[a.obs["celltype_refined"].astype(str) == "Macrophages"].copy()
macs = macs[macs.obs["Pt"].astype(str).isin(cm.SCRNA_COHORT)].copy()

avail = [g for g in M2_GENES if g in macs.var_names]
sc.tl.score_genes(macs, gene_list=avail, score_name="M2_score",
                  random_state=0, n_bins=25)

rows = []
for (pt, tp), idx in macs.obs.groupby(["Pt", "Tp"]).groups.items():
    rows.append({"Pt": pt, "Tp": tp,
                 "value": float(macs.obs.loc[idx, "M2_score"].mean())})
df = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(2.6, 3.4))
boxplot_panel(ax, df, value_col="value",
              ylabel="M2 module score",
              title="M2 module (alternative / IL-10–TGFβ)",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "SuppFigS4B")
plt.close(fig)
print("wrote SuppFigS4B — M2 macrophage module")
