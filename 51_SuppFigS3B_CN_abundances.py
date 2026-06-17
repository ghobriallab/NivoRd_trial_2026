"""Supp Fig S3B — Cellular neighborhood abundances (mega-panel, 2 sub-panels).

Purpose:      Render Supp Fig S3B: per-sample CN abundance summary (CN4, CN2).
Inputs:       pheno_8clus-sample_pct.csv.
Outputs:      figures/SuppFigS3B.pdf / .svg.
Dependencies: matplotlib, _00_common, _01_data_loaders, _02_boxplot.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _00_common as cm
from _01_data_loaders import load_imc_sample_pct
from _02_boxplot import boxplot_panel

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

PANELS = [
    ("CN4", "CN4 (% cells)", "CN4 — PC-rich\n(77% PC)"),
    ("CN2", "CN2 (% cells)", "CN2 — PC-Mac contact\n(30% PC + 24% Mac)"),
]

cn = load_imc_sample_pct()

fig, axes = plt.subplots(1, 2, figsize=(5.6, 3.6))
for ax, (col, ylab, title) in zip(axes, PANELS):
    agg = (cn.groupby(["Pt", "Tp"], as_index=False)[col].mean()
              .rename(columns={col: "value"}))
    boxplot_panel(ax, agg, "value",
                  ylabel=ylab, title=title, classifier=cm.DOR)

fig.tight_layout()
cm.save_fig(fig, "SuppFigS3B")
plt.close(fig)
print("wrote SuppFigS3B (2 sub-panels: CN4, CN2)")
