"""Supp Fig S3B — Cellular neighborhood abundances (mega-panel, 2 sub-panels).

Sub-panels (left to right):
  (i)  CN4 — plasma-cell-rich (77% PC)
  (ii) CN2 — PC-Mac contact zone (30% PC + 24% Mac)

CN7 (macrophage-rich, 71% Mac) is shown in main Figure 1D, so omitted here.
CN3 (47% Myeloid + 37% Mac) is myeloid-rich and not biologically relevant
to the PC compartment, so omitted.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _common as cm
from _data_loaders import load_imc_sample_pct
from _boxplot import boxplot_panel

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
