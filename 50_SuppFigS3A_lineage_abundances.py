"""Supp Fig S3A — IMC lineage abundances (mega-panel, 4 sub-panels).

Purpose:      Render Supp Fig S3A: IMC lineage abundance mega-panel (CD4, PC, B, Myeloid).
Inputs:       imc_extras_persample.csv.
Outputs:      figures/SuppFigS3A.pdf / .svg.
Dependencies: matplotlib, _00_common, _01_data_loaders, _02_boxplot.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _00_common as cm
from _01_data_loaders import load_imc_extras
from _02_boxplot import boxplot_panel

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

PANELS = [
    ("pct_T4",      "CD4+ T cells (%)",  "CD4"),
    ("pct_PC",      "Plasma cells (%)",  "PC"),
    ("pct_B",       "B cells (%)",       "B"),
    ("pct_Myeloid", "Myeloid (%)",       "Myeloid"),
]

df = load_imc_extras()

fig, axes = plt.subplots(1, 4, figsize=(11.0, 3.4))
for ax, (col, ylab, title) in zip(axes, PANELS):
    boxplot_panel(ax, df.rename(columns={col: "value"}),
                  "value", ylabel=ylab, title=title, classifier=cm.DOR)

fig.tight_layout()
cm.save_fig(fig, "SuppFigS3A")
plt.close(fig)
print("wrote SuppFigS3A (4 sub-panels: CD4, PC, B, Myeloid)")
