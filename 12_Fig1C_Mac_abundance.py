"""Fig 1C — Macrophage abundance (% lineage+) by IMC.

Purpose:      Render Figure 1C: Macrophage abundance (% of classified) DR vs P-NR x pre/post.
Inputs:       imc_extras_persample.csv.
Outputs:      figures/Fig1C.pdf / .svg.
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

df = load_imc_extras().rename(columns={"pct_Mac": "value"})
fig, ax = plt.subplots(figsize=(2.6, 3.4))
boxplot_panel(ax, df, "value",
              ylabel="Macrophages (% of classified)",
              title="Macrophage abundance",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "Fig1C")
plt.close(fig)
print("wrote Fig1C")
