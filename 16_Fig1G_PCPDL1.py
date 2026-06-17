"""Fig 1G — PD-L1+ plasma cells (% of plasma cells) by DR/P-NR × pre/post.

Purpose:      Render Figure 1G: Plasma-cell PD-L1+ % DR vs P-NR x pre/post.
Inputs:       imc_extras_persample.csv.
Outputs:      figures/Fig1G.pdf / .svg.
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

df = load_imc_extras().rename(columns={"pct_pdl1_in_PC": "value"})
fig, ax = plt.subplots(figsize=(2.6, 3.4))
boxplot_panel(ax, df, "value",
              ylabel="PD-L1+ plasma cells (%)",
              title="Plasma-cell PD-L1",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "Fig1G")
plt.close(fig)
print("wrote Fig1G")
