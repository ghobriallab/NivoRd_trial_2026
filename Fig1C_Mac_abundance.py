"""Fig 1C — Macrophage abundance (% lineage+) by IMC."""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _common as cm
from _data_loaders import load_imc_extras
from _boxplot import boxplot_panel

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
