"""Fig 1D — CN7 (macrophage-rich neighborhood, 71% Mac by composition) abundance."""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _common as cm
from _data_loaders import load_imc_sample_pct
from _boxplot import boxplot_panel

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

cn = load_imc_sample_pct()
agg = (cn.groupby(["Pt", "Tp"], as_index=False)["CN7"].mean()
          .rename(columns={"CN7": "value"}))
fig, ax = plt.subplots(figsize=(2.6, 3.4))
boxplot_panel(ax, agg, "value",
              ylabel="CN7 (% cells)",
              title="CN7 — Mac-rich neighborhood",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "Fig1D")
plt.close(fig)
print("wrote Fig1D")
