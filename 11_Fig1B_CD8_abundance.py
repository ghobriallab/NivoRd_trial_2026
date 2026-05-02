"""Fig 1B — CD8+ T cell abundance under DOR."""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _00_common as cm
from _01_data_loaders import load_imc_extras
from _02_boxplot import boxplot_panel, save_legend

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

df = load_imc_extras().rename(columns={"pct_T8": "value"})
fig, ax = plt.subplots(figsize=(2.4, 3.2))
boxplot_panel(ax, df, value_col="value",
              ylabel="CD8⁺ T cells (% of classified)",
              title="CD8⁺ abundance",
              classifier=cm.DOR)
fig.tight_layout()
cm.save_fig(fig, "Fig1B")
plt.close(fig)
save_legend(str(cm.FIG_DIR / "Fig1_legend"))
print("wrote Fig1B")
