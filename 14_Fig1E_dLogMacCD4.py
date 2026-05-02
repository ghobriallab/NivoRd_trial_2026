"""Fig 1E — paired Δlog(Mac/CD4) post − pre vs ordinal DOR.

Single pre-planned ordinal Spearman test for the IMC lineage axis.
n=7 patients with paired baseline+post-tx IMC (Pt03 baseline missing).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

import _00_common as cm
from _01_data_loaders import load_imc_extras

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

df = load_imc_extras().copy()
# Mac/CD4 ratio per patient × tp; log space; pair post − pre
df["mac_cd4"] = df["pct_Mac"] / df["pct_T4"]
piv = df.pivot_table(index="Pt", columns="Tp", values="mac_cd4",
                       aggfunc="first")
piv["delta_log"] = np.log(piv["PostNivo"]) - np.log(piv["Baseline"])
piv["DOR"] = piv.index.map(cm.DOR_ORDINAL)
piv["DOR_group"] = piv.index.map(cm.DOR)
piv = piv.dropna(subset=["delta_log", "DOR"]).reset_index()

rho, p, n = cm.spearman(piv["delta_log"], piv["DOR"])

fig, ax = plt.subplots(figsize=(3.8, 3.4))
texts = []
for _, r in piv.iterrows():
    color = cm.COLORS["R_PostNivo"] if r["DOR_group"] == "R" else cm.COLORS["NR_PostNivo"]
    ax.scatter(r["delta_log"], r["DOR"],
               s=110, alpha=0.9, edgecolor="black", linewidths=0.7,
               facecolor=color, zorder=3)
    texts.append(ax.text(r["delta_log"], r["DOR"], r["Pt"], fontsize=8))

slope, intercept = np.polyfit(piv["delta_log"], piv["DOR"], 1)
xs = np.linspace(piv["delta_log"].min() - 0.05,
                  piv["delta_log"].max() + 0.05, 50)
ax.plot(xs, slope * xs + intercept, c="black", lw=1.0, ls="--", alpha=0.7)
ax.axvline(0, c="grey", lw=0.5, ls=":")

ax.set_xlabel(r"$\Delta$ log(Mac/CD4)  post $-$ pre", fontsize=10)
ax.set_ylabel("Response (1=SD ... 4=CR)", fontsize=10)
ax.set_yticks([1, 2, 3, 4])
ax.set_yticklabels(["SD", "PR", "VGPR", "CR"], fontsize=9)
ax.set_title(rf"Spearman $\rho={rho:+.2f}$, $p={p:.4f}$",
              fontsize=10, fontweight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

adjust_text(texts, ax=ax, expand=(1.4, 1.6),
             arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5))
fig.tight_layout()
cm.save_fig(fig, "Fig1E")
plt.close(fig)

# Save the per-patient table for traceability
piv.to_csv(cm.OUT_DIR / "Fig1E_dLogMacCD4_perpt.csv", index=False)
print(f"wrote Fig1E  rho={rho:+.3f} p={p:.4f} n={n}")
