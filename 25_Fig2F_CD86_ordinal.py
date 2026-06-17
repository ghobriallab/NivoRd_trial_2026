"""Fig 2F — CD86 expression Δ post − pre vs ordinal DOR.

Purpose:      Render Figure 2F: post-pre delta-CD86 expression vs DOR_ordinal (Spearman).
Inputs:       CD138neg_immune_cells.h5ad (macrophages subset).
Outputs:      figures/Fig2F.pdf / .svg.
Dependencies: scanpy, matplotlib, scipy.stats, _00_common, _01_data_loaders.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text

import _00_common as cm
from _01_data_loaders import load_scrna_immune

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

a = load_scrna_immune(nivord_only=True)
mac = a[a.obs["Immune_cell_2"].astype(str) == "CD163+ Macrophages"].copy()
print(f"  Mac cells (broad def): {mac.n_obs}; "
      f"patients: {sorted(mac.obs['Pt'].dropna().unique())}", flush=True)

if "CD86" not in mac.var_names:
    raise SystemExit("CD86 not in var_names")
expr = mac[:, "CD86"].X
expr = expr.toarray().ravel() if hasattr(expr, "toarray") else np.asarray(expr).ravel()
mac.obs = mac.obs.copy()
mac.obs["CD86"] = expr

cell_count = mac.obs.groupby(["Pt", "Tp"]).size().rename("n").reset_index()
mean_expr = (mac.obs.groupby(["Pt", "DOR_ordinal", "Tp"], as_index=False)
                       ["CD86"].mean())
mean_expr = mean_expr.merge(cell_count, on=["Pt", "Tp"])
mean_expr = mean_expr[mean_expr["n"] >= 5].copy()

pre = mean_expr[mean_expr["Tp"] == "Baseline"].set_index("Pt")
post = mean_expr[mean_expr["Tp"] == "PostNivo"].set_index("Pt")
common = sorted(set(pre.index) & set(post.index))
delta = pd.DataFrame({
    "Pt": common,
    "DOR_ordinal": [pre.loc[p, "DOR_ordinal"] for p in common],
    "delta": [post.loc[p, "CD86"] - pre.loc[p, "CD86"] for p in common],
})
delta["DOR_group"] = delta["Pt"].map(cm.DOR)
rho, p, n = cm.spearman(delta["delta"], delta["DOR_ordinal"])
print(f"  CD86 Δ vs DOR: rho={rho:+.3f}, p={p:.4f}", flush=True)

fig, ax = plt.subplots(figsize=(3.8, 3.4))
texts = []
for _, r in delta.iterrows():
    color = cm.COLORS["DR_PostNivo"] if r["DOR_group"] == "DR" else cm.COLORS["P-NR_PostNivo"]
    ax.scatter(r["delta"], r["DOR_ordinal"],
                s=110, alpha=0.9, edgecolor="black", linewidths=0.7,
                facecolor=color, zorder=3)
    texts.append(ax.text(r["delta"], r["DOR_ordinal"], r["Pt"], fontsize=8))

slope, intercept = np.polyfit(delta["delta"], delta["DOR_ordinal"], 1)
xs = np.linspace(delta["delta"].min() - 0.05,
                  delta["delta"].max() + 0.05, 50)
ax.plot(xs, slope * xs + intercept, c="black", lw=1.0, ls="--", alpha=0.7)
ax.axvline(0, c="grey", lw=0.5, ls=":")

ax.set_xlabel(r"$\Delta$ CD86 expression  post $-$ pre", fontsize=10)
ax.set_ylabel("Response (1=SD ... 4=CR)", fontsize=10)
ax.set_yticks([1, 2, 3, 4])
ax.set_yticklabels(["SD", "PR", "VGPR", "CR"], fontsize=9)
ax.set_title(rf"Spearman $\rho={rho:+.2f}$, $p={p:.4f}$",
              fontsize=10, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
adjust_text(texts, ax=ax, expand=(1.4, 1.6),
             arrowprops=dict(arrowstyle="-", color="0.5", lw=0.5))
fig.tight_layout()
cm.save_fig(fig, "Fig2F")
plt.close(fig)
delta.to_csv(cm.OUT_DIR / "Fig2F_CD86delta_perpt.csv", index=False)
print("wrote Fig2F")
