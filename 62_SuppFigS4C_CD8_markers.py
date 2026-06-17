"""Supp Fig S4C — CD8 activation/exhaustion marker dot plot.

Per-cell expression of canonical CD8 activation markers (GZMA, GZMB, GZMK,
PRF1, IFNG, TNF, IL2) and exhaustion-associated markers (PDCD1, LAG3,
HAVCR2, TIGIT, CTLA4, TOX, EOMES) across NivoRd CD138⁻ scRNA-seq CD8 T cells
(Naive + GZMK⁺ Effector + GZMB⁺ Effector subsets pooled; 5 paired patients).

Dot SIZE  = % of cells with non-zero expression of the gene (detection rate).
Dot COLOR = mean log-normalised expression across all cells in the group.

Detection-aware visualisation (a Cohen's d heatmap would be misleading because
most exhaustion-associated transcripts are detected in <15% of CD8 T cells in
this BM cohort; this dot plot makes that detection-rate context explicit).
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

import _00_common as cm
from _01_data_loaders import load_scrna_immune

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})

ACTIVATION = ["GZMA", "GZMB", "GZMK", "PRF1", "IFNG", "TNF", "IL2"]
EXHAUSTION = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "EOMES"]
CD8_REFINED = {"CD8 Naive", "CD8 GZMK+ Effector", "CD8 GZMB+ Effector"}
GROUPS = [
    ("DR",   "Baseline", "DR pre"),
    ("DR",   "PostNivo", "DR post"),
    ("P-NR", "Baseline", "P/NR pre"),
    ("P-NR", "PostNivo", "P/NR post"),
]


def _compute(cd8):
    var_names = list(cd8.var_names)
    genes = [g for g in ACTIVATION + EXHAUSTION if g in var_names]
    X = cd8.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    group_keys = (cd8.obs["DOR_group"].astype(str) + "_"
                  + cd8.obs["Tp"].astype(str)).values
    rows = []
    for g in genes:
        gi = var_names.index(g)
        for dor, tp, label in GROUPS:
            mask = group_keys == f"{dor}_{tp}"
            x = X[mask, gi] if mask.sum() else np.array([0.0])
            rows.append({"gene": g, "group": label,
                         "pct_exp": 100.0 * (x > 0).sum() / len(x),
                         "mean_exp": float(x.mean())})
    return pd.DataFrame(rows), genes


a = load_scrna_immune(nivord_only=True, exclude_doublets=True)
cd8 = a[a.obs["celltype_refined"].astype(str).isin(CD8_REFINED)].copy()
cd8 = cd8[cd8.obs["Pt"].astype(str).isin(cm.SCRNA_COHORT)].copy()

summary, genes = _compute(cd8)
n_act = sum(1 for g in ACTIVATION if g in genes)

norm = Normalize(vmin=0.0, vmax=float(summary["mean_exp"].max()))
cmap = plt.get_cmap("Reds")
size_scale = 14   # dot area (points²) per 1% expressing

fig = plt.figure(figsize=(13, 5.0))
gs = fig.add_gridspec(1, 2, width_ratios=[6, 1], wspace=0.04)
ax = fig.add_subplot(gs[0, 0])
legend_ax = fig.add_subplot(gs[0, 1]); legend_ax.set_axis_off()

for gi, g in enumerate(genes):
    for yi, (_, _, label) in enumerate(GROUPS):
        row = summary[(summary["gene"] == g) & (summary["group"] == label)].iloc[0]
        ax.scatter(gi, yi, s=max(2, row["pct_exp"] * size_scale),
                   color=cmap(norm(row["mean_exp"])),
                   edgecolor="black", linewidth=0.4, zorder=3)

ax.axvline(n_act - 0.5, color="grey", lw=0.6, ls="--", alpha=0.6)
ax.set_xticks(range(len(genes)))
ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(len(GROUPS)))
ax.set_yticklabels([lbl for _, _, lbl in GROUPS], fontsize=9)
ax.invert_yaxis()
ax.set_xlim(-0.7, len(genes) - 0.3)
ax.set_ylim(len(GROUPS) - 0.5, -0.8)
ax.text((n_act - 1) / 2, -0.75, "Activation",
        ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.text((n_act + len(genes) - 1) / 2, -0.75, "Exhaustion",
        ha="center", va="bottom", fontsize=10, fontweight="bold")
ax.tick_params(axis="both", which="both", length=0)
ax.set_axisbelow(True)
ax.grid(True, axis="both", color="lightgrey", lw=0.4, alpha=0.5)
for spine in ax.spines.values():
    spine.set_visible(False)

# Size legend
for k, p in enumerate([1, 10, 25, 50]):
    y = 0.78 - k * 0.10
    legend_ax.scatter(0.10, y, s=max(2, p * size_scale),
                      color="lightgrey", edgecolor="black", linewidth=0.4,
                      transform=legend_ax.transAxes, clip_on=False)
    legend_ax.text(0.25, y, f"{p}%",
                   transform=legend_ax.transAxes, fontsize=8, va="center")
legend_ax.text(0.05, 0.88, "% cells\nexpressing",
               transform=legend_ax.transAxes, fontsize=9, fontweight="bold",
               va="bottom")

# Colorbar
cax = legend_ax.inset_axes([0.60, 0.20, 0.15, 0.60])
sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("Mean log-norm\nexpression", fontsize=8)
cbar.ax.tick_params(labelsize=7)

cm.save_fig(fig, "SuppFigS4C")
plt.close(fig)
print("wrote SuppFigS4C — CD8 activation/exhaustion marker dot plot")
