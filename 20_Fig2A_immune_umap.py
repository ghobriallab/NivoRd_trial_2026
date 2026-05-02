"""Fig 2A — CD138− BM immune cell UMAP, colored by cell type.

Adapted from 2-Scripts/03-Figure2A_immune_umap.py.
"""
from __future__ import annotations
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import anndata as ad
from scipy.stats import gaussian_kde
from adjustText import adjust_text

import _00_common as cm

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

CELLTYPE_COL = "celltype_refined"
EXCLUDE = ["Doublets", "Proliferating", "Hemoglobin", "ILC-like"]
RENAME = {"B": "B Cells", "DC": "Dendritic Cells"}
PALETTE = {
    "CD4 Naive": "#e41a1c", "CD4 Effector": "#ff7f00", "Tregs": "#b8860b",
    "CD8 Naive": "#6a3d9a", "CD8 GZMK+ Effector": "#e31a8f",
    "CD8 GZMB+ Effector": "#f781bf", "NK/NKT": "#1b9e77",
    "Classical Monocytes": "#33a02c", "Non-classical Monocytes": "#00bfc4",
    "Macrophages": "#1f4e9f", "B Cells": "#888888",
    "Dendritic Cells": "#d95f02",
}

a = ad.read_h5ad(cm.H5AD_IMMUNE)
obs = a.obs.copy()
umap = a.obsm["X_umap"].copy()
del a

mask = ~obs[CELLTYPE_COL].isin(EXCLUDE)
obs = obs[mask].copy()
umap = umap[mask.values]
obs[CELLTYPE_COL] = obs[CELLTYPE_COL].replace(RENAME)
counts = obs[CELLTYPE_COL].value_counts()
celltypes = counts[counts > 0].index.tolist()


def density_centroid(coords, percentile=70):
    if len(coords) < 50:
        return coords.mean(axis=0)
    idx = np.random.default_rng(42).choice(len(coords),
                                             min(2000, len(coords)),
                                             replace=False)
    kde = gaussian_kde(coords[idx].T)
    d = kde(coords.T)
    return coords[d >= np.percentile(d, percentile)].mean(axis=0)


fig, ax = plt.subplots(figsize=(5.5, 5))
for ct in reversed(celltypes):
    m = obs[CELLTYPE_COL] == ct
    ax.scatter(umap[m.values, 0], umap[m.values, 1],
                c=PALETTE.get(ct, "#888"), s=0.8, alpha=0.7,
                rasterized=True, linewidths=0)

ax.set_title(f"CD138$^-$ Bone Marrow Immune Cells (n={counts.sum():,})",
              fontsize=12, fontweight="bold")
ax.set_xlabel("UMAP 1", fontsize=11, fontweight="bold")
ax.set_ylabel("UMAP 2", fontsize=11, fontweight="bold")
ax.set_xticks([]); ax.set_yticks([])

texts = []
for ct in celltypes:
    m = obs[CELLTYPE_COL] == ct
    coords = umap[m.values]
    cx, cy = density_centroid(coords) if max(coords.std(0)) > 3.0 \
        else coords.mean(0)
    t = ax.text(cx, cy, ct, fontsize=10, fontweight="bold",
                 ha="center", va="center", color="black",
                 bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                            edgecolor=PALETTE.get(ct, "#888"),
                            alpha=0.95, linewidth=2))
    texts.append(t)

adjust_text(texts, ax=ax, expand=(1.05, 1.05),
             force_text=(0.3, 0.3), force_points=(0.8, 0.8),
             ensure_inside_axes=True,
             arrowprops=dict(arrowstyle="-", lw=0, alpha=0))

plt.tight_layout(pad=0.5)
cm.save_fig(fig, "Fig2A")
plt.close(fig)
print("wrote Fig2A")
