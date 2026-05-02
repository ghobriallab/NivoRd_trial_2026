"""Supp Fig S2C — Mean IMC marker expression by lineage class (heatmap).

Regenerated from the IMC h5ad (`IMC-import.h5ad`) with phenotype
assignments from `phenotyped_canonical.csv`. Mean expression is computed
per (cell_type, marker) and min-max scaled per column so each marker's
maximum cell-type signal = 1.0, matching scanpy's `matrixplot` "scale"
output convention used in the published version.

Markers shown match the published panel:
  CD163, CD11b, CD31, CD45, CD4, CD68, CD20, CD8a, CD56, CD138, CD3, PDL1
Cell types: B, MK, Mac, Myeloid, PC, T4, T8 (Unclassified excluded).

Output: figures/SuppFigS2C.{pdf,svg}
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

import _common as cm
from _data_loaders import load_imc

warnings.filterwarnings("ignore")
mpl.rcParams.update({"font.family": "Arial", "font.size": 10,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

MARKERS = ["CD163", "CD11b", "CD31", "CD45", "CD4", "CD68",
           "CD20", "CD8a", "CD56", "CD138", "CD3", "PDL1"]
CELL_TYPES = ["B", "MK", "Mac", "Myeloid", "PC", "T4", "T8"]
LABELS_X = ["CD163", "CD11b", "CD31", "CD45", "CD4", "CD68",
            "CD20", "CD8a", "CD56", "CD138", "CD3", "PD-L1"]


def main():
    a = load_imc(attach_pheno=True)
    obs = a.obs.copy()
    obs["pheno"] = obs["pheno"].astype(str)
    keep = obs["pheno"].isin(CELL_TYPES)
    a = a[keep.values].copy()

    var_idx = {m: list(a.var_names).index(m) for m in MARKERS}
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)

    rows = []
    for ct in CELL_TYPES:
        m = (a.obs["pheno"].astype(str).values == ct)
        if m.sum() == 0:
            continue
        means = {mk: float(X[m, var_idx[mk]].mean()) for mk in MARKERS}
        rows.append({"cell_type": ct, **means})
    df = pd.DataFrame(rows).set_index("cell_type")[MARKERS]

    # Min-max scale each column (marker) to [0, 1] across cell types.
    df_scaled = df.copy()
    for c in df_scaled.columns:
        col = df_scaled[c]
        df_scaled[c] = (col - col.min()) / (col.max() - col.min() + 1e-12)

    # White -> brick-red colormap matching the published scanpy matrixplot
    # (peak red ~#cb181d, not the very dark #67000d that matplotlib's
    # built-in 'Reds' saturates to).
    cmap = LinearSegmentedColormap.from_list(
        "matrixplot_red", ["#ffffff", "#fcbba1", "#fb6a4a", "#cb181d"], N=256)

    fig, ax = plt.subplots(figsize=(7.0, 3.6))
    im = ax.imshow(df_scaled.values, aspect="auto",
                    cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(np.arange(len(MARKERS)))
    ax.set_xticklabels(LABELS_X, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(np.arange(len(CELL_TYPES)))
    ax.set_yticklabels(CELL_TYPES, fontsize=10)
    # Cell borders
    ax.set_xticks(np.arange(len(MARKERS) + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(CELL_TYPES) + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="0.85", linewidth=0.6)
    ax.tick_params(which="minor", length=0)
    ax.spines[["top", "right", "bottom", "left"]].set_visible(True)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02,
                         ticks=[0, 0.5, 1])
    cbar.set_label("Mean expression\nin group (scaled)",
                    fontsize=9, labelpad=8)
    cbar.ax.tick_params(labelsize=9)

    fig.tight_layout()
    cm.save_fig(fig, "SuppFigS2C")
    plt.close(fig)
    df.to_csv(cm.OUT_DIR / "SuppFigS2C_marker_means_raw.csv")
    df_scaled.to_csv(cm.OUT_DIR / "SuppFigS2C_marker_means_scaled.csv")
    print(f"wrote SuppFigS2C ({len(df)} cell types × {len(df.columns)} markers)")


if __name__ == "__main__":
    main()
