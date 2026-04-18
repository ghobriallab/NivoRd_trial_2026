#!/usr/bin/env python3
"""Supplementary Figure S2C: IMC lineage marker expression across 7 cell types."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

INPUT_OBS = "phenotyped-obs.csv"
INPUT_MARKERS = "combat_all_markers.csv"
OUTPUT = "FigureS2C_imc_marker_expression"

LINEAGES = ["T4", "T8", "B", "Mac", "Myeloid", "MK", "PC"]
LINEAGE_LABELS = ["CD4+ T", "CD8+ T", "B", "Macrophage", "Myeloid", "Megakaryocyte", "Plasma Cell"]

obs = pd.read_csv(INPUT_OBS)
markers = pd.read_csv(INPUT_MARKERS)

merge_key = obs.columns.intersection(markers.columns).tolist()
if merge_key:
    dat = obs.merge(markers, on=merge_key, how="left")
else:
    dat = pd.concat([obs.reset_index(drop=True), markers.reset_index(drop=True)], axis=1)

marker_cols = [c for c in markers.columns if c not in merge_key and c != "Unnamed: 0"]
dat = dat[dat.pheno.isin(LINEAGES)]

mean_expr = dat.groupby("pheno")[marker_cols].mean().reindex(LINEAGES)
zscore = mean_expr.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=0)

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(zscore.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(len(marker_cols)))
ax.set_xticklabels(marker_cols, fontsize=9, rotation=45, ha="right")
ax.set_yticks(range(len(LINEAGES)))
ax.set_yticklabels(LINEAGE_LABELS, fontsize=11, fontweight="bold")
ax.set_title("IMC Marker Expression (Z-score)", fontsize=13, fontweight="bold", pad=10)

cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Z-score", fontsize=11, fontweight="bold")

plt.tight_layout()
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
