#!/usr/bin/env python3
"""Supplementary Figure S2F: Cytogenetic gene expression validation heatmap."""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd

import sys; sys.path.insert(0, "."); from _00_common import FIG_DIR, DATA_LOCAL
EXPR_CSV = str(DATA_LOCAL / "gene_expression_summary.csv")
VALID_CSV = str(DATA_LOCAL / "cytogenetic_validation.csv")
OUTPUT = str(FIG_DIR / "SuppFigS2F")

GENE_ORDER = ["RB1", "WHSC1", "CKS1B", "MCL1", "ANP32E", "CCND1", "CDKN2C", "FAF1"]
CASE_ORDER = ["Case #1", "Case #4", "Case #5", "Case #6", "Case #7"]

EXPECTED_PATTERNS = {
    "RB1":    {"high": [], "low": ["Case #4"]},
    "WHSC1":  {"high": ["Case #4", "Case #5"], "low": []},
    "CKS1B":  {"high": ["Case #5"], "low": []},
    "MCL1":   {"high": ["Case #5"], "low": []},
    "ANP32E": {"high": ["Case #5"], "low": []},
    "CCND1":  {"high": ["Case #6"], "low": []},
    "CDKN2C": {"high": [], "low": ["Case #7"]},
    "FAF1":   {"high": [], "low": ["Case #7"]},
}

def get_sig(p):
    if pd.isna(p): return ""
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

# Load
expr = pd.read_csv(EXPR_CSV)
valid = pd.read_csv(VALID_CSV)

mean_by_case = expr.pivot_table(index="Case", columns="Gene", values="Mean", aggfunc="mean")
if "Healthy Mean" in valid.columns:
    for gene in GENE_ORDER:
        hm = valid.loc[valid.Gene == gene, "Healthy Mean"]
        if len(hm): mean_by_case.loc["Healthy", gene] = hm.iloc[0]

available = [g for g in GENE_ORDER if g in mean_by_case.columns]
rows = ["Healthy"] + CASE_ORDER
mean_by_case = mean_by_case.reindex([r for r in rows if r in mean_by_case.index])[available]

zscore = mean_by_case.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=0)

pcol = "p-value" if "p-value" in valid.columns else "p_value"
sig_matrix = pd.DataFrame("", index=zscore.index, columns=available)
for _, row in valid.iterrows():
    if row.Gene in sig_matrix.columns and row.Case in sig_matrix.index:
        sig_matrix.loc[row.Case, row.Gene] = get_sig(row[pcol])

# Plot
fig, ax = plt.subplots(figsize=(8, 4.5))
im = ax.imshow(zscore.values, cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_xticks(range(len(available)))
ax.set_xticklabels(available, fontsize=12, fontweight="bold", rotation=45, ha="right")
ax.set_yticks(range(len(zscore.index)))
ax.set_yticklabels(zscore.index, fontsize=12, fontweight="bold")

for i, case in enumerate(zscore.index):
    for j, gene in enumerate(available):
        val = zscore.loc[case, gene]
        sig = sig_matrix.loc[case, gene] if case in sig_matrix.index else ""
        color = "white" if abs(val) > 1 else "black"
        text = f"{val:.1f}" + (f"\n{sig}" if sig and sig != "ns" else "")
        ax.text(j, i, text, ha="center", va="center", color=color, fontsize=10, fontweight="bold")
        info = EXPECTED_PATTERNS.get(gene, {})
        if case in info.get("high", []):
            ax.add_patch(patches.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor="red", lw=2.5))
        elif case in info.get("low", []):
            ax.add_patch(patches.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor="blue", lw=2.5))

cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
cbar.set_label("Z-score", fontsize=12, fontweight="bold")
ax.set_title("Cytogenetic Gene Expression Validation", fontsize=14, fontweight="bold", pad=12)
fig.text(0.5, 0.01,
         "* p<0.05  ** p<0.01  *** p<0.001 vs Healthy  |  Red = expected high  Blue = expected low",
         ha="center", va="bottom", fontsize=9, style="italic")

plt.tight_layout(rect=[0, 0.04, 1, 1])
for fmt in ["pdf", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
