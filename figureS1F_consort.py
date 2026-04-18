#!/usr/bin/env python3
"""Supplementary Figure S1F: CONSORT flow diagram."""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

OUTPUT = "FigureS1F_CONSORT"

def draw_box(ax, x, y, w, h, text, fc="white", ec="black", fontsize=10, bold_first=False):
    ax.add_patch(mpatches.FancyBboxPatch((x - w/2, y - h/2), w, h,
                 boxstyle="round,pad=0.02", facecolor=fc, edgecolor=ec, linewidth=1.5))
    if bold_first:
        lines = text.split("\n")
        ax.text(x, y + h/2 - 0.04, lines[0], ha="center", va="top",
                fontsize=fontsize, fontweight="bold")
        if len(lines) > 1:
            rest = "\n".join(lines[1:])
            ax.text(x, y + h/2 - 0.08, rest, ha="center", va="top",
                    fontsize=fontsize - 1, color="#444444")
    else:
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, fontweight="bold")

def arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5))

def side_arrow(ax, x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.2))

fig, ax = plt.subplots(figsize=(8, 10))
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_axis_off()

# Boxes
draw_box(ax, 0.42, 0.90, 0.32, 0.07, "Consented to study\n(n = 15)")
draw_box(ax, 0.42, 0.72, 0.38, 0.07, "Eligible and registered\n(n = 10)")
draw_box(ax, 0.42, 0.54, 0.38, 0.07, "Received NivoRd\n(n = 8)")
draw_box(ax, 0.42, 0.30, 0.50, 0.12,
         "Analyzed (n = 8; 5R, 3NR)\nscRNA-seq: n = 5 (3R, 2NR)\n3 excluded: insufficient sample quality",
         fc="#e8f5e9", bold_first=True)

# Exclusion boxes (right side)
draw_box(ax, 0.82, 0.81, 0.30, 0.09,
         "Screen failures (n = 5)\nDid not meet criteria (n = 3)\nWithdrew consent (n = 2)",
         fc="#fde8e8", bold_first=True)
draw_box(ax, 0.82, 0.63, 0.30, 0.09,
         "Withdrew before treatment (n = 2)\nWithdrew consent (n = 1)\nPursued other therapy (n = 1)",
         fc="#fde8e8", bold_first=True)

# Arrows (vertical)
arrow(ax, 0.42, 0.865, 0.42, 0.755)
arrow(ax, 0.42, 0.685, 0.42, 0.575)
arrow(ax, 0.42, 0.505, 0.42, 0.36)

# Arrows (to exclusion boxes)
side_arrow(ax, 0.58, 0.81, 0.67, 0.81)
side_arrow(ax, 0.61, 0.63, 0.67, 0.63)

# Data lock annotation
ax.text(0.42, 0.16, "Data lock: 20 December 2024",
        ha="center", fontsize=10, fontstyle="italic", color="#555555")

plt.tight_layout(pad=0.5)
for fmt in ["pdf", "png", "svg"]:
    fig.savefig(f"{OUTPUT}.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
