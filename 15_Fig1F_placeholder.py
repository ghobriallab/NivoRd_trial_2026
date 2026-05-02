"""Fig 1F — Representative IMC composites (placeholder).

Cannot be regenerated without raw 34-channel IMC TIFFs.
Source script: `_archive/pipeline/1_canonical_IMC/fig1D_imc_composite.py`.
"""
from __future__ import annotations
import matplotlib.pyplot as plt
import matplotlib as mpl
import _00_common as cm

mpl.rcParams.update({"font.family": "Arial", "font.size": 9,
                      "pdf.fonttype": 42, "svg.fonttype": "none"})

fig, ax = plt.subplots(figsize=(7, 3.6))
ax.text(0.5, 0.55,
         "Fig 1F — representative IMC composites (placeholder)",
         ha="center", va="center", fontsize=12, fontweight="bold",
         transform=ax.transAxes)
ax.text(0.5, 0.30,
         "Source: 34-channel IMC TIFF directory (rakaia_images).\n"
         "See _archive/pipeline/1_canonical_IMC/fig1D_imc_composite.py.",
         ha="center", va="center", fontsize=9, style="italic", color="#555",
         transform=ax.transAxes)
ax.set_xticks([]); ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_color("#888"); spine.set_linewidth(0.5)
cm.save_fig(fig, "Fig1F")
plt.close(fig)
print("wrote Fig1F (placeholder)")
