"""Fig 2B — Pre-treatment macrophage DGE volcano under DOR.

Purpose:      Render Figure 2B: Macrophage DGE volcano, pre-treatment DR vs P-NR.
Inputs:       macrophage_pretreatment_DGE_LMM_DOR.csv.
Outputs:      figures/Fig2B.pdf / .svg.
Dependencies: matplotlib, pandas, _00_common, _03_volcano_helper.
"""
import _00_common as cm
from _03_volcano_helper import make_volcano

n_DR, n_PNR = make_volcano(
    cm.OUT_DIR / "macrophage_pretreatment_DGE_LMM_DOR.csv",
    "Pre-treatment macrophages",
    "Fig2B")
print(f"wrote Fig2B  {n_DR + n_PNR} DEGs (up_DR={n_DR}, up_PNR={n_PNR})")
