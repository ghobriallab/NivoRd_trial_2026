"""Fig 2D — Pre-treatment macrophage Hallmark GSEA under DOR.

Purpose:      Render Figure 2D: Hallmark GSEA preranked on pre-treatment macrophage DGE.
Inputs:       macrophage_pretreatment_DGE_LMM_DOR.csv + MSigDB_Hallmark_2020.gmt.
Outputs:      figures/Fig2D.pdf / .svg + GSEA result CSV in OUT_DIR.
Dependencies: gseapy, matplotlib, pandas, _00_common, _04_gsea_helper.
"""
import _00_common as cm
from _04_gsea_helper import make_gsea_panel

make_gsea_panel(
    cm.OUT_DIR / "macrophage_pretreatment_DGE_LMM_DOR.csv",
    cm.OUT_DIR / "macrophage_pretreatment_GSEA_DOR.csv",
    "Pre-treatment macrophage GSEA",
    "Fig2D")
print("wrote Fig2D")
