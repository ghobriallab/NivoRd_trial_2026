"""Fig 2E — Post-treatment macrophage Hallmark GSEA under DOR.

Purpose:      Render Figure 2E: Hallmark GSEA preranked on post-treatment macrophage DGE.
Inputs:       macrophage_posttreatment_DGE_LMM_DOR.csv + MSigDB_Hallmark_2020.gmt.
Outputs:      figures/Fig2E.pdf / .svg + GSEA result CSV in OUT_DIR.
Dependencies: gseapy, matplotlib, pandas, _00_common, _04_gsea_helper.
"""
import _00_common as cm
from _04_gsea_helper import make_gsea_panel

make_gsea_panel(
    cm.OUT_DIR / "macrophage_posttreatment_DGE_LMM_DOR.csv",
    cm.OUT_DIR / "macrophage_posttreatment_GSEA_DOR.csv",
    "Post-treatment macrophage GSEA",
    "Fig2E")
print("wrote Fig2E")
