"""Fig 2B — Pre-treatment macrophage DGE volcano under DOR."""
import _00_common as cm
from _03_volcano_helper import make_volcano

n_R, n_NR = make_volcano(
    cm.OUT_DIR / "macrophage_pretreatment_DGE_LMM_DOR.csv",
    "Pre-treatment macrophages",
    "Fig2B")
print(f"wrote Fig2B  {n_R + n_NR} DEGs (up_R={n_R}, up_NR={n_NR})")
