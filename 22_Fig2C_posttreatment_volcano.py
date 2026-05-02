"""Fig 2C — Post-treatment macrophage DGE volcano under DOR."""
import _00_common as cm
from _03_volcano_helper import make_volcano

n_R, n_NR = make_volcano(
    cm.OUT_DIR / "macrophage_posttreatment_DGE_LMM_DOR.csv",
    "Post-treatment macrophages",
    "Fig2C")
print(f"wrote Fig2C  {n_R + n_NR} DEGs (up_R={n_R}, up_NR={n_NR})")
