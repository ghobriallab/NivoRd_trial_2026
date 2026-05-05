"""Fig 2C — Post-treatment macrophage DGE volcano under DOR."""
import _00_common as cm
from _03_volcano_helper import make_volcano

n_DR, n_PNR = make_volcano(
    cm.OUT_DIR / "macrophage_posttreatment_DGE_LMM_DOR.csv",
    "Post-treatment macrophages",
    "Fig2C")
print(f"wrote Fig2C  {n_DR + n_PNR} DEGs (up_DR={n_DR}, up_PNR={n_PNR})")
