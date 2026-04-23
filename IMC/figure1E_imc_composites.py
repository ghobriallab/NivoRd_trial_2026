#!/usr/bin/env python3
"""Figure 1E: IMC representative composites (Baseline vs Post-Nivo) for patients 1/4/6/9.

Usage:
    python figure1E.py --input path/to/rakaia_images --output path/to/out

Inputs: 34-channel IMC TIFFs named
    PatientN_Baseline_XXX-PatientN-Baseline.tiff
    PatientN_PostNivo_XXX-PatientN-PostNivo.tiff
Channel indices match the `Labels` list embedded in each TIFF's ImageJ IJMetadata.
"""

import argparse
import glob
import os

import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont

CHANNELS = {"Histone": 33, "CD138": 23, "CD68": 16, "CD8a": 19}
COLORS = {
    "Histone": (0.55, 0.55, 0.55),
    "CD138":   (0.95, 0.20, 0.85),
    "CD68":    (1.00, 0.80, 0.10),
    "CD8a":    (0.10, 0.90, 1.00),
}

# ROI pair per patient (largest shared crop with representative %CD8a+)
PAIRS = {
    1: ("Patient1_Baseline_001-Patient1-Baseline.tiff",
        "Patient1_PostNivo_001-Patient1-PostNivo.tiff"),
    4: ("Patient4_Baseline_003-Patient4-Baseline.tiff",
        "Patient4_PostNivo_002-Patient4-PostNivo.tiff"),
    6: ("Patient6_Baseline_001-Patient6-Baseline.tiff",
        "Patient6_PostNivo_006-Patient6-PostNivo.tiff"),
    9: ("Patient9_Baseline_002-Patient9-Baseline.tiff",
        "Patient9_PostNivo_002-Patient9-PostNivo.tiff"),
}

UPSCALE = 2
SCALE_BAR_UM = 100
UM_PER_PX = 1.0  # IMC Hyperion native resolution
PERCENTILES = (1, 99.5)

FONT_CANDIDATES = (
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    "C:/Windows/Fonts/arialbd.ttf",
)


def load_font(size):
    for path in FONT_CANDIDATES:
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            continue
    return ImageFont.load_default()


def shared_bounds(stacks):
    out = {}
    for name, idx in CHANNELS.items():
        pool = np.concatenate([s[idx].ravel() for s in stacks])
        lo, hi = np.percentile(pool, PERCENTILES)
        out[name] = (float(lo), float(max(hi, lo + 1)))
    return out


def composite(stack, bounds):
    H, W = stack.shape[1:]
    rgb = np.zeros((H, W, 3), np.float32)
    for name, idx in CHANNELS.items():
        lo, hi = bounds[name]
        x = np.clip((stack[idx] - lo) / (hi - lo), 0, 1)
        rgb += x[..., None] * np.array(COLORS[name], np.float32)[None, None, :]
    return np.clip(rgb, 0, 1)


def center_crop(stack, h, w):
    H, W = stack.shape[1:]
    y0 = (H - h) // 2
    x0 = (W - w) // 2
    return stack[:, y0:y0 + h, x0:x0 + w]


def add_scale_bar(img, length_um, um_per_px, upscale):
    bar_px = int(length_um / um_per_px) * upscale
    bar_h = max(10, upscale * 5)
    margin = 24
    label = f"{length_um} μm"

    lo, hi, best = 10, 400, None
    while lo <= hi:
        mid = (lo + hi) // 2
        font = load_font(mid)
        bb = font.getbbox(label)
        tw, th = bb[2] - bb[0], bb[3] - bb[1]
        if tw <= bar_px:
            best = (font, tw, th)
            lo = mid + 1
        else:
            hi = mid - 1
    font, tw, th = best

    x1 = img.width - margin
    x0 = x1 - bar_px
    y1 = img.height - margin
    y0 = y1 - bar_h
    d = ImageDraw.Draw(img)
    d.rectangle([x0 - 1, y0 - 1, x1 + 1, y1 + 1], fill="black")
    d.rectangle([x0, y0, x1, y1], fill="white")

    tx = x0 + (bar_px - tw) // 2
    ty = y0 - th - 10
    for dx in (-2, -1, 0, 1, 2):
        for dy in (-2, -1, 0, 1, 2):
            d.text((tx + dx, ty + dy), label, font=font, fill="black")
    d.text((tx, ty), label, font=font, fill="white")


def render_panel(stack, bounds, path):
    arr = composite(stack, bounds)
    img = Image.fromarray((arr * 255).astype(np.uint8)).convert("RGB")
    img = img.resize((img.width * UPSCALE, img.height * UPSCALE), Image.BILINEAR)
    add_scale_bar(img, SCALE_BAR_UM, UM_PER_PX, UPSCALE)
    img.save(path)


def save_legend(stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig = plt.figure(figsize=(2.0, 2.4), dpi=600)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.08, 0.92, "Legend", fontsize=16, fontweight="bold", va="top")
    for i, (name, col) in enumerate(COLORS.items()):
        y = 0.73 - i * 0.16
        ax.add_patch(Rectangle((0.08, y - 0.055), 0.22, 0.11,
                               facecolor=col, edgecolor="black", linewidth=0.8))
        ax.text(0.34, y, name, fontsize=14, va="center")
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=600, bbox_inches="tight",
                    pad_inches=0.05, facecolor="white")
    plt.close(fig)


def process_patient(pnum, input_dir, output_dir):
    out = os.path.join(output_dir, f"Patient{pnum}_CD8_BaselineVsPostNivo")
    os.makedirs(out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(input_dir, f"Patient{pnum}_*.tiff")))
    stacks = {os.path.basename(f): tifffile.imread(f).astype(np.float32) for f in files}
    bounds = shared_bounds(list(stacks.values()))

    bl_key, pn_key = PAIRS[pnum]
    bl, pn = stacks[bl_key], stacks[pn_key]
    h = min(bl.shape[1], pn.shape[1])
    w = min(bl.shape[2], pn.shape[2])

    render_panel(center_crop(bl, h, w), bounds,
                 os.path.join(out, f"Patient{pnum}_Baseline_representative.png"))
    render_panel(center_crop(pn, h, w), bounds,
                 os.path.join(out, f"Patient{pnum}_PostNivo_representative.png"))
    save_legend(os.path.join(out, "legend"))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Directory with 34-channel IMC TIFFs")
    ap.add_argument("--output", required=True, help="Output directory")
    ap.add_argument("--patients", nargs="+", type=int, default=[1, 4, 6, 9])
    args = ap.parse_args()
    os.makedirs(args.output, exist_ok=True)
    for p in args.patients:
        process_patient(p, args.input, args.output)


if __name__ == "__main__":
    main()
