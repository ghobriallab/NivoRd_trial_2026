#!/usr/bin/env python3
"""Download Zenodo data files for the NivoRd v4-DOR manuscript.

The Zenodo record version this script targets is 19988651. The public,
filename-based URL pattern is:

    https://zenodo.org/records/{record_id}/files/{filename}

These resolve only after the deposit is published (until then, the draft
files live behind authenticated bucket URLs). This script will give a
clear error if the URLs cannot be reached and tell you to either wait
for publication or fetch the files manually from
https://zenodo.org/uploads/{record_id} (draft) or
https://zenodo.org/records/{record_id} (published).

Usage:
    python 00_download_data.py
"""
from __future__ import annotations
import os
import sys
import urllib.request
from pathlib import Path

# Update RECORD_ID once the new Zenodo version is published. The
# concept-DOI URL https://doi.org/10.5281/zenodo.19430235 always
# resolves to the latest version.
RECORD_ID = 19988651

# (filename, expected size in bytes)
# README.md and DATA_DICTIONARY.md sizes are checked loosely (>0) since
# they are routinely re-uploaded with description tweaks.
FILES = [
    ("README.md", None),
    ("DATA_DICTIONARY.md", None),
    ("LICENSE", 835),
    ("IMC-import.h5ad", 16836901),
    # Nivo-only filtered (external NBM / ELO reference cohorts removed pre-deposit).
    ("CD138neg_immune_cells.h5ad", 151844414),
    ("CD138pos_tumor_cells.h5ad", 53967587),
    ("phenotyped_canonical.csv", 7152685),
    ("pheno_8clus-sample_pct.csv", 10711),
    ("pheno_8clus-cn_centers.csv", 753),
    ("imc_extras_persample.csv", 7508),
    ("IMC_antibody_panel.csv", 310),
    ("differential_cd4_macrophage_pretreatment.csv", 32061),
    ("differential_plasma_macrophage_pretreatment.csv", 36557),
    ("gene_expression_summary.csv", 4311),
    ("cytogenetic_validation.csv", 587),
    ("MSigDB_Hallmark_2020.gmt", 44235),
    ("patient_response_groups.csv", 148),
]

DEST = Path(os.environ.get("NIVO_DATA_DIR", str(Path(__file__).resolve().parent / "data")))
URL_TEMPLATE = "https://zenodo.org/records/{record_id}/files/{filename}"


def fetch(filename: str, expected_size: int | None) -> None:
    url = URL_TEMPLATE.format(record_id=RECORD_ID, filename=filename)
    target = DEST / filename
    if target.exists() and (
        expected_size is None and target.stat().st_size > 0
        or expected_size is not None and target.stat().st_size == expected_size
    ):
        size_text = f"{target.stat().st_size} bytes"
        print(f"[ok ] {filename} (already present, {size_text})")
        return
    print(f"[get] {filename} <- {url}")
    try:
        with urllib.request.urlopen(url) as resp, open(target, "wb") as out:
            while True:
                chunk = resp.read(1 << 16)
                if not chunk:
                    break
                out.write(chunk)
    except Exception as exc:  # noqa: BLE001
        print(f"[err] {filename}: {exc}")
        print("      If the deposit is still a draft, fetch from "
              f"https://zenodo.org/uploads/{RECORD_ID} via the web UI.")
        sys.exit(1)
    actual = target.stat().st_size
    if expected_size is not None and actual != expected_size:
        print(f"[warn] {filename}: size mismatch (got {actual}, expected {expected_size})")
    else:
        print(f"[ok ] {filename} ({actual} bytes)")


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {len(FILES)} files from Zenodo record {RECORD_ID} into {DEST}")
    for filename, size in FILES:
        fetch(filename, size)
    print("done.")


if __name__ == "__main__":
    main()
