# `data/` directory

Data files used by the figure scripts are hosted on Zenodo and are NOT
committed to this repository.

**Concept DOI:** https://doi.org/10.5281/zenodo.19430235
**Latest version:** v4-DOR (DOI assigned at publication; check the
concept-DOI page for the latest version).

## How to populate this directory

From the repository root, run:

```bash
python download_data.py
```

This will fetch every file in the deposit into `data/` (creating the
directory if it does not already exist) and verify file sizes against
the values in the script.

You can override the destination with the `NIVO_DATA_DIR` environment
variable; in that case, set it before importing any of the figure
scripts.

## What is on Zenodo (and what is not)

See the deposit's own README and DATA_DICTIONARY for the full list. Note:

- The CD138-negative immune-cell and CD138-positive tumor-cell
  scRNA-seq h5ads (`CD138neg_immune_cells.h5ad`,
  `CD138pos_tumor_cells.h5ad`) are **not** in the v4-DOR deposit pending
  scrubbing of dates from cell-barcode metadata. They are needed for
  Fig 2A–F and Supp S2A/S2B/S2E. Until they are added to a future
  Zenodo version, these scripts will not run end-to-end. Contact the
  authors for controlled-access copies.
- Master DGE / GSEA result CSVs from the Letter-to-Blood era are NOT
  in the v4-DOR deposit; they have been superseded by the per-figure
  outputs the scripts now generate locally.
