# tools/

Repository-housekeeping utilities (not part of the manuscript pipeline).

## `check_reproducibility.py`

Reproducibility scorer for this repo. Runs the same six checks the
GhobrialLab bot applies across lab repos:

| Category | Max | What it checks |
| --- | --- | --- |
| Step Ordering | 20 | Scripts use numeric prefixes; top-level analysis folders are ordered |
| Documentation | 25 | Each script has a Purpose / Inputs / Outputs / Dependencies header; each folder has a `README.md` |
| Path Hygiene | 20 | No hardcoded local paths (`/Users/`, `/home/`, `/mnt/`, …) |
| GCS Data Handling | 15 | `.gitignore` covers data extensions; READMEs mention `gs://` or `gsutil`; no large data files committed |
| Naming Conventions | 10 | No spaces or special chars in file/folder names |
| PHI / Credential Safety | 10 | No SSN/phone/MRN/email/credential patterns in code |

The scoring logic does NOT execute any pipeline code; it is a static-analysis
walk over file names, headers, `.gitignore`, and string patterns.

### Local run

```bash
python tools/check_reproducibility.py .          # report only
python tools/check_reproducibility.py . --json   # also write reproducibility_score.json
```

### Automated run

`.github/workflows/run-reproducibility-check.yml` runs the same scorer every
Saturday at 08:00 UTC (and on `workflow_dispatch`). When it completes it
prepends a new entry to `reproducibility.md` at the repo root and pushes the
update back to `main`.

The check uses only the Python stdlib (3.8+); no extra dependencies are
required. The scorer is identical to the one used by the IMPACT_2026 repo —
sourced from there and kept in sync.

Data for the manuscript pipeline is downloaded via `00_download_data.py` from
Zenodo or, for fast CI access, mirrored to a GCS bucket
(`gs://ghobrial-nivord-data/`).
