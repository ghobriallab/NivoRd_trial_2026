# GitHub Actions workflows

| Workflow | Trigger | What it does |
|----------|---------|--------------|
| `run-reproducibility-check.yml` | Weekly (Sat 08:00 UTC) + manual `workflow_dispatch` | Runs `tools/check_reproducibility.py` against the repo and appends the score to `reproducibility.md` at the repo root. Commits the updated file as `chore: reproducibility check <date> — <score>/100`. |

The reproducibility scorer is GCS-aware: it expects per-folder READMEs to mention a `gs://` mirror path (or `gsutil` / `gcloud storage`) so that anyone cloning the repo can copy the matching data directory from a bucket. The scorer does NOT mount or read any data, only verifies the documentation.

The workflow uses `actions/checkout@v4`, so it only sees tracked files; gitignored directories like `data/` and `figures/` never reach the runner.
