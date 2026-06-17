#!/usr/bin/env python3
"""
check_reproducibility.py

Purpose:    Evaluate the reproducibility of a computational biology project directory.
            Checks architecture, naming, documentation, and GCP data hygiene.
            Does NOT check code correctness — focus is on workflow structure.

Inputs:     PATH  — directory to evaluate (default: current working directory)

Outputs:    Human-readable report to stdout
            Optional JSON score file (--json flag)

Usage:
    pixi run python tools/check_reproducibility.py [PATH] [--json] [--min-score N]
    pixi run check-repro                                 # shorthand (evaluates current dir)

Dependencies: Python 3.8+ (stdlib only)
"""

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import List, Tuple

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT_EXTENSIONS = {".R", ".py", ".sh", ".bash", ".Rmd", ".ipynb"}
HEADER_KEYWORDS = {"purpose", "inputs", "outputs", "dependencies"}
DATA_EXTENSIONS = {".bam", ".cram", ".vcf", ".vcf.gz", ".fastq", ".fastq.gz",
                   ".fq", ".fq.gz", ".csv", ".tsv", ".rds", ".rda", ".h5",
                   ".h5ad", ".loom", ".bed", ".gtf", ".gff"}
LARGE_FILE_THRESHOLD_MB = 10
LOCAL_PATH_PATTERNS = re.compile(
    r'(?<!["\'])(/Users/|/home/|/mnt/|/tmp/|/var/|/opt/|/data/|/scratch/)'
)
# Patterns that are false-positives in the evaluation tool itself
SELF_EXEMPT_PATTERNS = re.compile(r"LOCAL_PATH_PATTERNS|check_reproducibility")
GCS_MENTION_PATTERNS = re.compile(
    r'gs://|gsutil|gcloud storage', re.IGNORECASE
)
PHI_PATTERNS = [
    (re.compile(r'\b\d{3}-\d{2}-\d{4}\b'), "SSN-like pattern"),
    (re.compile(r'\b\d{3}-\d{3}-\d{4}\b'), "Phone-like pattern"),
    (re.compile(r'(?i)\b(mrn|medical.record.number|date.of.birth|dob)\s*[=:]\s*\S'), "PHI field"),
    (re.compile(r'(?i)(password|passwd|secret|api_key|auth_token)\s*=\s*["\']?\S'), "Credential"),
    (re.compile(r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b'), "Email address"),
]
IGNORED_DIRS = {".git", "__pycache__", ".ipynb_checkpoints", "node_modules", ".venv", "venv", ".pixi"}
NUMERIC_PREFIX_RE = re.compile(r'^\d{2,3}[_\-]')
BAD_NAME_CHARS_RE = re.compile(r'[ ()\[\]&;!]')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CategoryResult:
    name: str
    score: float
    max_score: int
    issues: List[str] = field(default_factory=list)

    @property
    def int_score(self) -> int:
        return round(self.score)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def is_text_file(path: Path) -> bool:
    try:
        with open(path, "rb") as f:
            chunk = f.read(1024)
        return b"\x00" not in chunk
    except OSError:
        return False


def iter_scripts(root: Path):
    for p in root.rglob("*"):
        if any(part in IGNORED_DIRS for part in p.parts):
            continue
        if p.is_file() and p.suffix in SCRIPT_EXTENSIONS:
            yield p


def iter_dirs(root: Path):
    for p in root.rglob("*"):
        if any(part in IGNORED_DIRS for part in p.parts):
            continue
        if p.is_dir() and not p.name.startswith("."):
            yield p
    yield root


def read_head(path: Path, n_lines: int = 40) -> str:
    try:
        with open(path, errors="replace") as f:
            return "\n".join(f.readline() for _ in range(n_lines))
    except OSError:
        return ""


def bar(score: float, max_score: int, width: int = 10) -> str:
    filled = round((score / max_score) * width) if max_score else 0
    return "█" * filled + "░" * (width - filled)


def truncate_list(items: List[str], n: int = 3) -> str:
    if len(items) <= n:
        return ", ".join(items)
    return ", ".join(items[:n]) + f" (+{len(items) - n} more)"


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_step_ordering(root: Path) -> CategoryResult:
    """Scripts/folders should use numeric prefixes (01_, 02_) to show analysis order."""
    max_score = 20
    issues = []
    scripts = list(iter_scripts(root))

    if not scripts:
        return CategoryResult("Step Ordering", max_score, max_score,
                              ["No scripts found — nothing to evaluate"])

    missing_prefix = [s for s in scripts if not NUMERIC_PREFIX_RE.match(s.name)]
    ratio = len(missing_prefix) / len(scripts)

    # Deduct proportionally; full score if ≤20% lack prefix
    penalty = max(0.0, (ratio - 0.2) / 0.8) * max_score
    score = max_score - penalty

    if missing_prefix:
        examples = truncate_list([s.name for s in missing_prefix])
        issues.append(
            f"{len(missing_prefix)}/{len(scripts)} scripts lack a numeric prefix "
            f"(e.g. 01_, 02_): {examples}"
        )

    # Also check top-level analysis subdirs
    top_dirs = [d for d in root.iterdir()
                if d.is_dir() and not d.name.startswith(".") and d.name not in IGNORED_DIRS]
    unordered_dirs = [d for d in top_dirs if not NUMERIC_PREFIX_RE.match(d.name)]
    if len(top_dirs) > 2 and len(unordered_dirs) == len(top_dirs):
        issues.append(
            "Top-level analysis folders have no numeric ordering — "
            "consider renaming to 01_qc/, 02_align/, etc."
        )
        score = max(0.0, score - 5)

    return CategoryResult("Step Ordering", min(score, max_score), max_score, issues)


def check_documentation(root: Path) -> CategoryResult:
    """Each script needs a header block; each folder needs a README."""
    max_score = 25
    issues = []

    scripts = list(iter_scripts(root))
    missing_header = []
    for s in scripts:
        if s.suffix == ".ipynb":
            continue  # notebooks have their own structure
        head = read_head(s).lower()
        found = {kw for kw in HEADER_KEYWORDS if kw in head}
        if len(found) < 3:
            missing_header.append(s.relative_to(root))

    dirs = list(iter_dirs(root))
    missing_readme = []
    for d in dirs:
        if not (d / "README.md").exists() and not (d / "readme.md").exists():
            missing_readme.append(d.relative_to(root))

    header_ratio = 1.0 - (len(missing_header) / len(scripts)) if scripts else 1.0
    readme_ratio = 1.0 - (len(missing_readme) / len(dirs)) if dirs else 1.0

    score = (header_ratio * 15) + (readme_ratio * 10)

    if missing_header:
        examples = truncate_list([str(p) for p in missing_header])
        issues.append(
            f"{len(missing_header)} script(s) missing header block "
            f"(need Purpose, Inputs, Outputs, Dependencies): {examples}"
        )
    if missing_readme:
        examples = truncate_list([str(p) if str(p) != "." else "(root)" for p in missing_readme])
        issues.append(f"{len(missing_readme)} folder(s) missing README.md: {examples}")

    return CategoryResult("Documentation", score, max_score, issues)


def check_path_hygiene(root: Path) -> CategoryResult:
    """No hardcoded local absolute paths; use relative paths or gs:// URIs."""
    max_score = 20
    issues = []
    offending = []

    for s in iter_scripts(root):
        if not is_text_file(s):
            continue
        try:
            content = s.read_text(errors="replace")
        except OSError:
            continue
        if SELF_EXEMPT_PATTERNS.search(content):
            continue  # skip files that define/test path patterns (e.g. this script)
        matches = LOCAL_PATH_PATTERNS.findall(content)
        if matches:
            offending.append((s.relative_to(root), len(matches)))

    if offending:
        ratio = len(offending) / max(1, len(list(iter_scripts(root))))
        penalty = ratio * max_score
        score = max_score - penalty
        examples = truncate_list([f"{p} ({n} hit{'s' if n>1 else ''})" for p, n in offending])
        issues.append(
            f"{len(offending)} script(s) contain hardcoded local paths "
            f"(/Users/, /home/, /mnt/, etc.): {examples}"
        )
    else:
        score = max_score

    return CategoryResult("Path Hygiene", score, max_score, issues)


def check_gcs_handling(root: Path) -> CategoryResult:
    """READMEs mention GCS copy instructions; .gitignore covers data files; no large files committed."""
    max_score = 15
    issues = []
    score = max_score

    # 1. .gitignore presence and coverage
    gitignore = root / ".gitignore"
    if not gitignore.exists():
        issues.append(".gitignore is missing — add one covering data extensions (.bam, .vcf, .rds, .h5, .csv, etc.)")
        score -= 5
    else:
        content = gitignore.read_text(errors="replace").lower()
        missing_exts = [ext for ext in [".bam", ".vcf", ".rds", ".h5", ".csv", ".fastq"]
                        if ext not in content and ext.lstrip(".") not in content]
        if missing_exts:
            issues.append(
                f".gitignore doesn't cover these data extensions: {', '.join(missing_exts)}"
            )
            score -= 2

    # 2. READMEs mention GCS copy/mount instructions
    readmes = list(root.rglob("README.md")) + list(root.rglob("readme.md"))
    readmes = [r for r in readmes if not any(part in IGNORED_DIRS for part in r.parts)]
    readmes_without_gcs = []
    for readme in readmes:
        try:
            content = readme.read_text(errors="replace")
        except OSError:
            continue
        if not GCS_MENTION_PATTERNS.search(content):
            readmes_without_gcs.append(readme.relative_to(root))

    if readmes_without_gcs:
        examples = truncate_list([str(p) for p in readmes_without_gcs])
        issues.append(
            f"{len(readmes_without_gcs)} README(s) don't mention GCS (gs://, gsutil, or gcloud storage) "
            f"— add mount/copy instructions: {examples}"
        )
        deduction = min(5, len(readmes_without_gcs) * 2)
        score -= deduction
    elif not readmes:
        issues.append("No README.md found — cannot check GCS copy instructions")
        score -= 3

    # 3. Large data files present in the tree
    large_files = []
    for p in root.rglob("*"):
        if any(part in IGNORED_DIRS for part in p.parts):
            continue
        if p.is_file() and p.suffix.lower() in {ext.lower() for ext in DATA_EXTENSIONS}:
            size_mb = p.stat().st_size / (1024 * 1024)
            if size_mb > LARGE_FILE_THRESHOLD_MB:
                large_files.append((p.relative_to(root), size_mb))

    if large_files:
        examples = truncate_list([f"{p} ({s:.0f} MB)" for p, s in large_files])
        issues.append(
            f"{len(large_files)} large data file(s) found in repo (>{LARGE_FILE_THRESHOLD_MB} MB) — "
            f"move to GCS bucket: {examples}"
        )
        score -= min(3, len(large_files))

    return CategoryResult("GCS Data Handling", max(0.0, score), max_score, issues)


def check_naming(root: Path) -> CategoryResult:
    """File/folder names: no spaces or special chars; scripts should be snake_case."""
    max_score = 10
    issues = []
    bad_names = []

    for p in root.rglob("*"):
        if any(part in IGNORED_DIRS for part in p.parts):
            continue
        if p.name.startswith("."):
            continue
        if BAD_NAME_CHARS_RE.search(p.name):
            bad_names.append(p.relative_to(root))

    if bad_names:
        ratio = len(bad_names) / max(1, sum(1 for _ in root.rglob("*")))
        penalty = min(max_score, max(1, round(ratio * max_score * 5)))
        score = max(0.0, max_score - penalty)
        examples = truncate_list([str(p) for p in bad_names])
        issues.append(
            f"{len(bad_names)} file/folder name(s) contain spaces or special characters: {examples}"
        )
    else:
        score = max_score

    return CategoryResult("Naming Conventions", score, max_score, issues)


def check_phi_safety(root: Path) -> CategoryResult:
    """Scan for PHI-like patterns and credentials in code."""
    max_score = 10
    issues = []
    hits: List[Tuple[str, str, Path]] = []

    for s in root.rglob("*"):
        if any(part in IGNORED_DIRS for part in s.parts):
            continue
        if not s.is_file():
            continue
        if not is_text_file(s):
            continue
        if s.suffix.lower() in {".md", ".txt", ".rst"}:
            continue  # skip docs for PHI scan (they may have examples)
        try:
            content = s.read_text(errors="replace")
        except OSError:
            continue
        for pattern, label in PHI_PATTERNS:
            if pattern.search(content):
                hits.append((label, str(s.relative_to(root)), s))
                break  # one hit per file is enough

    if hits:
        penalty = min(max_score, len(hits) * 3)
        score = max(0.0, max_score - penalty)
        examples = truncate_list([f"{label} in {path}" for label, path, _ in hits])
        issues.append(
            f"Potential PHI or credentials found in {len(hits)} file(s): {examples}"
        )
    else:
        score = max_score

    return CategoryResult("PHI / Credential Safety", score, max_score, issues)


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------

def render_report(root: Path, results: List[CategoryResult], total: int) -> str:
    width = 52
    lines = []
    lines.append("=" * width)
    lines.append(f" Reproducibility Report")
    lines.append(f" {root}")
    lines.append("=" * width)

    for r in results:
        b = bar(r.score, r.max_score)
        status = "✓" if not r.issues else f"— {r.issues[0][:55]}" if r.issues else ""
        lines.append(
            f"{r.name:<20} [{b}] {r.int_score:>3}/{r.max_score:<3} {status}"
        )

    lines.append("-" * width)
    lines.append(f"TOTAL SCORE: {total}/100")
    lines.append("")

    # Collect all issues, ranked by score deficit
    all_issues = []
    for r in results:
        for issue in r.issues:
            all_issues.append((r.name, issue))

    if all_issues:
        top = all_issues[:5]
        lines.append("Top issues to fix (run again after addressing these):")
        for i, (cat, issue) in enumerate(top, 1):
            lines.append(f"  {i}. [{cat}] {issue}")
        if len(all_issues) > 5:
            lines.append(f"  ... and {len(all_issues) - 5} more issues (fix the above first)")
        lines.append("")

    if total >= 80:
        lines.append("Good reproducibility! Address remaining issues when convenient.")
    elif total >= 60:
        lines.append(f"Score {total}/100 — fix the issues above and run again.")
    else:
        lines.append(
            f"Critical reproducibility issues (score {total}/100). "
            "Fix the top issues before sharing this analysis."
        )

    lines.append("=" * width)
    return "\n".join(lines)


def build_json(root: Path, results: List[CategoryResult], total: int) -> dict:
    return {
        "score": total,
        "timestamp": date.today().isoformat(),
        "path": str(root),
        "categories": {
            r.name.lower().replace(" ", "_").replace("/", "_"): {
                "score": r.int_score,
                "max": r.max_score,
                "issues": r.issues,
            }
            for r in results
        },
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Evaluate reproducibility of a compbio project directory."
    )
    parser.add_argument(
        "path",
        nargs="?",
        default=".",
        help="Directory to evaluate (default: current directory)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Write reproducibility_score.json in the target directory",
    )
    parser.add_argument(
        "--min-score",
        type=int,
        default=60,
        metavar="N",
        help="Exit with code 1 if score is below N (default: 60)",
    )
    args = parser.parse_args()

    root = Path(args.path).resolve()
    if not root.is_dir():
        print(f"Error: {root} is not a directory", file=sys.stderr)
        sys.exit(2)

    results = [
        check_step_ordering(root),
        check_documentation(root),
        check_path_hygiene(root),
        check_gcs_handling(root),
        check_naming(root),
        check_phi_safety(root),
    ]

    total = sum(r.int_score for r in results)

    print(render_report(root, results, total))

    if args.json:
        out_path = root / "reproducibility_score.json"
        out_path.write_text(json.dumps(build_json(root, results, total), indent=2))
        print(f"JSON written to: {out_path}")

    if total < args.min_score:
        sys.exit(1)


if __name__ == "__main__":
    main()
