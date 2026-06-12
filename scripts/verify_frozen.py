"""Verify integrity of frozen_signatures/ against MANIFEST.json.

Modes:
  --strict        Verify every entry: TSVs + script + source files + auxiliary artifacts.
                  Requires source files (CCLE/OncoKB) present at recorded paths.
  --archive-only  Verify only what is bundled in the repo (TSVs + scripts + planning .md).
                  Skips external source files. Use this on a cloned/downloaded archive.

  --repo-root PATH  Override the repo root (default: parent of this script).

Exits 0 if all match, non-zero on mismatch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_entry(file_path: Path, expected_sha: str, label: str, errors: list):
    if not file_path.exists():
        errors.append(f"[{label}] missing: {file_path}")
        return
    actual = sha256_of(file_path)
    if actual != expected_sha:
        errors.append(f"[{label}] sha256 mismatch: {file_path}\n   expected={expected_sha}\n   actual  ={actual}")


def verify_pool(pool_manifest: dict, repo_root: Path, label: str, errors: list):
    """entry['file'] is a path relative to repo_root; do NOT prepend additional prefix."""
    for strategy_name, entries in pool_manifest.get("outputs", {}).items():
        for entry in entries:
            file_path = repo_root / entry["file"]
            verify_entry(file_path, entry["sha256"], f"{label}/{strategy_name}/{entry['drug_norm']}", errors)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=None,
                    help="Repo root (default: parent of scripts/)")
    ap.add_argument("--archive-only", action="store_true",
                    help="Skip external source-file checks (CCLE/OncoKB); verify only bundled artifacts.")
    ap.add_argument("--strict", action="store_true",
                    help="Verify everything including external source files (default if neither flag).")
    args = ap.parse_args()

    repo_root = args.repo_root if args.repo_root else Path(__file__).resolve().parent.parent.parent
    manifest_path = repo_root / "." / "frozen_signatures" / "MANIFEST.json"
    if not manifest_path.exists():
        print(f"FAIL: MANIFEST.json not found at {manifest_path}", file=sys.stderr)
        sys.exit(2)

    with open(manifest_path) as f:
        manifest = json.load(f)

    errors = []

    # Freeze script self-hash
    script_rel = Path(manifest["software"]["script_path"]).name
    # script_path is absolute; rebuild from repo_root
    script_path = repo_root / "." / "scripts" / script_rel
    verify_entry(script_path, manifest["software"]["script_sha256"], "freeze_script", errors)

    # Auxiliary artifacts (scripts + planning .md)
    for aux in manifest.get("auxiliary_artifacts", []):
        verify_entry(repo_root / aux["path"], aux["sha256"], "aux/" + Path(aux["path"]).name, errors)

    # External source files — only in strict mode
    if not args.archive_only:
        for name, meta in manifest.get("source_files", {}).items():
            path = Path(meta["path"])
            if not path.exists():
                errors.append(f"source file missing (skip with --archive-only): {name} @ {path}")
                continue
            actual = sha256_of(path)
            if actual != meta["sha256"]:
                errors.append(f"source file {name} sha256 mismatch:\n   expected={meta['sha256']}\n   actual  ={actual}")

    # Per-pool signature TSVs (entry["file"] is already repo-root-relative)
    verify_pool(manifest.get("pan_cancer", {}), repo_root, "pan_cancer", errors)
    verify_pool(manifest.get("solid_only", {}), repo_root, "solid_only", errors)

    if errors:
        print(f"\n{len(errors)} verification error(s):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    n_pan = sum(len(v) for v in manifest.get("pan_cancer", {}).get("outputs", {}).values())
    n_solid = sum(len(v) for v in manifest.get("solid_only", {}).get("outputs", {}).values())
    n_aux = len(manifest.get("auxiliary_artifacts", []))
    mode = "archive-only" if args.archive_only else "strict"
    print(f"OK ({mode}): {n_pan} pan_cancer + {n_solid} solid_only + {n_aux} auxiliary; all sha256 match.")


if __name__ == "__main__":
    main()
