"""Phase 0 — Score 5 additional drug signatures (ETOPOSIDE, TRAMETINIB,
LAPATINIB, FULVESTRANT, IMATINIB) on the 3 external BC cohorts using
the existing singscore deployment (Methods Eq. 1).

This script delegates to scripts/score_external_cohort.py
to avoid re-implementing the rank-based singscore and bootstrap pipeline.

Output: score_<DRUG>_<cohort>.tsv per drug × cohort = 15 files
(plus the per-run JSON sidecar emitted by score_external_cohort.py).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(".")
SCRIPT_TO_REUSE = REPO / "scripts" / "score_external_cohort.py"
SIG_DIR = REPO / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all"
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]
DRUGS = ["ETOPOSIDE", "TRAMETINIB", "LAPATINIB", "FULVESTRANT", "IMATINIB"]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(REPO / "results/phase_0"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    script_path = Path(__file__).resolve()
    script_sha = sha256_of(script_path)

    log_path = out_dir / "score_extra_drugs.log"
    log_entries: dict = {
        "script": str(script_path),
        "script_sha256": script_sha,
        "delegate_script": str(SCRIPT_TO_REUSE),
        "delegate_script_sha256": sha256_of(SCRIPT_TO_REUSE),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "runs": [],
    }

    for cohort in COHORTS:
        expr = REPO / "external_data" / cohort / f"{cohort}_expression.tsv"
        labels = REPO / "external_data" / cohort / f"{cohort}_pcr.tsv"
        for drug in DRUGS:
            sig = SIG_DIR / f"{drug}.tsv"
            out_tsv = out_dir / f"score_{drug}_{cohort}.tsv"
            run = {
                "cohort": cohort,
                "drug": drug,
                "signature": str(sig),
                "output": str(out_tsv),
            }
            if not sig.exists():
                run["status"] = "fail"
                run["error"] = f"signature missing: {sig}"
                print(f"[{drug}/{cohort}] FAIL: signature missing", file=sys.stderr)
                log_entries["runs"].append(run)
                continue

            # For GSE16446 (all ER−/TN basal): single-stratum — skip within-stratum
            # label permutation per pre-reg §6.2 fallback. For GSE25066/GSE22226
            # we still skip --stratum-col because PAM50 col is in pam50.tsv not pcr.tsv;
            # GATE diagnostic in Phase 1+ does the stratified analysis using axes_*.tsv.
            cmd = [
                sys.executable, str(SCRIPT_TO_REUSE),
                "--signature", str(sig),
                "--expression", str(expr),
                "--labels", str(labels),
                "--output", str(out_tsv),
                "--skip-perm1",
            ]
            t0 = _dt.datetime.utcnow().isoformat() + "Z"
            try:
                res = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)
                run["status"] = "ok"
                run["stdout"] = res.stdout.strip()
                run["started_utc"] = t0
                run["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
                # capture json sidecar summary
                json_side = out_tsv.with_suffix(".json")
                if json_side.exists():
                    with open(json_side) as f:
                        meta = json.load(f)
                    run["auroc"] = meta.get("auroc")
                    run["n_patients_scored"] = meta.get("n_patients_scored")
                    run["n_up_genes_present_in_cohort"] = meta.get("n_up_genes_present_in_cohort")
                    run["n_down_genes_present_in_cohort"] = meta.get("n_down_genes_present_in_cohort")
                print(f"[{drug}/{cohort}] OK -> {out_tsv.name}  ({res.stdout.strip()})")
            except subprocess.CalledProcessError as e:
                run["status"] = "fail"
                run["error"] = e.stderr.strip()[-1500:]
                print(f"[{drug}/{cohort}] FAIL: {e.stderr.strip()[:400]}", file=sys.stderr)
            except Exception as e:
                run["status"] = "fail"
                run["error"] = repr(e)
                print(f"[{drug}/{cohort}] FAIL: {e}", file=sys.stderr)
            log_entries["runs"].append(run)

    log_entries["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(log_path, "w") as f:
        json.dump(log_entries, f, indent=2)
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
