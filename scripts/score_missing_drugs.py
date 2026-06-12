"""Phase 1 prelude — Score missing drug signatures (PACLITAXEL, CARBOPLATIN,
CYCLOPHOSPHAMIDE, FLUOROURACIL, GEFITINIB) on the 3 external BC cohorts that
do not already exist.

Outputs are written to external_data/<cohort>/score_<DRUG>.tsv
(consistent with existing DOXORUBICIN naming convention).

Delegates per-drug-per-cohort scoring to scripts/score_external_cohort.py.
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(".")
SCRIPT = REPO / "scripts" / "score_external_cohort.py"
SIG_DIR = REPO / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all"

JOBS = [
    ("GSE16446", "PACLITAXEL"),
    ("GSE25066", "PACLITAXEL"),
    ("GSE16446", "CARBOPLATIN"),
    ("GSE25066", "CARBOPLATIN"),
    ("GSE22226", "CARBOPLATIN"),
    ("GSE16446", "CYCLOPHOSPHAMIDE"),
    ("GSE25066", "CYCLOPHOSPHAMIDE"),
    ("GSE22226", "CYCLOPHOSPHAMIDE"),
    ("GSE16446", "FLUOROURACIL"),
    ("GSE25066", "FLUOROURACIL"),
    ("GSE22226", "FLUOROURACIL"),
    ("GSE25066", "GEFITINIB"),
    ("GSE22226", "GEFITINIB"),
]


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "delegate": str(SCRIPT),
        "delegate_sha256": sha256_of(SCRIPT),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "runs": [],
    }
    for cohort, drug in JOBS:
        sig = SIG_DIR / f"{drug}.tsv"
        out_tsv = REPO / "external_data" / cohort / f"score_{drug}.tsv"
        expr = REPO / "external_data" / cohort / f"{cohort}_expression.tsv"
        labels = REPO / "external_data" / cohort / f"{cohort}_pcr.tsv"
        if out_tsv.exists():
            print(f"[{drug}/{cohort}] SKIP (exists)")
            log["runs"].append({"cohort": cohort, "drug": drug, "status": "skip_exists"})
            continue
        if not sig.exists():
            print(f"[{drug}/{cohort}] FAIL: signature missing {sig}", file=sys.stderr)
            log["runs"].append({"cohort": cohort, "drug": drug, "status": "fail_no_sig"})
            continue
        cmd = [
            sys.executable, str(SCRIPT),
            "--signature", str(sig),
            "--expression", str(expr),
            "--labels", str(labels),
            "--output", str(out_tsv),
            "--skip-perm1",
        ]
        try:
            res = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=900)
            print(f"[{drug}/{cohort}] OK -> {out_tsv.name} | {res.stdout.strip()}")
            log["runs"].append({"cohort": cohort, "drug": drug, "status": "ok",
                                 "stdout": res.stdout.strip()})
        except subprocess.CalledProcessError as e:
            print(f"[{drug}/{cohort}] FAIL: {e.stderr.strip()[:400]}", file=sys.stderr)
            log["runs"].append({"cohort": cohort, "drug": drug, "status": "fail",
                                 "stderr": e.stderr.strip()[-1500:]})
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    out_log = REPO / "." / "results" / "phase_1" / "score_missing_drugs.log"
    out_log.parent.mkdir(parents=True, exist_ok=True)
    with open(out_log, "w") as f:
        json.dump(log, f, indent=2)
    print(f"Log: {out_log}")


if __name__ == "__main__":
    main()
