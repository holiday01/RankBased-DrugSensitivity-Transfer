"""Stage 1: Download GEO Series Matrix + GPL annotation for external cohorts.

Cohorts (per notes.md D3):
  GSE16446 (single-arm epirubicin, ER-/TN BC, pCR, U133+2)
  GSE25066 (TFAC neoadjuvant ~508, pCR, U133A)
  GSE22226 (I-SPY1 AC-T ~150, pCR + RCB, U133A)

DNS rate-limit protocol (per memory feedback_network_rate_limit):
  - SERIAL only, no parallel downloads
  - 60-second sleep between cohorts
  - One download attempt per cohort; on failure log + continue
  - GEO uses NCBI FTP/HTTPS, no MCP alternative

Outputs:
  external_data/<gse>/series_matrix.txt.gz
  external_data/<gse>/<platform>_annot.csv (probe -> symbol map)
  external_data/<gse>/DOWNLOAD_LOG.json (timestamps, sha256, URLs)
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

OUT_ROOT = Path("./external_data")

# GEO platform annot URL convention: GPL{thousand_group}nnn/GPL{N}/annot/GPL{N}.annot.gz
#   GPL96 → GPLnnn (no thousands group since N < 1000)
#   GPL570 → GPL5nnn (thousands digit = 0; GEO uses leading-digit grouping)
#   GPL1708 → GPL1nnn
# GSE22226 (I-SPY1) is multi-platform; series matrix is split per platform:
#   GSE22226-GPL1708_series_matrix.txt.gz  (Agilent G4112A 22K, 130 patients)
COHORTS = [
    {
        "gse": "GSE16446",
        "platform": "GPL570",
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16446/matrix/GSE16446_series_matrix.txt.gz",
        # GPL570 < 1000 → GPLnnn group (same as GPL96)
        "platform_annot_url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz",
        "notes": "Single-arm epirubicin, ER-/TN BC, pCR; U133+2",
    },
    {
        "gse": "GSE25066",
        "platform": "GPL96",
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25066/matrix/GSE25066_series_matrix.txt.gz",
        "platform_annot_url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
        "notes": "TFAC neoadjuvant ~508, pCR; U133A",
    },
    {
        "gse": "GSE22226",
        "platform": "GPL1708",
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22226/matrix/GSE22226-GPL1708_series_matrix.txt.gz",
        "platform_annot_url": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1708/annot/GPL1708.annot.gz",
        "notes": "I-SPY1 AC-T ~150 prospective, pCR + RCB; Agilent G4112A (GPL1708 subseries)",
    },
]

SLEEP_BETWEEN = 60  # seconds; DNS rate-limit compliance


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    log = {
        "url": url,
        "dest": str(dest),
        "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        log["status"] = "ok"
        log["sha256"] = sha256_of(dest)
        log["size_bytes"] = dest.stat().st_size
    except Exception as e:
        log["status"] = "error"
        log["error"] = str(e)
    log["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return log


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    for i, c in enumerate(COHORTS):
        cohort_dir = OUT_ROOT / c["gse"]
        cohort_dir.mkdir(parents=True, exist_ok=True)
        # Resume mode: skip files that already exist with sha256
        prev_log_path = cohort_dir / "DOWNLOAD_LOG.json"
        prev_log = {}
        if prev_log_path.exists():
            try:
                prev_log = json.loads(prev_log_path.read_text())
            except Exception:
                prev_log = {}
        prev_ok = {f["label"] for f in prev_log.get("files", []) if f.get("status") == "ok"}

        download_log = {"cohort": c["gse"], "platform": c["platform"], "notes": c["notes"], "files": []}
        any_new = False
        for label, url in [
            ("series_matrix", c["series_matrix_url"]),
            ("platform_annot", c["platform_annot_url"]),
        ]:
            ext = ".txt.gz" if "matrix" in label else ".annot.gz"
            dest = cohort_dir / f"{label}{ext}"
            if label in prev_ok and dest.exists():
                # Carry over previous OK entry
                ok_entry = next((f for f in prev_log["files"] if f["label"] == label), None)
                if ok_entry:
                    download_log["files"].append(ok_entry)
                    print(f"[{c['gse']}] {label} already OK, skipping", flush=True)
                    continue
            any_new = True
            print(f"[{c['gse']}] downloading {label} ...", flush=True)
            log = download(url, dest)
            download_log["files"].append({"label": label, **log})
            print(f"  -> {log['status']}; {log.get('size_bytes', '?')} bytes", flush=True)
            time.sleep(5)

        with open(cohort_dir / "DOWNLOAD_LOG.json", "w") as f:
            json.dump(download_log, f, indent=2)

        if any_new and i < len(COHORTS) - 1:
            print(f"  sleeping {SLEEP_BETWEEN}s before next cohort (DNS rate-limit) ...", flush=True)
            time.sleep(SLEEP_BETWEEN)

    print("All cohorts processed.")


if __name__ == "__main__":
    main()
