"""Track D — download GSE9782 (bortezomib in multiple myeloma; Mulligan 2007).
Platform GPL96 (U133A). Series matrix is split per platform; we take GPL96.
Reuses the existing GPL96 .annot.gz already downloaded for GSE25066.

Rate-limit: single cohort, serial. NCBI FTP/HTTPS.
"""
from __future__ import annotations
import hashlib, json, shutil, urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/external_data")
DEST = ROOT / "GSE9782"
DEST.mkdir(parents=True, exist_ok=True)

SERIES_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9782/matrix/GSE9782-GPL96_series_matrix.txt.gz"
GPL96_ANNOT_SRC = ROOT / "GSE25066" / "platform_annot.annot.gz"  # GPL96, already local
GPL96_ANNOT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def download(url, dest):
    log = {"url": url, "dest": str(dest), "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
        log["status"] = "ok"; log["sha256"] = sha(dest); log["size_bytes"] = dest.stat().st_size
    except Exception as e:
        log["status"] = "error"; log["error"] = str(e)
    log["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return log


def main():
    out_log = {"cohort": "GSE9782", "platform": "GPL96", "files": []}
    sm = DEST / "series_matrix.txt.gz"
    print("downloading GSE9782-GPL96 series matrix ...", flush=True)
    out_log["files"].append({"label": "series_matrix", **download(SERIES_URL, sm)})
    print("  ->", out_log["files"][-1]["status"], out_log["files"][-1].get("size_bytes"), flush=True)

    # GPL96 annot: reuse local copy if present, else download
    annot_dest = DEST / "platform_annot.annot.gz"
    if GPL96_ANNOT_SRC.exists():
        shutil.copy(GPL96_ANNOT_SRC, annot_dest)
        out_log["files"].append({"label": "platform_annot", "status": "copied_local",
                                 "src": str(GPL96_ANNOT_SRC), "sha256": sha(annot_dest)})
        print("  GPL96 annot copied from local GSE25066", flush=True)
    else:
        out_log["files"].append({"label": "platform_annot", **download(GPL96_ANNOT_URL, annot_dest)})

    (DEST / "DOWNLOAD_LOG.json").write_text(json.dumps(out_log, indent=2))

    # quick peek at characteristics lines to learn response-label format
    import gzip
    chars = []
    with gzip.open(sm, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_characteristics") or line.startswith("!Sample_title"):
                chars.append(line.rstrip("\n")[:400])
            if line.startswith("!series_matrix_table_begin"):
                break
    (DEST / "characteristics_preview.txt").write_text("\n".join(chars))
    print("\n=== characteristics preview (first 400 chars/line) ===")
    for c in chars[:15]:
        print(c[:240])
    print("DONE.")


if __name__ == "__main__":
    main()
