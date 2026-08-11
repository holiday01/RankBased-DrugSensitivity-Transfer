"""Track D replication — download GSE19784 (HOVON-65) and GSE2658 (UAMS TT),
both GPL570, serial with 60s gap (rate-limit). Reuse local GPL570 annot.
Peeks at characteristics to learn treatment-arm + survival field formats.
"""
from __future__ import annotations
import gzip, hashlib, json, shutil, time, urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/external_data")
GPL570_SRC = ROOT / "GSE16446" / "platform_annot.annot.gz"
COHORTS = [
    ("GSE19784", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19784/matrix/GSE19784_series_matrix.txt.gz"),
    ("GSE2658",  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE2nnn/GSE2658/matrix/GSE2658_series_matrix.txt.gz"),
]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def dl(url, dest):
    log = {"url": url, "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    try:
        with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as f:
            while True:
                ch = r.read(1 << 20)
                if not ch:
                    break
                f.write(ch)
        log.update(status="ok", sha256=sha(dest), size_bytes=dest.stat().st_size)
    except Exception as e:
        log.update(status="error", error=str(e))
    log["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    return log


def main():
    for i, (gse, url) in enumerate(COHORTS):
        d = ROOT / gse; d.mkdir(parents=True, exist_ok=True)
        sm = d / "series_matrix.txt.gz"
        out = {"cohort": gse, "platform": "GPL570", "files": []}
        print(f"[{gse}] downloading series matrix ...", flush=True)
        out["files"].append({"label": "series_matrix", **dl(url, sm)})
        print("  ->", out["files"][-1]["status"], out["files"][-1].get("size_bytes"), flush=True)
        annot = d / "platform_annot.annot.gz"
        if GPL570_SRC.exists():
            shutil.copy(GPL570_SRC, annot)
            out["files"].append({"label": "platform_annot", "status": "copied_local", "sha256": sha(annot)})
        (d / "DOWNLOAD_LOG.json").write_text(json.dumps(out, indent=2))
        # peek metadata
        if sm.exists() and out["files"][0]["status"] == "ok":
            chars = []
            with gzip.open(sm, "rt", errors="replace") as f:
                for line in f:
                    if line.startswith("!Sample_characteristics") or line.startswith("!Sample_title") \
                       or line.startswith("!Sample_source"):
                        chars.append(line.rstrip("\n")[:300])
                    if line.startswith("!series_matrix_table_begin"):
                        break
            (d / "characteristics_preview.txt").write_text("\n".join(chars))
            print(f"  [{gse}] characteristics lines: {len(chars)}")
            for c in chars[:20]:
                print("   ", c[:160])
        if i < len(COHORTS) - 1:
            print("  sleeping 60s (rate-limit) ...", flush=True)
            time.sleep(60)
    print("DONE.")


if __name__ == "__main__":
    main()
