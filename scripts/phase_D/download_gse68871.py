"""Track D — download GSE68871 (Terragna 2016; bortezomib-containing VTD in
newly-diagnosed multiple myeloma; n~118; GPL570 U133 Plus 2.0) as the 2nd
held-out cohort for n=2 replication. Reuses local GPL570 annot (GSE16446).
Caveat: VTD is a bortezomib COMBINATION (thalidomide+dex), not mono — weaker
than GSE9782 PS341 mono, but bortezomib-containing and independent.
"""
from __future__ import annotations
import gzip, hashlib, json, shutil, urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("./external_data")
DEST = ROOT / "GSE68871"; DEST.mkdir(parents=True, exist_ok=True)
SERIES_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68871/matrix/GSE68871_series_matrix.txt.gz"
GPL570_SRC = ROOT / "GSE16446" / "platform_annot.annot.gz"
GPL570_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def download(url, dest):
    log = {"url": url, "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    try:
        with urllib.request.urlopen(url, timeout=300) as r, open(dest, "wb") as f:
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
    sm = DEST / "series_matrix.txt.gz"
    out = {"cohort": "GSE68871", "platform": "GPL570", "files": []}
    print("downloading GSE68871 series matrix ...", flush=True)
    out["files"].append({"label": "series_matrix", **download(SERIES_URL, sm)})
    print("  ->", out["files"][-1]["status"], out["files"][-1].get("size_bytes"), flush=True)
    annot = DEST / "platform_annot.annot.gz"
    if GPL570_SRC.exists():
        shutil.copy(GPL570_SRC, annot)
        out["files"].append({"label": "platform_annot", "status": "copied_local", "sha256": sha(annot)})
    else:
        out["files"].append({"label": "platform_annot", **download(GPL570_URL, annot)})
    (DEST / "DOWNLOAD_LOG.json").write_text(json.dumps(out, indent=2))
    # peek characteristics
    chars = []
    if sm.exists() and out["files"][0]["status"] == "ok":
        with gzip.open(sm, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("!Sample_characteristics") or line.startswith("!Sample_title"):
                    chars.append(line.rstrip("\n")[:300])
                if line.startswith("!series_matrix_table_begin"):
                    break
        (DEST / "characteristics_preview.txt").write_text("\n".join(chars))
        print("=== characteristics preview ===")
        for c in chars[:18]:
            print(c[:200])
    print("DONE.")


if __name__ == "__main__":
    main()
