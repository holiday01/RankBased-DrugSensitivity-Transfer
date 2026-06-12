"""Track D replication — download GSE136337 (MM bulk cohort, ~426, survival) and
peek platform + data-row format (probe vs gene) + survival/treatment fields.
Generic: detects !Series_platform_id so we can fetch the right annot next.
"""
from __future__ import annotations
import gzip, hashlib, json, urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("./external_data")
DEST = ROOT / "GSE136337"; DEST.mkdir(parents=True, exist_ok=True)
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136337/matrix/GSE136337_series_matrix.txt.gz"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    sm = DEST / "series_matrix.txt.gz"
    log = {"cohort": "GSE136337", "url": URL, "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    try:
        with urllib.request.urlopen(URL, timeout=600) as r, open(sm, "wb") as f:
            while True:
                ch = r.read(1 << 20)
                if not ch:
                    break
                f.write(ch)
        log.update(status="ok", sha256=sha(sm), size_bytes=sm.stat().st_size)
    except Exception as e:
        log.update(status="error", error=str(e))
    (DEST / "DOWNLOAD_LOG.json").write_text(json.dumps(log, indent=2))
    print("download:", log["status"], log.get("size_bytes"))
    if log["status"] != "ok":
        return

    platform = None; chars = []; first_data = []; in_body = False; nrow = 0
    with gzip.open(sm, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                platform = line.split("\t", 1)[1].strip().strip('"')
            elif line.startswith("!Sample_characteristics") or line.startswith("!Sample_title"):
                chars.append(line.rstrip("\n")[:260])
            elif line.startswith("!series_matrix_table_begin"):
                in_body = True; continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_body:
                first_data.append(line.split("\t")[0].strip('"'))
                nrow += 1
                if nrow >= 6:
                    break
    print("PLATFORM:", platform)
    print("first row IDs (probe vs gene?):", first_data[:6])
    print("=== characteristics ===")
    for c in chars[:16]:
        print("  ", c[:200])
    (DEST / "peek.txt").write_text(f"platform={platform}\nfirst_ids={first_data}\n" + "\n".join(chars))
    print("DONE.")


if __name__ == "__main__":
    main()
