"""Self-lock the protocol document.

Reads:
  protocol.md
  frozen_signatures/MANIFEST.json
Writes:
  protocol.LOCK   (one-line file with sha256 + UTC + manifest sha256)

The LOCK file is the third-party-witnessable anchor: after this script writes
it, both files (protocol.md and protocol.LOCK) should be
uploaded together to Zenodo / OSF as a single deposit. The deposit's DOI +
upload timestamp become the immutable witness.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

# Portable: resolve relative to this script's location (scripts/)
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent  # 
PRE_REG = ROOT / "protocol.md"
LOCK = ROOT / "protocol.LOCK"
MANIFEST = ROOT / "frozen_signatures" / "MANIFEST.json"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    h.update(p.read_bytes())
    return h.hexdigest()


def main():
    if not PRE_REG.exists():
        raise SystemExit(f"Missing {PRE_REG}")
    if not MANIFEST.exists():
        raise SystemExit(f"Missing {MANIFEST}")

    prereg_sha = sha256_of(PRE_REG)
    manifest_sha = sha256_of(MANIFEST)
    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")

    payload = {
        "protocol_md_sha256": prereg_sha,
        "manifest_json_sha256": manifest_sha,
        "lock_timestamp_utc": ts,
        "instructions": (
            "Upload BOTH protocol.md and protocol.LOCK to Zenodo or OSF as a "
            "single deposit. The deposit's DOI + upload timestamp serve as the immutable "
            "third-party witness. To verify integrity later: re-hash protocol.md "
            "and MANIFEST.json, compare against this LOCK file."
        ),
    }
    LOCK.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Locked protocol.md")
    print(f"  protocol.md sha256: {prereg_sha}")
    print(f"  MANIFEST.json sha256:      {manifest_sha}")
    print(f"  Timestamp (UTC):           {ts}")
    print(f"  Written to:                {LOCK}")


if __name__ == "__main__":
    main()
