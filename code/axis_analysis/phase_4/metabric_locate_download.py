"""Phase 4 / M1 — Locate or download METABRIC (Curtis 2012 / Pereira 2016).

Strategy:
 1. Search the local filesystem for an existing METABRIC copy first
    (/mnt/10t/...; /home/holiday01/...).
 2. If not found, download from the cBioPortal datahub mirror on GitHub
    (LFS media URLs are public, no auth needed).
 3. Parse:
      - data_mrna_illumina_microarray.txt        (raw log-intensity)
      - data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt (z-scores)
      - data_clinical_sample.txt
      - data_clinical_patient.txt
 4. Emit:
      - METABRIC_expression.tsv          (samples × HGNC symbols, raw log-intensity)
      - METABRIC_expression_zscore.tsv   (samples × HGNC symbols, z-scored)
      - METABRIC_clinical.tsv            (sample-level merged with patient-level)
      - METABRIC_DOWNLOAD_LOG.json       (URLs, sha256, n samples, etc.)

cBioPortal METABRIC schema (sample-level CHEMOTHERAPY column = YES/NO).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import shutil
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
PHASE_4 = REPO / "resubmission_v2" / "results" / "phase_4" / "metabric"

# cBioPortal datahub LFS media URLs (public, no auth)
DATAHUB_BASE = "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_metabric"
FILES = {
    "data_mrna_illumina_microarray.txt": f"{DATAHUB_BASE}/data_mrna_illumina_microarray.txt",
    "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt": f"{DATAHUB_BASE}/data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",
    "data_clinical_sample.txt": f"{DATAHUB_BASE}/data_clinical_sample.txt",
    "data_clinical_patient.txt": f"{DATAHUB_BASE}/data_clinical_patient.txt",
}

LOCAL_SEARCH_ROOTS = [
    Path("/mnt/10t"),
    Path("/home/holiday01"),
]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def search_local() -> dict:
    """Look for any pre-downloaded METABRIC copy.
    Returns mapping {filename: local_path} for any files found."""
    found: dict[str, Path] = {}
    needles = list(FILES.keys())
    for root in LOCAL_SEARCH_ROOTS:
        if not root.exists():
            continue
        # Quick scan of likely sub-trees
        for sub in [root / "brca_metabric", root / "METABRIC", root / "metabric"]:
            if sub.exists() and sub.is_dir():
                for name in needles:
                    p = sub / name
                    if p.exists() and name not in found:
                        found[name] = p
    return found


def http_download(url: str, dest: Path, timeout: int = 120) -> int:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (compatible; metabric-download/1.0)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        data = r.read()
    dest.write_bytes(data)
    return len(data)


def parse_expression(path: Path) -> pd.DataFrame:
    """Parse cBioPortal-format expression matrix.
    First col = Hugo_Symbol, second col = Entrez_Gene_Id, rest = samples.
    Returns DataFrame with samples as rows and HGNC symbols as columns.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "Hugo_Symbol" not in df.columns:
        raise RuntimeError(f"Expected 'Hugo_Symbol' column in {path}; got {list(df.columns)[:5]}")
    sample_cols = [c for c in df.columns if c not in ("Hugo_Symbol", "Entrez_Gene_Id")]
    # Collapse duplicate Hugo_Symbols by max-IQR
    df["Hugo_Symbol"] = df["Hugo_Symbol"].astype(str).str.upper().str.strip()
    df = df[df["Hugo_Symbol"].notna() & (df["Hugo_Symbol"] != "")]
    mat = df[["Hugo_Symbol"] + sample_cols].copy()
    if mat["Hugo_Symbol"].duplicated().any():
        # IQR-based selection of representative probe row per gene
        vals = mat[sample_cols]
        iqr = (vals.quantile(0.75, axis=1) - vals.quantile(0.25, axis=1)).fillna(-np.inf)
        mat = mat.assign(_iqr=iqr).sort_values("_iqr", ascending=False).drop_duplicates(
            "Hugo_Symbol", keep="first"
        ).drop(columns="_iqr")
    mat = mat.set_index("Hugo_Symbol")[sample_cols]
    # Samples × genes
    out = mat.T
    # Coerce to float; drop genes that are all NaN
    out = out.apply(pd.to_numeric, errors="coerce")
    out = out.dropna(axis=1, how="all")
    return out


def parse_clinical(sample_path: Path, patient_path: Path) -> pd.DataFrame:
    """cBioPortal clinical files have 4 metadata header lines starting with '#'."""

    def read_skip_meta(p: Path) -> pd.DataFrame:
        with open(p) as f:
            lines = f.readlines()
        skip = 0
        for line in lines:
            if line.startswith("#"):
                skip += 1
            else:
                break
        return pd.read_csv(p, sep="\t", skiprows=skip, low_memory=False)

    sample = read_skip_meta(sample_path)
    patient = read_skip_meta(patient_path)
    # Use SAMPLE_ID + PATIENT_ID as keys
    merged = sample.merge(patient, on="PATIENT_ID", how="left", suffixes=("", "_PAT"))
    return merged


def main():
    PHASE_4.mkdir(parents=True, exist_ok=True)
    log = {
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "out_dir": str(PHASE_4),
        "files": {},
    }

    # Step 1: search local
    local_hits = search_local()
    log["local_hits"] = {k: str(v) for k, v in local_hits.items()}

    # Step 2: download anything missing into PHASE_4/raw/
    raw_dir = PHASE_4 / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        dest = raw_dir / name
        if name in local_hits:
            shutil.copyfile(local_hits[name], dest)
            log["files"][name] = {
                "source": "local",
                "local_path": str(local_hits[name]),
                "dest": str(dest),
                "size_bytes": dest.stat().st_size,
                "sha256": sha256_of(dest),
            }
            print(f"[M1] {name}: copied from local {local_hits[name]}")
        elif dest.exists() and dest.stat().st_size > 1024:
            log["files"][name] = {
                "source": "preexisting",
                "dest": str(dest),
                "size_bytes": dest.stat().st_size,
                "sha256": sha256_of(dest),
            }
            print(f"[M1] {name}: pre-existing in raw/, skipping download")
        else:
            print(f"[M1] {name}: downloading from {url}", flush=True)
            n = http_download(url, dest)
            log["files"][name] = {
                "source": "download",
                "url": url,
                "dest": str(dest),
                "size_bytes": n,
                "sha256": sha256_of(dest),
            }
            print(f"        -> {n} bytes")

    # Step 3: parse
    print("[M1] parsing raw expression (log-intensity)...", flush=True)
    expr_log = parse_expression(raw_dir / "data_mrna_illumina_microarray.txt")
    expr_log_out = PHASE_4 / "METABRIC_expression.tsv"
    expr_log.to_csv(expr_log_out, sep="\t", float_format="%.6g")
    log["expression_log_shape"] = list(expr_log.shape)
    log["expression_log_path"] = str(expr_log_out)
    log["expression_log_sha256"] = sha256_of(expr_log_out)

    print("[M1] parsing z-score expression...", flush=True)
    expr_z = parse_expression(raw_dir / "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
    expr_z_out = PHASE_4 / "METABRIC_expression_zscore.tsv"
    expr_z.to_csv(expr_z_out, sep="\t", float_format="%.6g")
    log["expression_zscore_shape"] = list(expr_z.shape)
    log["expression_zscore_path"] = str(expr_z_out)
    log["expression_zscore_sha256"] = sha256_of(expr_z_out)

    print("[M1] parsing clinical...", flush=True)
    clin = parse_clinical(raw_dir / "data_clinical_sample.txt", raw_dir / "data_clinical_patient.txt")
    clin_out = PHASE_4 / "METABRIC_clinical.tsv"
    clin.to_csv(clin_out, sep="\t", index=False)
    log["clinical_shape"] = list(clin.shape)
    log["clinical_columns"] = list(clin.columns)
    log["clinical_path"] = str(clin_out)
    log["clinical_sha256"] = sha256_of(clin_out)

    # Quick QC counts
    if "CHEMOTHERAPY" in clin.columns:
        log["chemotherapy_value_counts"] = clin["CHEMOTHERAPY"].value_counts(dropna=False).to_dict()
    log["n_samples_in_expression_log"] = expr_log.shape[0]
    log["n_samples_in_clinical"] = int(clin["SAMPLE_ID"].nunique()) if "SAMPLE_ID" in clin.columns else None

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    out_log = PHASE_4 / "METABRIC_DOWNLOAD_LOG.json"
    with open(out_log, "w") as f:
        json.dump(log, f, indent=2, default=str)
    print(f"[M1] log written: {out_log}")
    print(f"[M1] expression(log) shape: {expr_log.shape}; expression(z) shape: {expr_z.shape}; clinical shape: {clin.shape}")


if __name__ == "__main__":
    main()
