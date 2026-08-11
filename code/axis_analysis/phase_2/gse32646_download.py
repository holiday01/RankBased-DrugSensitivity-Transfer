"""Phase 2 / B1 — GSE32646 download + parse (Miyake et al. 2012, PMID 22320227).

Series Matrix already downloaded into external_data/GSE32646/series_matrix.txt.gz.
GPL570 annotation reused from GSE16446 directory (same platform).

Outputs (into resubmission_v2/results/phase_2/):
  - GSE32646_expression.tsv  (rows=patient_id, cols=HGNC_symbol; max-IQR probe)
  - GSE32646_phenotype.tsv   (per-sample metadata)
  - GSE32646_pcr.tsv         (patient_id, pcr 0/1)
  - GSE32646_DOWNLOAD_LOG.json
"""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
COHORT = "GSE32646"
COHORT_DIR = REPO / "external_data" / COHORT
PHASE_2 = REPO / "resubmission_v2" / "results" / "phase_2"
ALIAS_PATH = REPO / "external_data" / "hgnc_alias_table.tsv"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_series_matrix(path: Path):
    opener = gzip.open if path.name.endswith(".gz") else open
    metadata = {"sample_chars": {"characteristics": []}}
    body_lines = []
    in_body = False
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_body = True
                continue
            if line == "!series_matrix_table_end":
                in_body = False
                continue
            if in_body:
                body_lines.append(line)
            else:
                if line.startswith("!Sample_geo_accession"):
                    metadata["sample_ids"] = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_characteristics_ch1"):
                    metadata["sample_chars"]["characteristics"].append(
                        [x.strip('"') for x in line.split("\t")[1:]]
                    )
                elif line.startswith("!Sample_title"):
                    metadata["sample_titles"] = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Series_platform_id"):
                    metadata["platform_id"] = line.split("\t", 1)[1].strip('"')
    body = "\n".join(body_lines)
    expr = pd.read_csv(StringIO(body), sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    expr.index.name = "probe_id"
    return expr, metadata


def read_gpl_annot(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.name.endswith(".gz") else open
    rows = []
    header = None
    in_data = False
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                in_data = True
                continue
            if line.startswith("!platform_table_end"):
                in_data = False
                continue
            if in_data:
                if header is None:
                    header = line.split("\t")
                else:
                    rows.append(line.split("\t"))
    df = pd.DataFrame(rows, columns=header)
    sym_col = None
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("gene symbol", "gene_symbol", "genesymbol", "symbol"):
            sym_col = c
            break
    if sym_col is None:
        for c in df.columns:
            cl = c.lower().strip()
            if cl.startswith("gene ") and "symbol" in cl and "unigene" not in cl:
                sym_col = c
                break
    if sym_col is None:
        for c in df.columns:
            cl = c.lower()
            if "symbol" in cl and "unigene" not in cl and "platform" not in cl:
                sym_col = c
                break
    if sym_col is None:
        raise ValueError(f"Cannot find gene symbol column: {df.columns.tolist()}")
    df = df.rename(columns={sym_col: "gene_symbol"})
    if "ID" in df.columns:
        df = df.rename(columns={"ID": "probe_id"})
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.upper().str.strip()
    return df[["probe_id", "gene_symbol"]]


def load_hgnc_aliases() -> dict:
    if not ALIAS_PATH.exists():
        return {}
    df = pd.read_csv(ALIAS_PATH, sep="\t")
    return dict(zip(df.iloc[:, 0].astype(str).str.upper(), df.iloc[:, 1].astype(str).str.upper()))


def collapse_probes_max_iqr(expr: pd.DataFrame, probe_to_symbol: dict):
    expr_t = expr.T  # sample x probe
    valid_probes = [p for p in expr_t.columns if probe_to_symbol.get(p, "") not in {"", "NAN", "NONE"}]
    diag = {
        "n_probes_total": int(expr_t.shape[1]),
        "n_probes_with_symbol": len(valid_probes),
    }
    expr_t = expr_t[valid_probes]
    q75 = expr_t.quantile(0.75)
    q25 = expr_t.quantile(0.25)
    iqr = (q75 - q25).fillna(0)
    probe_df = pd.DataFrame(
        {"probe": valid_probes, "symbol": [probe_to_symbol[p] for p in valid_probes], "iqr": iqr.values}
    )
    idx = probe_df.groupby("symbol")["iqr"].idxmax()
    picked = probe_df.loc[idx]
    picked = picked[picked["symbol"].notna() & (picked["symbol"] != "")]
    sym_to_probe = dict(zip(picked["symbol"], picked["probe"]))
    collapsed = expr_t[list(sym_to_probe.values())].copy()
    collapsed.columns = list(sym_to_probe.keys())
    diag["n_unique_symbols_after_collapse"] = int(collapsed.shape[1])
    diag["n_samples"] = int(collapsed.shape[0])
    return collapsed, diag


def parse_phenotype(metadata: dict) -> pd.DataFrame:
    sample_ids = metadata.get("sample_ids", [])
    chars_lists = metadata.get("sample_chars", {}).get("characteristics", [])
    rows = []
    for i, sid in enumerate(sample_ids):
        row = {"patient_id": sid}
        for line in chars_lists:
            val = line[i] if i < len(line) else ""
            if ":" in val:
                k, v = val.split(":", 1)
                row[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def extract_pcr(phen: pd.DataFrame) -> pd.DataFrame:
    """GSE32646 stores pCR as `pathologic_response_pcr_ncr` ∈ {pCR, nCR}."""
    # Find pcr column (case-insensitive, contains 'pcr')
    pcr_col = None
    for c in phen.columns:
        if "pcr" in c.lower() and "ncr" in c.lower():
            pcr_col = c
            break
    if pcr_col is None:
        # fallback: any column with 'pcr'
        for c in phen.columns:
            if "pcr" in c.lower():
                pcr_col = c
                break
    if pcr_col is None:
        raise ValueError(f"pCR column not found in phenotype columns: {phen.columns.tolist()}")

    pcr_map = {"pcr": 1, "ncr": 0}
    pcr = phen[pcr_col].astype(str).str.strip().str.lower().map(pcr_map)
    out = pd.DataFrame({"patient_id": phen["patient_id"], "pcr": pcr})

    # Stratification covariates
    er_col = next((c for c in phen.columns if "er_status" in c.lower()), None)
    pr_col = next((c for c in phen.columns if "pr_status" in c.lower()), None)
    her2_col = next((c for c in phen.columns if "her2_status" in c.lower()), None)
    out["er_status"] = phen[er_col].astype(str).str.strip().str.lower() if er_col else "unknown"
    out["pr_status"] = phen[pr_col].astype(str).str.strip().str.lower() if pr_col else "unknown"
    out["her2_status"] = phen[her2_col].astype(str).str.strip().str.lower() if her2_col else "unknown"
    out["regimen"] = "P-FEC"
    return out


def main():
    matrix_path = COHORT_DIR / "series_matrix.txt.gz"
    annot_path = COHORT_DIR / "platform_annot.annot.gz"

    if not matrix_path.exists():
        print(f"ERROR: series matrix missing at {matrix_path}", file=sys.stderr)
        sys.exit(1)
    if not annot_path.exists():
        print(f"ERROR: GPL570 annotation missing at {annot_path}", file=sys.stderr)
        sys.exit(1)

    PHASE_2.mkdir(parents=True, exist_ok=True)

    log: dict = {
        "cohort": COHORT,
        "platform": "GPL570",
        "started_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "series_matrix_path": str(matrix_path),
        "series_matrix_sha256": sha256_of(matrix_path),
        "series_matrix_size_bytes": matrix_path.stat().st_size,
        "platform_annot_path": str(annot_path),
        "platform_annot_sha256": sha256_of(annot_path),
        "source_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32646/matrix/GSE32646_series_matrix.txt.gz",
        "reference": "Miyake et al. 2012, PMID 22320227 (Osaka University, Noguchi group); n=115 BC neoadjuvant P-FEC",
    }

    print(f"[GSE32646] Parsing series matrix ...", flush=True)
    expr, metadata = read_series_matrix(matrix_path)
    log["raw_probe_x_sample_shape"] = list(expr.shape)
    print(f"  raw shape (probes x samples): {expr.shape}")

    print(f"[GSE32646] Parsing GPL570 annotation ...", flush=True)
    annot = read_gpl_annot(annot_path)
    probe_to_symbol = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    aliases = load_hgnc_aliases()
    log["hgnc_aliases_loaded"] = len(aliases)
    if aliases:
        probe_to_symbol = {p: aliases.get(s, s) for p, s in probe_to_symbol.items()}

    collapsed, diag = collapse_probes_max_iqr(expr, probe_to_symbol)
    log["probe_collapse"] = diag
    print(f"  collapsed (sample x symbol): {collapsed.shape}")

    phen = parse_phenotype(metadata)
    log["phenotype_columns"] = list(phen.columns)
    print(f"  phenotype rows: {len(phen)}; cols: {list(phen.columns)}")

    pcr_df = extract_pcr(phen)
    n_total = len(pcr_df)
    pcr_df = pcr_df.dropna(subset=["pcr"]).copy()
    pcr_df["pcr"] = pcr_df["pcr"].astype(int)
    log["n_phenotype_total"] = n_total
    log["n_with_pcr_label"] = len(pcr_df)
    log["n_responders"] = int((pcr_df["pcr"] == 1).sum())
    log["n_non_responders"] = int((pcr_df["pcr"] == 0).sum())
    log["pcr_rate"] = float((pcr_df["pcr"] == 1).mean()) if len(pcr_df) else float("nan")
    log["er_status_counts"] = pcr_df["er_status"].value_counts().to_dict()
    log["her2_status_counts"] = pcr_df["her2_status"].value_counts().to_dict()

    # Detect ESR1 / PGR / ERBB2
    expr_cols = set(collapsed.columns)
    for gene in ["ESR1", "PGR", "ERBB2"]:
        log[f"{gene}_present"] = gene in expr_cols

    # Save outputs
    expr_out = PHASE_2 / "GSE32646_expression.tsv"
    phen_out = PHASE_2 / "GSE32646_phenotype.tsv"
    pcr_out = PHASE_2 / "GSE32646_pcr.tsv"
    log_out = PHASE_2 / "GSE32646_DOWNLOAD_LOG.json"

    collapsed.index.name = "patient_id"
    collapsed.to_csv(expr_out, sep="\t", float_format="%.10g")
    phen.to_csv(phen_out, sep="\t", index=False)
    pcr_df.to_csv(pcr_out, sep="\t", index=False)

    log["expression_output"] = str(expr_out)
    log["phenotype_output"] = str(phen_out)
    log["pcr_output"] = str(pcr_out)
    log["finished_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with open(log_out, "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n[GSE32646] Summary:")
    print(f"  n samples with pCR: {log['n_with_pcr_label']} (responders: {log['n_responders']}, non: {log['n_non_responders']}, rate: {log['pcr_rate']:.3f})")
    print(f"  ESR1: {log['ESR1_present']}, PGR: {log['PGR_present']}, ERBB2: {log['ERBB2_present']}")
    print(f"  Outputs: {expr_out.name}, {phen_out.name}, {pcr_out.name}")


if __name__ == "__main__":
    main()
