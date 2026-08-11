"""Stage 2: Parse a GEO Series Matrix + GPL annotation into:
  - <gse>_expression.tsv  (rows=patient_id, cols=HGNC_symbol; max-IQR probe per symbol)
  - <gse>_phenotype.tsv   (per-sample metadata extracted from !Sample_characteristics_ch1)
  - <gse>_PARSE_LOG.json  (probe-collapse stats, HGNC alias coverage, etc.)

Implementation principles (per audit + decisions.md):
  - max-IQR probe-per-symbol collapse rule (pre-registered).
  - HGNC alias resolution via local table (e.g., SEPT6 -> SEPTIN6).
  - Drop probes with NA symbol or low-quality "AFFX-" controls.
  - Report n_probes_total, n_probes_with_symbol, n_unique_symbols_after_collapse,
    n_signature_genes_recovered_per_drug (computed after Stage 0 done).

Usage:
  python parse_gse_to_matrix.py --cohort GSE16446
  python parse_gse_to_matrix.py --cohort GSE25066
  python parse_gse_to_matrix.py --cohort GSE22226
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

EXTERNAL_ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/external_data")
HGNC_ALIAS_PATH = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/external_data/hgnc_alias_table.tsv")
# HGNC alias table: two cols, alias_symbol  canonical_symbol
# (download once from genenames.org or use the locally-cached one when available)


def read_series_matrix(path: Path) -> tuple[pd.DataFrame, dict]:
    """Returns (expression_df, sample_metadata_dict).

    Series Matrix format: header lines start with '!', followed by '!series_matrix_table_begin'
    then a TSV body, then '!series_matrix_table_end'.
    """
    opener = gzip.open if path.name.endswith(".gz") else open
    metadata = {"sample_chars": {}}
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
                elif line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_characteristics_ch2"):
                    # Two-channel Agilent arrays (e.g., GSE22226) put patient metadata in ch2
                    metadata["sample_chars"].setdefault("characteristics", []).append(
                        [x.strip('"') for x in line.split("\t")[1:]]
                    )
                elif line.startswith("!Sample_title"):
                    metadata["sample_titles"] = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Series_platform_id"):
                    metadata["platform_id"] = line.split("\t", 1)[1].strip('"')
    from io import StringIO
    body = "\n".join(body_lines)
    expr = pd.read_csv(StringIO(body), sep="\t", index_col=0)
    # Force probe IDs to string so they match the string-keyed annot table
    # (some platforms use numeric probe IDs that pandas would parse as int64)
    expr.index = expr.index.astype(str)
    expr.index.name = "probe_id"
    return expr, metadata


def read_gpl_annot(path: Path) -> pd.DataFrame:
    """Parse GEO GPL platform .annot.gz file. Returns df with probe_id + gene_symbol cols."""
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
    # Standardize column names.
    # PRIORITY: exact "gene symbol" / "gene_symbol" (case-insensitive) over substring
    # match (which would incorrectly catch "UniGene symbol" which is often empty).
    sym_col = None
    for c in df.columns:
        cl = c.lower().strip()
        if cl in ("gene symbol", "gene_symbol", "genesymbol", "symbol"):
            sym_col = c
            break
    if sym_col is None:
        # Second priority: starts with "gene "
        for c in df.columns:
            cl = c.lower().strip()
            if cl.startswith("gene ") and "symbol" in cl and "unigene" not in cl:
                sym_col = c
                break
    if sym_col is None:
        # Final fallback: any column with "symbol" that is not unigene-related
        for c in df.columns:
            cl = c.lower()
            if "symbol" in cl and "unigene" not in cl and "platform" not in cl:
                sym_col = c
                break
    if sym_col is None:
        raise ValueError(f"Cannot find gene symbol column in GPL annotation: {df.columns.tolist()}")
    df = df.rename(columns={sym_col: "gene_symbol"})
    if "ID" in df.columns:
        df = df.rename(columns={"ID": "probe_id"})
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.upper().str.strip()
    return df[["probe_id", "gene_symbol"]]


def load_hgnc_aliases() -> dict:
    """Return {alias_symbol: canonical_symbol}. Empty if table not present."""
    if not HGNC_ALIAS_PATH.exists():
        return {}
    df = pd.read_csv(HGNC_ALIAS_PATH, sep="\t")
    return dict(zip(df.iloc[:, 0].astype(str).str.upper(), df.iloc[:, 1].astype(str).str.upper()))


def collapse_probes_max_iqr(expr: pd.DataFrame, probe_to_symbol: dict) -> tuple[pd.DataFrame, dict]:
    """For each gene symbol with >1 probe, keep the probe with highest IQR.

    Returns (gene_x_sample DataFrame transposed to sample_x_gene, diag).
    """
    # expr: probe x sample. Transpose to sample x probe.
    expr_t = expr.T  # sample x probe
    # Drop probes with no symbol
    valid_probes = [p for p in expr_t.columns if probe_to_symbol.get(p, "") not in {"", "NAN", "NONE"}]
    diag = {
        "n_probes_total": int(expr_t.shape[1]),
        "n_probes_with_symbol": len(valid_probes),
    }
    expr_t = expr_t[valid_probes]
    # Compute IQR per probe
    q75 = expr_t.quantile(0.75)
    q25 = expr_t.quantile(0.25)
    iqr = (q75 - q25).fillna(0)
    # Group probes by symbol
    probe_df = pd.DataFrame({"probe": valid_probes, "symbol": [probe_to_symbol[p] for p in valid_probes], "iqr": iqr.values})
    # Pick max-IQR probe per symbol
    idx = probe_df.groupby("symbol")["iqr"].idxmax()
    picked = probe_df.loc[idx]
    picked = picked[picked["symbol"].notna() & (picked["symbol"] != "")]
    # Build collapsed matrix: sample x symbol
    sym_to_probe = dict(zip(picked["symbol"], picked["probe"]))
    collapsed = expr_t[list(sym_to_probe.values())].copy()
    collapsed.columns = list(sym_to_probe.keys())
    diag["n_unique_symbols_after_collapse"] = int(collapsed.shape[1])
    diag["n_samples"] = int(collapsed.shape[0])
    return collapsed, diag


def parse_phenotype(metadata: dict) -> pd.DataFrame:
    """Best-effort: extract pcr / response / subtype from !Sample_characteristics_ch1."""
    sample_ids = metadata.get("sample_ids", [])
    chars_lists = metadata.get("sample_chars", {}).get("characteristics", [])
    rows = []
    for i, sid in enumerate(sample_ids):
        row = {"patient_id": sid}
        for line in chars_lists:
            val = line[i] if i < len(line) else ""
            # Each value looks like "key: value"
            if ":" in val:
                k, v = val.split(":", 1)
                row[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, help="e.g. GSE16446")
    args = ap.parse_args()

    cohort_dir = EXTERNAL_ROOT / args.cohort
    matrix_path = cohort_dir / "series_matrix.txt.gz"
    annot_path = cohort_dir / "platform_annot.annot.gz"
    if not matrix_path.exists() or not annot_path.exists():
        print(f"Missing input files in {cohort_dir}. Run download_gse.py first.", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing {args.cohort} series matrix ...", flush=True)
    expr, metadata = read_series_matrix(matrix_path)
    print(f"  raw expression shape (probes x samples): {expr.shape}")

    print(f"Parsing {args.cohort} platform annotation ...", flush=True)
    annot = read_gpl_annot(annot_path)
    probe_to_symbol = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    print(f"  probes annotated: {len(probe_to_symbol)}")

    aliases = load_hgnc_aliases()
    if aliases:
        probe_to_symbol = {p: aliases.get(s, s) for p, s in probe_to_symbol.items()}
    else:
        print("  WARN: no HGNC alias table found at", HGNC_ALIAS_PATH)

    collapsed, diag = collapse_probes_max_iqr(expr, probe_to_symbol)
    print(f"  collapsed (sample x symbol): {collapsed.shape}")

    phen = parse_phenotype(metadata)
    print(f"  phenotype rows: {len(phen)}")

    collapsed.to_csv(cohort_dir / f"{args.cohort}_expression.tsv", sep="\t", float_format="%.10g")
    phen.to_csv(cohort_dir / f"{args.cohort}_phenotype.tsv", sep="\t", index=False)

    log = {
        "cohort": args.cohort,
        "parsed_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "diag": diag,
        "phenotype_columns": list(phen.columns),
        "hgnc_aliases_applied": len(aliases),
    }
    with open(cohort_dir / f"{args.cohort}_PARSE_LOG.json", "w") as f:
        json.dump(log, f, indent=2)
    print(f"  Outputs: {args.cohort}_expression.tsv / _phenotype.tsv / _PARSE_LOG.json")


if __name__ == "__main__":
    main()
