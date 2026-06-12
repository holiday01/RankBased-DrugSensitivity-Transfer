"""Phase 3 / C1 — Locate (and stage) the TCGA-BRCA expression + clinical inputs.

Strategy:
  1. Look for an existing TCGA-BRCA expression matrix on disk
     (>= 500 samples). If found, symlink into phase_3/.
  2. Otherwise fall back to the UCSC Xena dump at
     external_data/tcga_brca/{HiSeqV2.gz, BRCA_clinicalMatrix}
     (downloaded once by the Phase 3 runner; SHA256 logged).
  3. Decompress + transpose expression matrix into the Phase 3 working format:
     rows = patient_id (TCGA aliquot barcode truncated to participant-level),
     cols = HGNC gene symbol. Saved as TSV.
  4. Merge Xena clinical with `TCGA/2016_tcga_drug.csv` (which carries the
     per-drug treatment records the Xena clinical matrix lacks) on
     participant ID. The merged table is what Phase 3 / C2 will filter on.

Outputs (under phase_3/):
  - tcga_brca_expression.tsv         : staged expression (samples x genes, log2 RSEM+1)
  - tcga_brca_clinical_merged.tsv    : Xena clinical + 2016_tcga_drug treatment rows
  - locate_tcga_brca.log             : provenance JSON
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(".")
EXT = REPO / "." / "external_data" / "tcga_brca"
PHASE3 = REPO / "." / "." / "results" / "phase_3"
DRUG_TABLE = REPO / "TCGA" / "2016_tcga_drug.csv"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def participant_id(barcode: str) -> str:
    """TCGA-AR-A5QQ-01 -> TCGA-AR-A5QQ (participant level, drops sample/aliquot)."""
    parts = str(barcode).split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else str(barcode)


def main() -> int:
    PHASE3.mkdir(parents=True, exist_ok=True)

    expr_gz = EXT / "HiSeqV2.gz"
    clin_xena = EXT / "BRCA_clinicalMatrix"
    if not expr_gz.exists() or not clin_xena.exists():
        print(f"FATAL: Xena dump missing in {EXT}", file=sys.stderr)
        return 2
    if not DRUG_TABLE.exists():
        print(f"FATAL: {DRUG_TABLE} not found", file=sys.stderr)
        return 2

    # 1) Read expression (genes x samples) and transpose to (samples x genes).
    print(f"[expr] reading {expr_gz} ...", flush=True)
    with gzip.open(expr_gz, "rt") as fh:
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    # rows = gene symbols, cols = TCGA aliquot barcodes
    expr = expr.T  # samples x genes
    expr.index.name = "sample_aliquot"

    n_samples_raw, n_genes_raw = expr.shape
    # Keep tumor primary (01A/01B etc) — TCGA "-01" sample code
    is_primary = expr.index.to_series().str.contains(r"-0[12](?:[A-Z]?)?$", regex=True)
    expr = expr[is_primary]
    n_samples_primary = expr.shape[0]

    # Drop duplicate samples (multiple aliquots per participant -> keep first)
    expr["_pid"] = expr.index.to_series().map(participant_id)
    expr = expr.drop_duplicates(subset=["_pid"], keep="first")
    expr.index = expr["_pid"].values
    expr.index.name = "patient_id"
    expr = expr.drop(columns=["_pid"])
    n_samples_dedup = expr.shape[0]

    # 2) Clinical merge.
    print(f"[clin] reading {clin_xena}", flush=True)
    cl_xena = pd.read_csv(clin_xena, sep="\t", low_memory=False)
    cl_xena["patient_id"] = cl_xena["sampleID"].map(participant_id)
    # Keep one clinical row per participant (Xena rows are sample-level)
    cl_xena = cl_xena.drop_duplicates(subset=["patient_id"], keep="first")

    print(f"[drug] reading {DRUG_TABLE}", flush=True)
    drugs = pd.read_csv(DRUG_TABLE, low_memory=False)
    drugs = drugs[drugs["Cancer"].str.contains("BRCA", na=False, regex=False)].copy()
    drugs["patient_id"] = drugs["bcr_patient_barcode"].astype(str)
    # collapse per-patient drug list
    drug_agg = (drugs.groupby("patient_id")
                .agg(drug_name_list=("drug_name", lambda s: "|".join(sorted(set(str(x) for x in s if pd.notna(x))))),
                     n_drug_records=("drug_name", "size"),
                     response_list=("measure_of_response", lambda s: "|".join(sorted(set(str(x) for x in s if pd.notna(x)))))
                     )
                .reset_index())

    merged = cl_xena.merge(drug_agg, on="patient_id", how="left")
    n_with_drugs = int(merged["drug_name_list"].notna().sum())

    # 3) Persist outputs (staged) — symlink original gz + write TSVs.
    out_expr = PHASE3 / "tcga_brca_expression.tsv"
    out_clin = PHASE3 / "tcga_brca_clinical_merged.tsv"
    out_log = PHASE3 / "locate_tcga_brca.log.json"

    expr.to_csv(out_expr, sep="\t", float_format="%.4f")
    merged.to_csv(out_clin, sep="\t", index=False)

    log = {
        "stage": "phase_3_C1_locate_tcga_brca",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "expression_gz": str(expr_gz),
            "expression_gz_sha256": sha256_of(expr_gz),
            "clinical_xena": str(clin_xena),
            "clinical_xena_sha256": sha256_of(clin_xena),
            "drug_table": str(DRUG_TABLE),
            "drug_table_sha256": sha256_of(DRUG_TABLE),
        },
        "expression_qc": {
            "n_samples_raw": int(n_samples_raw),
            "n_samples_primary_tumor": int(n_samples_primary),
            "n_samples_after_dedup_participant": int(n_samples_dedup),
            "n_genes": int(expr.shape[1]),
        },
        "clinical_qc": {
            "n_participants_xena": int(cl_xena.shape[0]),
            "n_brca_drug_participants_2016_table": int(drug_agg.shape[0]),
            "n_merged_rows": int(merged.shape[0]),
            "n_merged_with_drug_record": n_with_drugs,
        },
        "outputs": {
            "tcga_brca_expression.tsv": str(out_expr),
            "tcga_brca_clinical_merged.tsv": str(out_clin),
        },
        "notes": [
            "Expression source: UCSC Xena TCGA.BRCA HiSeqV2 (log2 RSEM+1, gene-level).",
            "Clinical source: Xena BRCA_clinicalMatrix; drug source: TCGA/2016_tcga_drug.csv.",
            "Sample-level rows collapsed to participant-level (first primary tumor).",
        ],
    }
    with open(out_log, "w") as f:
        json.dump(log, f, indent=2)
    print(json.dumps(log["expression_qc"] | log["clinical_qc"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
