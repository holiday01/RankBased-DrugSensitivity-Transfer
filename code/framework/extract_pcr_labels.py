"""Extract standardized pCR labels + stratification covariates per cohort.

Reads cohort-specific phenotype.tsv (parse_gse_to_matrix.py output) and writes:
  external_data/<cohort>/<cohort>_pcr.tsv   columns:
    patient_id, pcr (0/1), pam50_basal (0/1, where derivable),
    er_status, her2_status, regimen, rcb_class (where available)

Per cohort label semantics:
  GSE16446: column `pcr` ∈ {YES, NO, NA}
  GSE25066: column `pathologic_response_pcr_rd` ∈ {pCR, RD}
  GSE22226: column `pathological_complete_response_(pcr)` ∈ {Yes, No}
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EXTERNAL_ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/external_data")


def find_col(df, *candidates):
    """Case-insensitive contains lookup."""
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        for cl, orig in cols_lower.items():
            if cand.lower() in cl:
                return orig
    return None


def extract_gse16446(phen: pd.DataFrame) -> pd.DataFrame:
    pcr_col = find_col(phen, "pcr")
    out = pd.DataFrame({"patient_id": phen["patient_id"]})
    # GSE16446 stores pcr as numeric (0.0 / 1.0 / NaN); coerce to int via numeric
    out["pcr"] = pd.to_numeric(phen[pcr_col], errors="coerce")
    # ER-/HER2- cohort (single stratum); PAM50 not in metadata — compute via PAM50 later
    out["pam50_basal"] = np.nan  # will be filled by PAM50 computation in scoring step
    out["er_status"] = "negative"  # uniform per cohort (ER- pre-selected)
    out["her2_status"] = "negative"  # uniform per cohort (TN pre-selected)
    out["regimen"] = "epirubicin_monotherapy"
    out["rcb_class"] = np.nan
    return out


def extract_gse25066(phen: pd.DataFrame) -> pd.DataFrame:
    pcr_col = find_col(phen, "pathologic_response_pcr_rd")
    er_col = find_col(phen, "er_status_ihc")
    her2_col = find_col(phen, "her2_status")
    rcb_col = find_col(phen, "pathologic_response_rcb_class")
    out = pd.DataFrame({"patient_id": phen["patient_id"]})
    pcr_map = {"pcr": 1, "rd": 0}
    out["pcr"] = phen[pcr_col].astype(str).str.strip().str.lower().map(pcr_map)
    out["pam50_basal"] = np.nan  # not in this cohort
    out["er_status"] = phen[er_col].astype(str).str.strip().str.upper().map({"P": "positive", "N": "negative", "I": "indeterminate"}).fillna("unknown") if er_col else "unknown"
    out["her2_status"] = phen[her2_col].astype(str).str.strip().str.upper().map({"P": "positive", "N": "negative", "I": "indeterminate"}).fillna("unknown") if her2_col else "unknown"
    out["regimen"] = "TFAC"
    out["rcb_class"] = phen[rcb_col] if rcb_col else np.nan
    return out


def extract_gse22226(phen: pd.DataFrame) -> pd.DataFrame:
    pcr_col = find_col(phen, "pathological_complete_response")
    pam50_col = find_col(phen, "intrinsic_subtype_by_pam50")
    er_col = find_col(phen, "er_(0=negative")  # "er_(0=negative;_1=positive;)"
    her2_col = find_col(phen, "her2_(0=negative")
    rcb_col = find_col(phen, "rcb_class")
    regimen_col = find_col(phen, "neoadjuvant_chemotherapy")

    out = pd.DataFrame({"patient_id": phen["patient_id"]})
    pcr_map = {"yes": 1, "no": 0}
    out["pcr"] = phen[pcr_col].astype(str).str.strip().str.lower().map(pcr_map)
    # PAM50 from metadata
    out["pam50_subtype"] = phen[pam50_col].astype(str).str.strip() if pam50_col else "unknown"
    out["pam50_basal"] = (out["pam50_subtype"].str.lower().str.contains("basal", na=False)).astype(int)
    er_num = pd.to_numeric(phen[er_col], errors="coerce") if er_col else None
    her2_num = pd.to_numeric(phen[her2_col], errors="coerce") if her2_col else None
    out["er_status"] = er_num.map({1.0: "positive", 0.0: "negative"}).fillna("unknown") if er_num is not None else "unknown"
    out["her2_status"] = her2_num.map({1.0: "positive", 0.0: "negative"}).fillna("unknown") if her2_num is not None else "unknown"
    out["regimen"] = phen[regimen_col].astype(str).str.strip() if regimen_col else "unknown"
    out["rcb_class"] = phen[rcb_col] if rcb_col else np.nan
    return out


EXTRACTORS = {
    "GSE16446": extract_gse16446,
    "GSE25066": extract_gse25066,
    "GSE22226": extract_gse22226,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=list(EXTRACTORS.keys()))
    args = ap.parse_args()

    phen_path = EXTERNAL_ROOT / args.cohort / f"{args.cohort}_phenotype.tsv"
    out_path = EXTERNAL_ROOT / args.cohort / f"{args.cohort}_pcr.tsv"
    log_path = EXTERNAL_ROOT / args.cohort / f"{args.cohort}_PCR_EXTRACT_LOG.json"

    phen = pd.read_csv(phen_path, sep="\t")
    extractor = EXTRACTORS[args.cohort]
    out = extractor(phen)
    # Drop patients with missing pCR
    n_total = len(out)
    out = out.dropna(subset=["pcr"]).copy()
    out["pcr"] = out["pcr"].astype(int)
    n_kept = len(out)
    out.to_csv(out_path, sep="\t", index=False)

    log = {
        "cohort": args.cohort,
        "n_total_in_phenotype": n_total,
        "n_with_pcr_label": n_kept,
        "n_responders": int((out["pcr"] == 1).sum()),
        "n_non_responders": int((out["pcr"] == 0).sum()),
        "pcr_rate": float((out["pcr"] == 1).mean()),
        "columns_in_output": list(out.columns),
    }
    if "pam50_subtype" in out.columns:
        log["pam50_subtype_counts"] = out["pam50_subtype"].value_counts().to_dict()
    if "er_status" in out.columns:
        log["er_status_counts"] = out["er_status"].value_counts().to_dict()
    if "her2_status" in out.columns:
        log["her2_status_counts"] = out["her2_status"].value_counts().to_dict()
    if "regimen" in out.columns:
        log["regimen_counts"] = out["regimen"].value_counts().to_dict()

    with open(log_path, "w") as f:
        json.dump(log, f, indent=2)
    print(json.dumps(log, indent=2))


if __name__ == "__main__":
    main()
