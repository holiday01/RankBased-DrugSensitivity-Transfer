"""Phase 3 / C2 — Build the TCGA-BRCA cohort tables for the §4.9 Mammaprint
worked example.

Two cohorts are written:

  (A) Full TCGA-BRCA cohort (1101 participants) — the PRIMARY worked example
      cohort because (1) Mammaprint is a *prognostic* signature developed and
      validated on a broad BC cohort, not a chemo-response predictor, and (2)
      the TCGA-BRCA chemo subset is severely under-powered for OS prognosis
      (median follow-up ~2.5y, only 10 deaths among 108 chemo-treated patients).
      DEVIATION from the original §4.9 plan: documented below + in
      PHASE_3_LOG.md.

  (B) Chemo-treated subset (sensitivity; defined as: ≥1 record of a canonical
      chemo agent in TCGA/2016_tcga_drug.csv OR
      history_of_neoadjuvant_treatment == 'YES').

Outcome:
  5-year and 3-year binarized OS-or-recurrence (label=1 = poor outcome).
  Per the data inspection in this folder, the 3-year window retains
  meaningfully more samples (n_labeled ≈ 446 full / 46 chemo) than 5-year
  (338 / 27) while keeping clinical relevance. We compute BOTH and primarily
  report 3-year.

  Construction:
    poor (label=1) =
        (vital_status=='DECEASED' AND days_to_death <= window)
      OR (new_tumor_event_after_initial_treatment=='YES' AND
          days_to_new_tumor_event_after_initial_treatment <= window)
    good (label=0) =
        (vital_status=='LIVING' AND days_to_last_followup >= window AND
         NOT recurrence-within-window)
    indeterminate (label=-1) = dropped from GATE.

PAM50 stratifier: Xena `PAM50Call_RNAseq`. Samples missing PAM50 are kept in
the marginal AUROC but dropped from within-stratum AUROC.

Outputs (phase_3/):
  - tcga_brca_full_cohort.tsv     : full cohort (1101) with labels + PAM50
  - tcga_brca_chemo_subset.tsv    : chemo-treated subset
  - extract_chemo_subset.log.json : QC + provenance
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(".")
PHASE3 = REPO / "." / "." / "results" / "phase_3"
CLIN = PHASE3 / "tcga_brca_clinical_merged.tsv"

CHEMO_DRUGS = {
    "DOXORUBICIN", "ADRIAMYCIN",
    "PACLITAXEL", "TAXOL",
    "CYCLOPHOSPHAMIDE", "CYTOXAN",
    "5-FLUOROURACIL", "FLUOROURACIL", "5FU",
    "EPIRUBICIN",
    "DOCETAXEL", "TAXOTERE",
    "DOXIL",
}

FIVE_YEARS = 5 * 365.25
THREE_YEARS = 3 * 365.25


def normalize_drug_token(t: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", t.upper())


def any_chemo(drug_list: str) -> bool:
    if not isinstance(drug_list, str) or not drug_list:
        return False
    tokens = [normalize_drug_token(x) for x in drug_list.split("|")]
    needles = {normalize_drug_token(x) for x in CHEMO_DRUGS}
    return any(t in needles for t in tokens)


def matched_chemo(drug_list: str) -> str:
    if not isinstance(drug_list, str):
        return ""
    needles = {normalize_drug_token(x) for x in CHEMO_DRUGS}
    out = []
    for tok in drug_list.split("|"):
        if normalize_drug_token(tok) in needles:
            out.append(tok)
    return "|".join(sorted(set(out)))


def parse_float(v):
    try:
        if pd.isna(v):
            return np.nan
        return float(v)
    except Exception:
        return np.nan


def build_os(row) -> tuple[int, float]:
    vital = str(row.get("vital_status", "")).strip().upper()
    if vital == "DECEASED":
        return 1, parse_float(row.get("days_to_death"))
    elif vital == "LIVING":
        return 0, parse_float(row.get("days_to_last_followup"))
    return -1, np.nan  # unknown


def build_nte(row) -> tuple[int, float]:
    nte = str(row.get("new_tumor_event_after_initial_treatment", "")).strip().upper()
    days = parse_float(row.get("days_to_new_tumor_event_after_initial_treatment"))
    if nte == "YES":
        return 1, days
    elif nte == "NO":
        return 0, np.nan
    return -1, np.nan


def label_window(row, window_days: float) -> int:
    """Binarized OS-or-recurrence outcome at `window_days`. Returns -1 if indeterminate."""
    if row["OS_event"] == 1 and np.isfinite(row["OS_time"]) and row["OS_time"] <= window_days:
        return 1
    if row["NTE_event"] == 1 and np.isfinite(row["NTE_time"]) and row["NTE_time"] <= window_days:
        return 1
    if row["OS_event"] == 0 and np.isfinite(row["OS_time"]) and row["OS_time"] >= window_days:
        # alive past window AND no within-window recurrence on record
        if not (row["NTE_event"] == 1 and np.isfinite(row["NTE_time"]) and row["NTE_time"] <= window_days):
            return 0
    return -1


def annotate_outcomes(df: pd.DataFrame) -> pd.DataFrame:
    """Add OS / NTE / 3yr+5yr labels and PAM50 column to a clinical subframe."""
    df = df.copy()
    os_event_time = df.apply(lambda r: build_os(r), axis=1, result_type="expand")
    df["OS_event"] = os_event_time[0]
    df["OS_time"] = os_event_time[1]
    nte = df.apply(lambda r: build_nte(r), axis=1, result_type="expand")
    df["NTE_event"] = nte[0]
    df["NTE_time"] = nte[1]
    df["label_3yr"] = df.apply(lambda r: label_window(r, THREE_YEARS), axis=1)
    df["label_5yr"] = df.apply(lambda r: label_window(r, FIVE_YEARS), axis=1)
    df["PAM50"] = df["PAM50Call_RNAseq"].astype(str).where(
        df["PAM50Call_RNAseq"].notna(), other=np.nan
    )
    return df


def select_cols(df: pd.DataFrame, extras: list[str]) -> pd.DataFrame:
    base = [
        "patient_id",
        "OS_event", "OS_time", "NTE_event", "NTE_time",
        "label_3yr", "label_5yr", "PAM50",
        "vital_status", "days_to_last_followup", "days_to_death",
        "new_tumor_event_after_initial_treatment",
        "days_to_new_tumor_event_after_initial_treatment",
    ]
    cols = [c for c in (base + extras) if c in df.columns]
    return df[cols].copy()


def main() -> int:
    if not CLIN.exists():
        print(f"FATAL: {CLIN} missing — run locate_tcga_brca.py first", file=sys.stderr)
        return 2

    clin = pd.read_csv(CLIN, sep="\t", low_memory=False)
    n_total = len(clin)

    # Annotate the full cohort.
    full = annotate_outcomes(clin)
    full_out = select_cols(full, ["drug_name_list"])
    full_path = PHASE3 / "tcga_brca_full_cohort.tsv"
    full_out.to_csv(full_path, sep="\t", index=False)

    # Build chemo flag (drug record + neoadjuvant).
    full["is_chemo_drug"] = full["drug_name_list"].apply(any_chemo)
    full["matched_chemo_drugs"] = full["drug_name_list"].apply(matched_chemo)
    full["is_chemo_neoadj"] = (
        full.get("history_of_neoadjuvant_treatment", pd.Series([np.nan]*len(full)))
        .astype(str).str.strip().str.upper() == "YES"
    )
    full["is_chemo"] = full["is_chemo_drug"] | full["is_chemo_neoadj"]
    chemo = full[full["is_chemo"]].copy()
    chemo_out = select_cols(chemo, ["matched_chemo_drugs", "drug_name_list",
                                     "is_chemo_drug", "is_chemo_neoadj"])
    chemo_path = PHASE3 / "tcga_brca_chemo_subset.tsv"
    chemo_out.to_csv(chemo_path, sep="\t", index=False)

    # QC summary.
    def label_breakdown(df: pd.DataFrame, lbl: str) -> dict:
        return df[lbl].value_counts(dropna=False).to_dict()

    def labeled_pam50(df: pd.DataFrame, lbl: str) -> dict:
        sub = df[df[lbl] != -1]
        return sub["PAM50"].value_counts(dropna=False).to_dict()

    def crosstab(df: pd.DataFrame, lbl: str) -> dict:
        sub = df[df[lbl] != -1]
        ct = pd.crosstab(sub["PAM50"].fillna("NaN"), sub[lbl])
        return ct.to_dict()

    qc = {
        "n_total_participants": int(n_total),
        "full_cohort": {
            "n": int(len(full_out)),
            "label_3yr_distribution": label_breakdown(full_out, "label_3yr"),
            "label_5yr_distribution": label_breakdown(full_out, "label_5yr"),
            "label_3yr_pam50_distribution": labeled_pam50(full_out, "label_3yr"),
            "label_5yr_pam50_distribution": labeled_pam50(full_out, "label_5yr"),
            "label_3yr_pam50_x_label": crosstab(full_out, "label_3yr"),
        },
        "chemo_subset": {
            "n": int(len(chemo_out)),
            "label_3yr_distribution": label_breakdown(chemo_out, "label_3yr"),
            "label_5yr_distribution": label_breakdown(chemo_out, "label_5yr"),
            "label_3yr_pam50_distribution": labeled_pam50(chemo_out, "label_3yr"),
            "matched_chemo_drug_counts": chemo_out["matched_chemo_drugs"].fillna("").value_counts().head(15).to_dict(),
        },
    }
    log = {
        "stage": "phase_3_C2_extract_chemo_subset",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(CLIN),
        "outputs": {
            "tcga_brca_full_cohort.tsv": str(full_path),
            "tcga_brca_chemo_subset.tsv": str(chemo_path),
        },
        "chemo_drug_filter": sorted(CHEMO_DRUGS),
        "windows_days": {"three_years": THREE_YEARS, "five_years": FIVE_YEARS},
        "qc": qc,
        "deviation_note": (
            "Original Phase 3 plan called for chemo-only subset. TCGA-BRCA chemo "
            "subset (n=120) is severely under-powered for OS/recurrence: only 10 "
            "events at 3y, median follow-up ~2.5y. PRIMARY worked example is "
            "therefore the FULL TCGA-BRCA cohort at 3-yr OS-or-recurrence "
            "(n_labeled=446, with each PAM50 stratum ≥31). Chemo subset retained "
            "as a sensitivity in the same file."
        ),
    }
    with open(PHASE3 / "extract_chemo_subset.log.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    print(json.dumps(qc, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
