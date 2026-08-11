"""Phase 6 / Sig 2 — OncotypeDX RS-like 21-gene × GSE25066 HR+ subset.

Signature reconstruction
------------------------
Paik 2004 (NEJM 351:2817) 21-gene Recurrence Score panel. Exact ddCt weights
not implementable without raw probe-level data, so we use a SIGNED Z-SUM
approximation:

  Proliferation group (UP, +1):  MKI67, AURKA, BIRC5, CCNB1, MYBL2
  ER group         (DOWN, -1):   ESR1, PGR, BCL2, SCUBE2
  HER2 group       (UP, +1):     ERBB2, GRB7
  Invasion group   (UP, +1):     MMP11, CTSV (alias CTSL2)
  Other (UP, +1):                CD68, BAG1, GSTM1

Reference housekeeping (ACTB, RPLP0, GUSB, GAPDH, TFRC) are NOT used directly —
we use full-transcriptome-rank singscore so HK normalization is implicit.

The RS formula is approximated as a singscore with UP and DOWN groups as
listed. This captures the dominant signed direction of RS (proliferation -
ER) used in clinical decision-making.

Cohort: GSE25066 filtered to HR+ stratum (er_status_ihc == 'P').
"""
from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    sha256_of, load_alias_map, resolve_genes, compute_singscore,
    auroc_with_boot_ci, fit_residual, verdict_3x3,
)

REPO = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
PHASE_0 = REPO / "resubmission_v2" / "results" / "phase_0"
PHASE_6 = REPO / "resubmission_v2" / "results" / "phase_6"
EXPR_TSV = REPO / "external_data" / "GSE25066" / "GSE25066_expression.tsv"
PHENO_TSV = REPO / "external_data" / "GSE25066" / "GSE25066_phenotype.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

RS_UP = [
    # Proliferation (+1)
    "MKI67", "AURKA", "BIRC5", "CCNB1", "MYBL2",
    # HER2 (+1, with implicit 0.47 weight but treated as +1 in z-sum)
    "ERBB2", "GRB7",
    # Invasion (+1)
    "MMP11", "CTSV",
    # Other (+1)
    "CD68", "BAG1", "GSTM1",
]
RS_DOWN = [
    # ER group (-1 in formula, so DOWN for high-score = high recurrence)
    "ESR1", "PGR", "BCL2", "SCUBE2",
]


def main():
    PHASE_6.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "OncotypeDX RS-like (Paik 2004 21-gene approximation)",
        "signature_caveat": (
            "Signed singscore approximation; exact ddCt weights not used. "
            "HK genes not applied directly (full-transcriptome rank is implicit normalization). "
            "HER2 weight (0.47x in true RS) treated as +1 in z-sum."
        ),
        "n_up": len(RS_UP),
        "n_down": len(RS_DOWN),
        "cohort_filter": "GSE25066 HR+ subset (er_status_ihc == 'P')",
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 → CONFIRMED; >0.60 AND ci_lo>0.55 → FALSIFIED; else GRAY",
    }

    # Phenotype + filter to HR+. Note: expression matrix is indexed by GSM
    # accession which corresponds to `patient_id` in the phenotype TSV
    # (column `sample_id` there is the Hatzis study ID, not the GEO GSM).
    pheno = pd.read_csv(PHENO_TSV, sep="\t")
    hr_pos_samples = set(pheno.loc[pheno["er_status_ihc"] == "P", "patient_id"].astype(str))
    log["n_HRpos_total_phenotype"] = len(hr_pos_samples)

    # Expression — subset rows to HR+
    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    expr.index = expr.index.astype(str)
    expr_hr = expr.loc[expr.index.isin(hr_pos_samples)].copy()
    log["n_HRpos_with_expression"] = int(expr_hr.shape[0])

    alias_map = load_alias_map(ALIAS)
    up_present, _, up_missing = resolve_genes(RS_UP, set(expr_hr.columns), alias_map)
    down_present, _, down_missing = resolve_genes(RS_DOWN, set(expr_hr.columns), alias_map)
    overlap = set(up_present) & set(down_present)
    if overlap:
        down_present = [g for g in down_present if g not in overlap]
    log["up_resolved"] = up_present
    log["down_resolved"] = down_present
    log["up_missing"] = up_missing
    log["down_missing"] = down_missing

    score = compute_singscore(expr_hr, up_present, down_present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_RSlike": score.values})
    score_df.to_csv(PHASE_6 / "score_RSlike_GSE25066_HRpos.tsv", sep="\t", index=False,
                    float_format="%.6f")

    # Axes — use existing axes_GSE25066.tsv but restrict to HR+ subset
    axes = pd.read_csv(PHASE_0 / "axes_GSE25066.tsv", sep="\t")
    axes = axes[axes["sample_id"].astype(str).isin(hr_pos_samples)].copy()
    log["n_HRpos_axes"] = int(len(axes))

    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    sub_raw = df.dropna(subset=["sig_RSlike", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_RSlike"].to_numpy(), sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_RSlike",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_6 / "RSlike_gse25066_HRpos_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_RSlike ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_6 / "RSlike_gse25066_HRpos_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 6 / Sig 2 RS-like × GSE25066 HR+ (n={len(sub)}) =====")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
