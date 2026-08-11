"""Phase 7 / Sig 3 — DLDA30-like (Hatzis 2011) x GSE16446.

Hypothesis: The Hatzis 2011 DLDA30 TFAC-pCR predictor was trained as a
multi-class, drug-specific (TFAC) classifier on MDACC patients. Pre-spec:
residual AUROC point > 0.60 AND CI low > 0.55 -> FALSIFIED (drug-specific
signal orthogonal to standard P/H/B axes).

Gene list - HONEST CAVEAT
-------------------------
The exact 30 DLDA-selected probesets from Hatzis 2011 (JAMA 305:1873-81;
Suppl. Table 2) are not memorized verbatim. The gene set below is a
26-gene approximation distilled from the most reliably recalled members
of the DLDA30 panel. It mixes:
  - proliferation/mitotic core (PROL+; UP)
  - ER-axis differentiation (ER+; DOWN)
matching the published DLDA30 signed-direction structure (proliferation
positive for pCR, ER-axis negative for pCR in TFAC).

The DLDA30 was trained on an Affymetrix U133A MDACC cohort independent of
GSE16446 (Bonnefoi-EORTC FEC100 single-agent doxorubicin TNBC). GSE16446
is a single-stratum-fallback cohort (all ER-), which is also a stress test
for any ER-axis component of the score.

Approximate gene list (26 genes, +/- direction)
------------------------------------------------
PROL+ / mitotic (UP, weight +1):
  BTG3, CCNE2, CDC2 (CDK1), FOXM1, MELK, MKI67, MYBL2, RACGAP1, RRM2, TYMS,
  ZWINT, STMN1, AURKA, TRIP13, NUSAP1, KIF11, TPX2, BIRC5, EZH2, UHRF1,
  DLGAP5, GINS1, MCM10
ER+ / luminal differentiation (DOWN, weight -1):
  STC2, TFF3, IL6ST
Score = mean(UP z) - mean(DOWN z)  [signed-z sum form]

This is documented in the JSON output as "DLDA30-like approximation".
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
    sha256_of, load_alias_map, resolve_genes,
    auroc_with_boot_ci, fit_residual, verdict_3x3,
)

REPO = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
PHASE_0 = REPO / "resubmission_v2" / "results" / "phase_0"
PHASE_7 = REPO / "resubmission_v2" / "results" / "phase_7"
EXPR_TSV = REPO / "external_data" / "GSE16446" / "GSE16446_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

DLDA30_UP = [
    # proliferation + mitotic
    "BTG3", "CCNE2", "CDK1", "FOXM1", "MELK", "MKI67", "MYBL2", "RACGAP1",
    "RRM2", "TYMS", "ZWINT", "STMN1", "AURKA", "TRIP13", "NUSAP1", "KIF11",
    "TPX2", "BIRC5", "EZH2", "UHRF1", "DLGAP5", "GINS1", "MCM10",
]
DLDA30_DOWN = [
    # ER-axis luminal differentiation
    "STC2", "TFF3", "IL6ST",
]


def mean_zscore(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr[genes].copy()
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def main():
    PHASE_7.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "DLDA30-like (Hatzis JAMA 2011) - 26-gene approximation",
        "signature_confidence": "LOW-MEDIUM (gene list reconstructed from memory; not verbatim from Hatzis 2011 Suppl. Table 2)",
        "signature_caveat": (
            "23 UP (PROL+/mitotic) + 3 DOWN (ER+/luminal). True DLDA30 has "
            "30 probesets with DLDA weights; we use signed-z-sum with +/-1. "
            "This captures the dominant proliferation-vs-ER-direction signal "
            "but loses drug-specific weights."
        ),
        "scoring": "mean(UP z) - mean(DOWN z); signed-z-sum form",
        "n_up_requested": len(DLDA30_UP),
        "n_down_requested": len(DLDA30_DOWN),
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 -> CONFIRMED; >0.60 AND ci_lo>0.55 -> FALSIFIED; else GRAY",
        "expected_outcome": "MEDIUM-HIGH chance of FALSIFICATION; drug-specific (TFAC) trained on different MDACC cohort",
        "cohort_caveat": "GSE16446 is all-ER- (single-stratum-fallback); ER-axis DOWN component may be flat",
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    alias_map = load_alias_map(ALIAS)

    up_present, up_map, up_missing = resolve_genes(DLDA30_UP, cohort_genes, alias_map)
    down_present, down_map, down_missing = resolve_genes(DLDA30_DOWN, cohort_genes, alias_map)
    overlap = set(up_present) & set(down_present)
    if overlap:
        down_present = [g for g in down_present if g not in overlap]
    log["n_up_resolved"] = len(up_present)
    log["n_down_resolved"] = len(down_present)
    log["up_resolved"] = up_present
    log["down_resolved"] = down_present
    log["up_missing"] = up_missing
    log["down_missing"] = down_missing

    if len(up_present) < 10:
        log["error"] = f"too few UP genes resolved ({len(up_present)})"
        with open(PHASE_7 / "DLDA30_gse16446_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    up_score = mean_zscore(expr, up_present)
    if down_present:
        down_score = mean_zscore(expr, down_present)
    else:
        down_score = pd.Series(0.0, index=expr.index)
    score = up_score - down_score
    score_df = pd.DataFrame({"sample_id": score.index, "sig_DLDA30": score.values})
    score_df.to_csv(PHASE_7 / "score_DLDA30_GSE16446.tsv", sep="\t", index=False, float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE16446.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    if df["pcr"].notna().sum() > 0:
        log["pcr_rate"] = float((df["pcr"] == 1).mean())

    sub_raw = df.dropna(subset=["sig_DLDA30", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_DLDA30"].to_numpy(), sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_DLDA30",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_7 / "DLDA30_gse16446_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_DLDA30 ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_7 / "DLDA30_gse16446_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 7 / Sig 3 DLDA30-like x GSE16446 (n={len(sub)}) =====")
    print(f"UP resolved: {len(up_present)}/{len(DLDA30_UP)}; DOWN resolved: {len(down_present)}/{len(DLDA30_DOWN)}")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
