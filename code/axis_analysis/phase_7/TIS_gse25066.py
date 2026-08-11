"""Phase 7 / Sig 1 — TIS (Ayers 2017 18-gene IFN-gamma immune signature) x GSE25066.

Hypothesis (orthogonal-signature falsification test):
    The Ayers TIS encodes adaptive-immune / IFN-gamma activity (CD8 T-cell,
    antigen presentation, checkpoint axis) that is biologically orthogonal to
    the BC standard axes P (proliferation) / H (ER-PR) / B (ERBB2) /
    H_PAM50basal. Pre-spec: residual AUROC point > 0.60 AND CI low > 0.55
    → FALSIFIED (signal survives axis deconfounding = demonstrates the test
    has discrimination power).

Gene list (Ayers et al., J Clin Invest 2017, 127(8):2930-40):
    CCL5, CD27, CD274 (PD-L1), CD276 (B7-H3), CD8A, CMKLR1, CXCL9, CXCR6,
    GZMA, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2 (PD-L2),
    PSMB10, STAT1
    All same direction (immune-active = high score). High confidence; this is
    the published TIS reference list.

Score = mean z-score across the 18 genes after log expression (same
formula as Phase 6 _common compute_singscore with up-only set, but for
direct interpretability we use mean-z to match Ayers 2017).
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
EXPR_TSV = REPO / "external_data" / "GSE25066" / "GSE25066_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

TIS_GENES = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "GZMA", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7",
    "PDCD1LG2", "PSMB10", "STAT1",
]


def mean_zscore(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr[genes].copy()
    # per-gene z-score across samples
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def main():
    PHASE_7.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "TIS / Ayers 18-gene IFN-gamma immune signature (Ayers JCI 2017)",
        "signature_confidence": "HIGH (published verbatim 18-gene list)",
        "scoring": "mean per-gene z-score across 18 genes; all up (immune-active = high)",
        "n_requested": len(TIS_GENES),
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 -> CONFIRMED; >0.60 AND ci_lo>0.55 -> FALSIFIED; else GRAY",
        "expected_outcome": "HIGH chance of FALSIFICATION; TIS is biologically orthogonal to P/H/B/H_PAM50basal",
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    alias_map = load_alias_map(ALIAS)

    present, mapping, missing = resolve_genes(TIS_GENES, cohort_genes, alias_map)
    log["n_resolved"] = len(present)
    log["resolved_genes"] = present
    log["missing_genes"] = missing
    log["alias_mapping"] = mapping

    if len(present) < 10:
        log["error"] = f"too few TIS genes resolved ({len(present)}/{len(TIS_GENES)})"
        with open(PHASE_7 / "TIS_gse25066_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    score = mean_zscore(expr, present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_TIS": score.values})
    score_df.to_csv(PHASE_7 / "score_TIS_GSE25066.tsv", sep="\t", index=False, float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE25066.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    sub_raw = df.dropna(subset=["sig_TIS", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_TIS"].to_numpy(), sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_TIS",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_7 / "TIS_gse25066_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_TIS ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_7 / "TIS_gse25066_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 7 / Sig 1 TIS x GSE25066 (n={len(sub)}) =====")
    print(f"Resolved TIS genes: {len(present)}/{len(TIS_GENES)}")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
