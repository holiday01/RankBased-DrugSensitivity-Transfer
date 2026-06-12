"""Phase 7 / Sig 2 — GZMA + PRF1 cytolytic 2-gene (Rooney 2015) x GSE22226.

Hypothesis: The Rooney cytolytic score (geometric mean of GZMA + PRF1) captures
CD8 T-cell + NK effector activity, orthogonal to BC P/H/B/H_PAM50basal axes.
Pre-spec: residual AUROC point > 0.60 AND CI low > 0.55 -> FALSIFIED.

Gene list (Rooney et al., Cell 2015, 160(1-2):48-61):
    GZMA, PRF1
Score = mean z-score across the 2 genes after log expression. (Rooney
defines cytolytic activity as geometric mean of TPM expression; on
log-normalized microarray data the closest analogue is mean of log values,
implemented here as mean z-score for cross-sample comparability.)
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

REPO = Path(".")
PHASE_0 = REPO / "." / "results" / "phase_0"
PHASE_7 = REPO / "." / "results" / "phase_7"
EXPR_TSV = REPO / "external_data" / "GSE22226" / "GSE22226_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

CYT_GENES = ["GZMA", "PRF1"]


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
        "signature": "GZMA + PRF1 cytolytic 2-gene (Rooney Cell 2015)",
        "signature_confidence": "HIGH (canonical 2-gene definition)",
        "scoring": "mean per-gene z-score across the 2 genes (log proxy for Rooney geometric-mean cytolytic activity)",
        "n_requested": len(CYT_GENES),
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 -> CONFIRMED; >0.60 AND ci_lo>0.55 -> FALSIFIED; else GRAY",
        "expected_outcome": "MEDIUM chance of FALSIFICATION; immune-effector axis orthogonal to P/H/B",
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    alias_map = load_alias_map(ALIAS)

    present, mapping, missing = resolve_genes(CYT_GENES, cohort_genes, alias_map)
    log["n_resolved"] = len(present)
    log["resolved_genes"] = present
    log["missing_genes"] = missing
    log["alias_mapping"] = mapping

    if len(present) < 2:
        log["error"] = f"need both GZMA and PRF1; got {len(present)}"
        with open(PHASE_7 / "cytolytic_gse22226_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    score = mean_zscore(expr, present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_cytolytic": score.values})
    score_df.to_csv(PHASE_7 / "score_cytolytic_GSE22226.tsv", sep="\t", index=False, float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE22226.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    sub_raw = df.dropna(subset=["sig_cytolytic", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_cytolytic"].to_numpy(), sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_cytolytic",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_7 / "cytolytic_gse22226_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_cytolytic ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_7 / "cytolytic_gse22226_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 7 / Sig 2 Cytolytic x GSE22226 (n={len(sub)}) =====")
    print(f"Resolved genes: {present}")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
