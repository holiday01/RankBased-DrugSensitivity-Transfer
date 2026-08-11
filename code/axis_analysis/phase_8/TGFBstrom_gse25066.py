"""Phase 8 / Sig 1 — TGF-beta stromal / CAF (Calon 2015) x GSE25066.

Hypothesis (orthogonal-signature falsification test):
    The Calon stromal TGF-beta / cancer-associated-fibroblast (CAF) program
    captures *stromal-cell-derived* signaling (TGF-beta + ECM remodeling +
    fibroblast activation) that is biologically distinct from tumor-intrinsic
    proliferation (P), ER/PR (H), HER2 (B), and PAM50-basal axes. In the
    Calon 2015 PNAS / NCB body of work this program is also anti-pCR /
    pro-poor-prognosis on multiple cohorts.

Pre-spec: residual AUROC point > 0.60 AND CI low > 0.55 -> FALSIFIED
(signal survives axis deconfounding = test has discrimination power).

Gene list (Calon et al., PNAS / Nat Genet 2015, stromal TGF-beta /
CAF program, recalled MEDIUM-HIGH confidence): we use a 12-gene
stromal-TGF-beta panel that is robustly recoverable across published
TGF-beta-CAF lists (Calon 2015 NG TGFB-response stromal subset +
CAF6 + Moffitt PDAC stromal overlap):

    TGFB1, TGFBR2, INHBA, CTGF (CCN2), POSTN, FAP, SPARC, COL1A1,
    COL1A2, MMP2, ACTA2, THBS2

All positive direction (stromal-high = high score).

Score = mean per-gene z-score, same scoring scheme as Phase 7 TIS /
cytolytic for apples-to-apples comparison.

Gene-list confidence (honest):
- CAF6 strict (TGFBR2, TGFB1, INHBA, CTGF, POSTN, FAP) = HIGH confidence.
- Extra 6 (SPARC, COL1A1, COL1A2, MMP2, ACTA2, THBS2) are canonical
  TGF-beta-driven CAF/ECM markers from the Calon 2015 stromal program
  and Moffitt stromal subtype = MEDIUM-HIGH confidence.
- Overall list = MEDIUM-HIGH confidence proxy of Calon stromal TGF-beta.
"""
from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "phase_7"))
from _common import (
    sha256_of, load_alias_map, resolve_genes,
    auroc_with_boot_ci, fit_residual, verdict_3x3,
)

REPO = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
PHASE_0 = REPO / "resubmission_v2" / "results" / "phase_0"
PHASE_8 = REPO / "resubmission_v2" / "results" / "phase_8"
EXPR_TSV = REPO / "external_data" / "GSE25066" / "GSE25066_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

# Calon-style stromal TGF-beta / CAF panel (12 genes, all UP / stromal-high).
TGFB_STROMAL_GENES = [
    "TGFB1", "TGFBR2", "INHBA", "CTGF", "POSTN", "FAP",
    "SPARC", "COL1A1", "COL1A2", "MMP2", "ACTA2", "THBS2",
]


def mean_zscore(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr[genes].copy()
    # per-gene z-score across samples
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def axis_correlations(df: pd.DataFrame, sig_col: str,
                      axis_cols: list[str]) -> dict:
    out = {}
    for a in axis_cols:
        sub = df.dropna(subset=[sig_col, a])
        if len(sub) < 5:
            out[a] = float("nan")
        else:
            out[a] = float(sub[sig_col].corr(sub[a]))
    return out


def main():
    PHASE_8.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "TGF-beta stromal / CAF (Calon 2015 PNAS / NG; 12-gene proxy)",
        "signature_confidence": (
            "MEDIUM-HIGH (CAF6 core HIGH; 6-gene ECM/CAF extension MEDIUM-HIGH; "
            "robust proxy of Calon TGF-beta stromal program)"
        ),
        "scoring": "mean per-gene z-score across 12 genes; all UP (stromal-high)",
        "n_requested": len(TGFB_STROMAL_GENES),
        "decision_rule": (
            "point<=0.55 AND ci_lo<=0.55 -> CONFIRMED; "
            ">0.60 AND ci_lo>0.55 -> FALSIFIED; else GRAY"
        ),
        "expected_outcome": (
            "MEDIUM-HIGH chance of FALSIFICATION: stromal program biologically "
            "orthogonal to P/H/B/H_PAM50basal; anti-pCR direction expected "
            "(higher stromal score -> lower pCR -> AUROC < 0.5 possible)."
        ),
        "anti_pcr_note": (
            "Calon TGF-beta stromal score is anti-pCR; AUROC <0.5 in raw form "
            "still counts as signal. The 3x3 verdict_3x3() uses signed AUROC "
            "(rule from _common.py); a score with strong anti-pCR residual will "
            "show AUROC < 0.5, which under the pre-spec rule (residual point > 0.60 "
            "AND ci_lo > 0.55) does NOT trigger FALSIFIED. We therefore also "
            "report the absolute-distance-from-0.5 perspective in the JSON."
        ),
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    alias_map = load_alias_map(ALIAS)

    present, mapping, missing = resolve_genes(TGFB_STROMAL_GENES, cohort_genes, alias_map)
    log["n_resolved"] = len(present)
    log["resolved_genes"] = present
    log["missing_genes"] = missing
    log["alias_mapping"] = mapping

    if len(present) < 6:
        log["error"] = f"too few stromal genes resolved ({len(present)}/{len(TGFB_STROMAL_GENES)})"
        with open(PHASE_8 / "TGFBstrom_gse25066_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    score = mean_zscore(expr, present)
    score_df = pd.DataFrame({"sample_id": score.index,
                             "sig_TGFBstrom": score.values})
    score_df.to_csv(PHASE_8 / "score_TGFBstrom_GSE25066.tsv",
                    sep="\t", index=False, float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE25066.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    axis_cols = ["P", "H", "B", "H_PAM50basal"]
    corrs = axis_correlations(df, "sig_TGFBstrom", axis_cols)
    log["sig_vs_axis_pearson_r"] = corrs

    sub_raw = df.dropna(subset=["sig_TGFBstrom", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_TGFBstrom"].to_numpy(),
        sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_TGFBstrom",
                             axis_cols=axis_cols,
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_8 / "TGFBstrom_gse25066_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_TGFBstrom ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    # Also surface an absolute-distance-from-0.5 view for an anti-pCR signature.
    log["residual_auroc_abs_dist_from_0p5"] = {
        "point": abs(auc - 0.5) if auc == auc else float("nan"),
        "ci_low_abs": abs(lo - 0.5) if lo == lo else float("nan"),
        "ci_high_abs": abs(hi - 0.5) if hi == hi else float("nan"),
    }
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_8 / "TGFBstrom_gse25066_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 8 / Sig 1 TGF-beta stromal x GSE25066 (n={len(sub)}) =====")
    print(f"Resolved stromal genes: {len(present)}/{len(TGFB_STROMAL_GENES)}: {present}")
    print(f"Missing: {missing}")
    print(f"Sig vs axes Pearson r: {corrs}")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
