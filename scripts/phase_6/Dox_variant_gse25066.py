"""Phase 6 / Sig 3 — Doxorubicin SIGNATURE-VARIANT × GSE25066.

We cannot run actual pRRophetic (R package, CCLE + IC50 ridge regression)
within this autonomous Python pipeline. As a substitute, we use the
**alternate frozen Doxorubicin signature** from
`frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv`. This signature is
also CCLE-Pearson-derived but uses a DIFFERENT gene-universe filter
(quantile_20 = top-quantile expression cutoff; protein-coding only)
relative to the headline framework used in Phase 2 (which used the
`solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv` signature). So the
analysis here is best described as **"Dox-signature variant comparison"**
rather than a true pRRophetic comparison.

Pre-specified pipeline (mirrors Phase 2):
    sig_DoxVariant ~ P + H + B + H_PAM50basal  →  residual vs pCR.
3×3 verdict.
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

REPO = Path(".")
PHASE_0 = REPO / "." / "results" / "phase_0"
PHASE_6 = REPO / "." / "results" / "phase_6"
EXPR_TSV = REPO / "external_data" / "GSE25066" / "GSE25066_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"
SIG_TSV = REPO / "frozen_signatures" / "quantile_20_protein" / "DOXORUBICIN.tsv"


def main():
    PHASE_6.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature_path": str(SIG_TSV),
        "signature_sha256": sha256_of(SIG_TSV),
        "signature_label": "DoxVariant_quantile20_protein (pRRophetic stand-in)",
        "signature_caveat": (
            "Not actual pRRophetic. Same CCLE-Pearson scaffold but quantile_20 "
            "expression-cutoff + protein-coding gene-universe filter."
        ),
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 → CONFIRMED; >0.60 AND ci_lo>0.55 → FALSIFIED; else GRAY",
    }

    # Load signature, partition by direction
    sig = pd.read_csv(SIG_TSV, sep="\t")
    up_genes = sig.loc[sig["direction"] == "sensitivity", "gene_symbol"].astype(str).tolist()
    down_genes = sig.loc[sig["direction"] == "resistance", "gene_symbol"].astype(str).tolist()
    log["n_sensitivity_in_sig"] = len(up_genes)
    log["n_resistance_in_sig"] = len(down_genes)

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    expr.index = expr.index.astype(str)
    alias_map = load_alias_map(ALIAS)

    up_present, _, up_missing = resolve_genes(up_genes, set(expr.columns), alias_map)
    down_present, _, down_missing = resolve_genes(down_genes, set(expr.columns), alias_map)
    overlap = set(up_present) & set(down_present)
    if overlap:
        down_present = [g for g in down_present if g not in overlap]
    log["n_up_resolved"] = len(up_present)
    log["n_down_resolved"] = len(down_present)
    log["up_missing"] = up_missing
    log["down_missing"] = down_missing

    score = compute_singscore(expr, up_present, down_present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_DoxVariant": score.values})
    score_df.to_csv(PHASE_6 / "score_DoxVariant_GSE25066.tsv", sep="\t", index=False,
                    float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE25066.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    sub_raw = df.dropna(subset=["sig_DoxVariant", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_DoxVariant"].to_numpy(), sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_DoxVariant",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_6 / "DoxVariant_gse25066_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_DoxVariant ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_6 / "DoxVariant_gse25066_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 6 / Sig 3 DoxVariant × GSE25066 (n={len(sub)}) =====")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
