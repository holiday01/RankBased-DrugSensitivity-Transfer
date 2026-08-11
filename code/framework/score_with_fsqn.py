"""Pre-reg D1 sensitivity analysis: FSQN-aligned scoring.

Applies Feature-Specific Quantile Normalization (FSQN) to align the cohort's
per-gene expression distribution to a CCLE reference (solid-only) BEFORE
computing within-sample rank-based singscore.

This tests whether rank-only's platform-invariance claim holds: if FSQN gives
substantially different AUROCs from rank-only, the platform-invariance claim
must be qualified in the manuscript.

CCLE reference pool = solid-only (heme cell lines excluded, matches signature
derivation pool).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


REPO_ROOT = Path("/home/holiday01/2026_ISMB_code")
CCLE_EXPR_PATH = REPO_ROOT / "CCLE25Q3" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_MODEL_PATH = REPO_ROOT / "CCLE25Q3" / "Model.csv"
CCLE_DRUG_PATH = REPO_ROOT / "enrichment_test" / "ccle_drug_common_with_tcga.csv"


def clean_gene(g):
    return re.sub(r"\s*\(.*?\)", "", str(g)).upper().strip()


def load_ccle_solid_reference() -> pd.DataFrame:
    """Load CCLE solid-only expression (samples x genes)."""
    ccle_drug = pd.read_csv(CCLE_DRUG_PATH)
    common_models = ccle_drug.iloc[:, 0].astype(str)
    ccle_exp = pd.read_csv(CCLE_EXPR_PATH)
    expr = ccle_exp[ccle_exp.ModelID.isin(common_models)].copy().set_index("ModelID")
    meta_cols = ["Unnamed: 0", "SequencingID", "IsDefaultEntryForModel",
                 "ModelConditionID", "IsDefaultEntryForMC"]
    expr = expr.drop(columns=[c for c in meta_cols if c in expr.columns], errors="ignore")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr[~expr.index.duplicated(keep="first")]
    expr.columns = [clean_gene(g) for g in expr.columns]
    expr = expr.T.groupby(level=0).mean().T
    # Restrict to solid only
    model = pd.read_csv(CCLE_MODEL_PATH, low_memory=False)
    heme_keys = ["myeloid", "lymphoid", "plasma cell"]
    heme_dis = ["leukemia", "lymphoma", "myeloma", "myelodysplastic"]
    lin = model["OncotreeLineage"].astype(str).str.lower().fillna("")
    dis = model["OncotreePrimaryDisease"].astype(str).str.lower().fillna("")
    heme_ids = set(model.loc[
        lin.apply(lambda s: any(k in s for k in heme_keys)) |
        dis.apply(lambda s: any(k in s for k in heme_dis)),
        "ModelID"].astype(str))
    return expr.loc[~expr.index.isin(heme_ids)]


def fsqn_transform(ref: np.ndarray, target: np.ndarray) -> np.ndarray:
    """For one gene: map target values to ref's empirical quantile distribution.

    For each target value, find its empirical quantile in target distribution,
    then map to the same quantile in ref distribution.
    """
    ref_sorted = np.sort(ref[~np.isnan(ref)])
    if len(ref_sorted) == 0:
        return np.full_like(target, np.nan, dtype=float)
    target = np.asarray(target, dtype=float)
    out = np.full_like(target, np.nan, dtype=float)
    valid = ~np.isnan(target)
    if not valid.any():
        return out
    tv = target[valid]
    # Empirical CDF of target -> quantile in [0, 1]
    order = np.argsort(tv, kind="mergesort")
    quantiles = (np.arange(len(tv)) + 0.5) / len(tv)
    # Map quantile to ref's empirical quantile
    mapped = np.percentile(ref_sorted, quantiles * 100.0)
    out_valid = np.empty_like(tv)
    out_valid[order] = mapped
    out[valid] = out_valid
    return out


def fsqn_align(ref_expr: pd.DataFrame, target_expr: pd.DataFrame) -> pd.DataFrame:
    """Align target_expr (sample x gene) to ref_expr's per-gene distributions."""
    common_genes = sorted(set(ref_expr.columns) & set(target_expr.columns))
    aligned = pd.DataFrame(index=target_expr.index, columns=common_genes, dtype=float)
    for g in common_genes:
        aligned[g] = fsqn_transform(ref_expr[g].to_numpy(), target_expr[g].to_numpy())
    return aligned


def singscore(expr: pd.DataFrame, up: set, down: set) -> pd.Series:
    if expr.shape[1] < 3000:
        raise ValueError(f"Too few genes: {expr.shape[1]}")
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    up_p = sorted(up & set(expr.columns))
    down_p = sorted(down & set(expr.columns))
    if len(up_p) < 3 or len(down_p) < 3:
        return pd.Series(index=expr.index, dtype=float)

    def n(m, k):
        mn = (1 + k) / (2 * n_total); mx = 1 - mn
        return 2 * (m - mn) / (mx - mn) - 1

    raw = n(ranks[up_p].mean(1), len(up_p)) - n((1 - ranks[down_p]).mean(1), len(down_p))
    return raw - raw.median()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signature", required=True, type=Path)
    ap.add_argument("--expression", required=True, type=Path)
    ap.add_argument("--labels", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    sig = pd.read_csv(args.signature, sep="\t")
    up = set(sig[sig["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    down = set(sig[sig["direction"] == "resistance"]["gene_symbol"].str.upper())

    target = pd.read_csv(args.expression, sep="\t", index_col=0)
    target.columns = [c.upper() for c in target.columns]
    labels_df = pd.read_csv(args.labels, sep="\t")
    common = sorted(set(target.index.astype(str)) & set(labels_df["patient_id"].astype(str)))
    target = target.loc[common]
    y = labels_df.set_index("patient_id").loc[common]["pcr"].astype(int).to_numpy()

    print(f"Loading CCLE solid-only reference ...", flush=True)
    ref = load_ccle_solid_reference()
    print(f"  CCLE solid: {ref.shape}", flush=True)

    print(f"FSQN-aligning {args.expression.parent.name} to CCLE reference ...", flush=True)
    aligned = fsqn_align(ref, target)
    print(f"  Aligned shape: {aligned.shape}", flush=True)

    score_rankonly = singscore(target, up, down)
    score_fsqn = singscore(aligned, up, down)

    valid = score_rankonly.dropna().index.intersection(score_fsqn.dropna().index)
    s_rank = score_rankonly.loc[valid].to_numpy()
    s_fsqn = score_fsqn.loc[valid].to_numpy()
    y_valid = labels_df.set_index("patient_id").loc[valid]["pcr"].astype(int).to_numpy()

    auc_rank = roc_auc_score(y_valid, s_rank) if len(np.unique(y_valid)) >= 2 else float("nan")
    auc_fsqn = roc_auc_score(y_valid, s_fsqn) if len(np.unique(y_valid)) >= 2 else float("nan")
    pearson = float(np.corrcoef(s_rank, s_fsqn)[0, 1])
    spearman_rho = float(pd.Series(s_rank).corr(pd.Series(s_fsqn), method="spearman"))

    out_df = pd.DataFrame({
        "patient_id": list(valid),
        "score_rankonly": s_rank,
        "score_fsqn": s_fsqn,
        "pcr": y_valid,
    })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False, float_format="%.10g")

    summary = {
        "signature": str(args.signature),
        "cohort": args.expression.parent.name,
        "n_patients": int(len(valid)),
        "auc_rankonly": float(auc_rank),
        "auc_fsqn": float(auc_fsqn),
        "auc_difference": float(auc_fsqn - auc_rank),
        "pearson_rank_vs_fsqn": pearson,
        "spearman_rank_vs_fsqn": spearman_rho,
        "interpretation": (
            "If auc_difference is large (>0.05), platform-invariance claim of rank-only is qualified. "
            "If small (<0.02) and Pearson > 0.9, rank-only is justified."
        ),
    }
    with open(args.output.with_suffix(".json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
