"""Head-to-head benchmark: oncoPredict/pRRophetic-style ridge regression.

Both oncoPredict and pRRophetic share the same core algorithm:
  1. Train ridge regression: drug_AUC ~ gene_expression on CCLE/GDSC training cell lines
  2. Apply to target cohort expression to get predicted drug_AUC per patient
  3. Lower predicted AUC = more sensitive = predicted responder

oncoPredict adds power-transformation; pRRophetic adds ComBat batch correction.
Since both packages are unmaintained (R 4.5 incompatible / github 404), we
implement the core ridge algorithm in Python (sklearn).

Two variants:
  (a) rank_features: per-patient rank-based features (platform-invariant; matches
      our framework's design philosophy)
  (b) zscore_features: per-gene z-scoring across each platform separately, then
      ridge — pRRophetic-like (assumes log-scale absolute expression)

Output: per-cohort predicted_auc.tsv (rows=patient_id, predicted_auc).
Lower predicted AUC → predicted sensitive → predicted pCR=1.
For comparison with pCR labels, we use -predicted_auc as the "score" so that
high score = predicted sensitive (matches singscore convention).
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge, RidgeCV
from sklearn.metrics import roc_auc_score

CCLE_DIR = Path("./CCLE25Q3")
REVISE_ROOT = Path(".")
EXT = REVISE_ROOT / "external_data"

CCLE_EXPR = CCLE_DIR / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_AUC = CCLE_DIR / "Drug_sensitivity_AUC_(CTD^2)_subsetted.csv"
DRUG_COL = "DOXORUBICIN (CTRP:36599)"


def load_ccle_expression():
    """Load CCLE expression. Rows=cell line, cols=genes. Strip 'GENE (ENTREZ)' to 'GENE'."""
    df = pd.read_csv(CCLE_EXPR, index_col=0)
    # First several columns are metadata; expression starts at "TSPAN6 (7105)"
    meta_cols = ["SequencingID", "ModelID", "IsDefaultEntryForModel", "ModelConditionID", "IsDefaultEntryForMC"]
    model_id = df["ModelID"]
    expr = df.drop(columns=[c for c in meta_cols if c in df.columns])
    # Strip "(NCBIID)" from gene names
    expr.columns = [c.split(" (")[0].upper() for c in expr.columns]
    # De-duplicate columns (keep first if duplicates exist)
    expr = expr.loc[:, ~expr.columns.duplicated()]
    expr.index = model_id.values
    expr.index.name = "ModelID"
    # Some ModelIDs have multiple expression entries (replicates); collapse via mean
    expr = expr.groupby(level=0).mean()
    return expr


def load_ccle_auc():
    df = pd.read_csv(CCLE_AUC, index_col=0)
    if DRUG_COL not in df.columns:
        raise SystemExit(f"Drug column not found: {DRUG_COL}")
    return df[DRUG_COL].dropna()


def per_sample_rank(expr):
    """Within-sample rank-percentile across genes (rows = samples)."""
    return expr.rank(axis=1, method="average", pct=True)


def per_gene_zscore(expr):
    """Per-gene z-score (col-wise) — pRRophetic style."""
    mu = expr.mean(axis=0)
    sd = expr.std(axis=0).replace(0, 1)
    return (expr - mu) / sd


def train_ridge(X, y, alphas=(0.01, 0.1, 1, 10, 100, 1000)):
    """RidgeCV with built-in CV alpha selection."""
    model = RidgeCV(alphas=alphas, store_cv_results=False)
    model.fit(X, y)
    return model


def predict_per_cohort(method_name, ccle_expr, ccle_auc, transform_fn, cohort_expr):
    # Find common genes
    common = sorted(set(ccle_expr.columns) & set(cohort_expr.columns))
    if len(common) < 2000:
        raise ValueError(f"Too few common genes: {len(common)}")
    # Restrict to common genes
    ccle_e = ccle_expr[common]
    cohort_e = cohort_expr[common]
    # Drop CCLE cell lines without AUC
    train_idx = sorted(set(ccle_e.index) & set(ccle_auc.index))
    X_train = transform_fn(ccle_e.loc[train_idx])
    y_train = ccle_auc.loc[train_idx].values

    # Impute NaN with column median (per-gene)
    X_train = X_train.fillna(X_train.median(axis=0))
    X_test = transform_fn(cohort_e)
    # Use training median for held-out NaN imputation
    train_median = X_train.median(axis=0)
    X_test = X_test.fillna(train_median)
    # Any remaining NaN (gene fully missing) -> fill with neutral value
    fallback = 0.5 if transform_fn is per_sample_rank else 0.0
    X_train = X_train.fillna(fallback)
    X_test = X_test.fillna(fallback)

    model = train_ridge(X_train.values, y_train)
    pred = model.predict(X_test.values)
    return pd.Series(pred, index=cohort_e.index, name=f"predicted_auc_{method_name}"), {
        "n_train_cell_lines": len(train_idx),
        "n_common_genes": len(common),
        "alpha_selected": float(model.alpha_),
        "intercept": float(model.intercept_),
        "predicted_auc_mean": float(np.mean(pred)),
        "predicted_auc_sd": float(np.std(pred)),
    }


COHORTS = ["GSE16446", "GSE25066", "GSE22226"]


def main():
    print("Loading CCLE expression + AUC ...")
    ccle_expr = load_ccle_expression()
    ccle_auc = load_ccle_auc()
    print(f"  CCLE expr: {ccle_expr.shape[0]} cell lines × {ccle_expr.shape[1]} genes")
    print(f"  CCLE Doxorubicin AUC: {len(ccle_auc)} cell lines with data")

    out = {}
    for cohort in COHORTS:
        print(f"\n=== {cohort} ===")
        expr = pd.read_csv(EXT / cohort / f"{cohort}_expression.tsv", sep="\t", index_col=0)
        expr.columns = [c.upper() for c in expr.columns]
        pcr_df = pd.read_csv(EXT / cohort / f"{cohort}_pcr.tsv", sep="\t").dropna(subset=["pcr"])
        common_patients = sorted(set(expr.index.astype(str)) & set(pcr_df["patient_id"].astype(str)))
        expr = expr.loc[common_patients]
        pcr = pcr_df.set_index("patient_id").loc[common_patients]["pcr"].astype(int).to_numpy()

        cohort_out = {"n": len(common_patients), "pcr_rate": float(pcr.mean()), "methods": {}}

        for method_name, transform in [
            ("ridge_rank", per_sample_rank),
            ("ridge_zscore", per_gene_zscore),
        ]:
            print(f"  Method: {method_name} ...")
            pred, info = predict_per_cohort(method_name, ccle_expr, ccle_auc, transform, expr)
            # Score = -predicted_auc (high score = predicted sensitive = predicted pCR=1)
            score = -pred.values
            auc_pcr = float(roc_auc_score(pcr, score)) if len(np.unique(pcr)) >= 2 else float("nan")
            # Save predictions
            df_out = pd.DataFrame({"patient_id": common_patients, "predicted_auc": pred.values, "score": score, "pcr": pcr})
            df_out.to_csv(EXT / cohort / f"score_ridge_{method_name}.tsv", sep="\t", index=False, float_format="%.6g")
            info["auroc_predicted_auc_vs_pcr"] = auc_pcr
            cohort_out["methods"][method_name] = info
            print(f"    n_train={info['n_train_cell_lines']} alpha={info['alpha_selected']:.3g} "
                  f"pred_auc_mean={info['predicted_auc_mean']:.3f}  AUROC_vs_pCR={auc_pcr:.3f}")
        out[cohort] = cohort_out

    out_path = REVISE_ROOT / "primary_test_results" / "headtohead_ridge.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    # Summary table
    print("\n" + "=" * 80)
    print("HEAD-TO-HEAD SUMMARY: ridge regression vs framework vs proliferation")
    print("=" * 80)
    print(f"{'Cohort':<10} {'Method':<18} {'AUROC vs pCR':>15}")
    for c, info in out.items():
        for m, mi in info["methods"].items():
            print(f"{c:<10} {m:<18} {mi['auroc_predicted_auc_vs_pcr']:>15.3f}")
    print(f"\nResults: {out_path}")


if __name__ == "__main__":
    main()
