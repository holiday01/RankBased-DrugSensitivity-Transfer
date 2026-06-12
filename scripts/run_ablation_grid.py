"""B5 ablation robustness grid (pre-reg §8b).

For one external cohort, derive ablation signatures (varying quantile threshold,
top-N, gene universe, up:down asymmetry, FDR gate), score the cohort with each
config, and verify that the pre-registered primary config ranks in the TOP
QUARTILE of AUROC across all grid configs.

Grid (per pre-reg §8b):
  - quantile threshold q ∈ {0.20, 0.25, 0.30}
  - top-N per direction ∈ {10, 20, 30, 50}
  - gene universe ∈ {OncoKB-all, OncoKB-onco-only, protein-coding}
  - up:down asymmetry ∈ {30:30, 30:10, 10:30}    (only for top-N>=30)
  - FDR gate ∈ {none, FDR<0.10, FDR<0.05}        (only for quantile path)

NOTE: this is a robustness check — not used for primary inference.

Usage:
  python run_ablation_grid.py \
      --drug DOXORUBICIN \
      --cohort GSE16446 \
      --expression external_data/GSE16446/GSE16446_expression.tsv \
      --labels external_data/GSE16446/GSE16446_pcr.tsv \
      --output external_data/GSE16446/ablation_grid_DOXORUBICIN.tsv
"""

from __future__ import annotations

import argparse
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind
from sklearn.metrics import roc_auc_score


REPO_ROOT = Path(".")
CCLE_DRUG_PATH = REPO_ROOT / "enrichment_test" / "ccle_drug_common_with_tcga.csv"
CCLE_EXPR_PATH = REPO_ROOT / "CCLE25Q3" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_MODEL_PATH = REPO_ROOT / "CCLE25Q3" / "Model.csv"
ONCOKB_PATH = REPO_ROOT / "OncoKB" / "cancerGeneList.tsv"


def clean_gene(g):
    return re.sub(r"\s*\(.*?\)", "", str(g)).upper().strip()


def norm_drug(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def load_solid_expr_and_drug(drug_norm: str):
    ccle_drug = pd.read_csv(CCLE_DRUG_PATH)
    ccle_exp = pd.read_csv(CCLE_EXPR_PATH)
    model = pd.read_csv(CCLE_MODEL_PATH, low_memory=False)
    oncokb = pd.read_csv(ONCOKB_PATH, sep="\t")

    heme_keys = ["myeloid", "lymphoid", "plasma cell"]
    heme_dis = ["leukemia", "lymphoma", "myeloma", "myelodysplastic"]
    lineage = model["OncotreeLineage"].astype(str).str.lower().fillna("")
    disease = model["OncotreePrimaryDisease"].astype(str).str.lower().fillna("")
    heme_ids = set(model.loc[
        lineage.apply(lambda s: any(k in s for k in heme_keys)) |
        disease.apply(lambda s: any(k in s for k in heme_dis)),
        "ModelID"].astype(str))

    common_models = ccle_drug.iloc[:, 0].astype(str)
    expr = ccle_exp[ccle_exp.ModelID.isin(common_models)].copy().set_index("ModelID")
    meta_cols = ["Unnamed: 0", "SequencingID", "IsDefaultEntryForModel",
                 "ModelConditionID", "IsDefaultEntryForMC"]
    expr = expr.drop(columns=[c for c in meta_cols if c in expr.columns], errors="ignore")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr[~expr.index.duplicated(keep="first")]
    expr.columns = [clean_gene(g) for g in expr.columns]
    expr = expr.T.groupby(level=0).mean().T
    expr = expr.loc[~expr.index.isin(heme_ids)]  # solid-only

    sym_col = "Hugo Symbol" if "Hugo Symbol" in oncokb.columns else "Gene Symbol"
    oncogenes_all = set(oncokb[sym_col].astype(str).str.upper())
    if "Is Oncogene" in oncokb.columns:
        oncogenes_only = set(oncokb[oncokb["Is Oncogene"] == "Yes"][sym_col].astype(str).str.upper())
    elif "Gene Type" in oncokb.columns:
        oncogenes_only = set(oncokb[oncokb["Gene Type"].str.upper().isin(["ONCOGENE", "ONCO"])][sym_col].astype(str).str.upper())
    else:
        oncogenes_only = oncogenes_all

    ccle_drug_map = {}
    for col in ccle_drug.columns[1:]:
        key = norm_drug(col.split("(")[0])
        if key not in ccle_drug_map:
            ccle_drug_map[key] = col
    if drug_norm not in ccle_drug_map:
        raise SystemExit(f"Drug {drug_norm} not in CCLE")
    auc = ccle_drug.set_index(ccle_drug.columns[0])[ccle_drug_map[drug_norm]].astype(float)
    return expr, auc, oncogenes_all, oncogenes_only


def derive_quantile_sig(expr, auc, q, top_n_up, top_n_down, gene_universe_set, fdr_gate):
    """Returns (up_genes, down_genes) or None."""
    common = expr.index.intersection(auc.index)
    auc_sub = auc.loc[common].dropna()
    if len(auc_sub) < 100:
        return None
    expr_sub = expr.loc[auc_sub.index]
    if gene_universe_set is not None:
        cols = [c for c in expr_sub.columns if c in gene_universe_set]
        if len(cols) < 50:
            return None
        expr_sub = expr_sub[cols]
    low = auc_sub.quantile(q)
    high = auc_sub.quantile(1 - q)
    sens_idx = auc_sub[auc_sub <= low].index
    res_idx = auc_sub[auc_sub >= high].index
    if len(sens_idx) < 50 or len(res_idx) < 50:
        return None
    sens_e = expr_sub.loc[sens_idx]
    res_e = expr_sub.loc[res_idx]
    diff = (sens_e.mean(0) - res_e.mean(0)).dropna()
    if fdr_gate is not None:
        # Compute per-gene Welch t-test p-value, BH-FDR-filter to <fdr_gate, then top-N within survivors
        from statsmodels.stats.multitest import multipletests
        p_vals = []
        genes = diff.index.tolist()
        for g in genes:
            try:
                p_vals.append(ttest_ind(sens_e[g].dropna(), res_e[g].dropna(),
                                        equal_var=False, nan_policy="omit").pvalue)
            except Exception:
                p_vals.append(1.0)
        p_arr = np.array(p_vals)
        rejected = multipletests(p_arr, method="fdr_bh", alpha=fdr_gate)[0]
        survivors = set(np.array(genes)[rejected])
        diff = diff.loc[diff.index.isin(survivors)]
        if len(diff) < max(top_n_up, top_n_down):
            return None
    up = sorted(diff.nlargest(top_n_up).index.tolist())
    down = sorted(diff.nsmallest(top_n_down).index.tolist())
    return set(up), set(down)


def singscore(expr_cohort: pd.DataFrame, up_set: set, down_set: set) -> pd.Series:
    if expr_cohort.shape[1] < 5000:
        raise ValueError(f"need full transcriptome (got {expr_cohort.shape[1]})")
    n_total = expr_cohort.shape[1]
    ranks = expr_cohort.rank(axis=1, method="average", pct=True)
    up_p = sorted(up_set & set(expr_cohort.columns))
    down_p = sorted(down_set & set(expr_cohort.columns))
    if len(up_p) < 3 or len(down_p) < 3:
        return pd.Series(index=expr_cohort.index, dtype=float)
    def n(m, k):
        mn = (1 + k) / (2 * n_total); mx = 1 - mn
        return 2 * (m - mn) / (mx - mn) - 1
    raw = n(ranks[up_p].mean(1), len(up_p)) - n((1 - ranks[down_p]).mean(1), len(down_p))
    return raw - raw.median()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--expression", required=True, type=Path)
    ap.add_argument("--labels", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()

    drug_norm = norm_drug(args.drug)
    print(f"Loading CCLE solid-only expression + {drug_norm} AUC ...", flush=True)
    expr_solid, auc, oncokb_all, oncokb_only = load_solid_expr_and_drug(drug_norm)

    cohort_expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    labels_df = pd.read_csv(args.labels, sep="\t")
    common = sorted(set(cohort_expr.index.astype(str)) & set(labels_df["patient_id"].astype(str)))
    cohort_expr = cohort_expr.loc[common]
    y = labels_df.set_index("patient_id").loc[common]["pcr"].astype(int).to_numpy()
    print(f"Cohort: {len(common)} patients, {sum(y)} pCR", flush=True)

    universes = {
        "oncokb_all": oncokb_all,
        "oncokb_only": oncokb_only,
        "protein_coding": None,
    }
    grid = []
    for q in [0.20, 0.25, 0.30]:
        for top_n in [10, 20, 30, 50]:
            for asym_label, (up_n, down_n) in [("sym", (top_n, top_n)), ("up_heavy", (top_n, max(1, top_n // 3))), ("down_heavy", (max(1, top_n // 3), top_n))]:
                if top_n < 30 and asym_label != "sym":
                    continue  # asymmetry only meaningful for larger top-N
                for univ_name, univ_set in universes.items():
                    for fdr_label, fdr_gate in [("none", None), ("fdr_10", 0.10), ("fdr_05", 0.05)]:
                        if fdr_label != "none" and univ_name == "protein_coding":
                            continue  # too expensive
                        grid.append({
                            "q": q, "top_n_up": up_n, "top_n_down": down_n,
                            "asym": asym_label, "universe": univ_name,
                            "fdr_gate": fdr_label, "fdr_value": fdr_gate,
                        })
    print(f"Grid size: {len(grid)} configs", flush=True)

    rows = []
    for i, cfg in enumerate(grid):
        try:
            res = derive_quantile_sig(expr_solid, auc, cfg["q"], cfg["top_n_up"], cfg["top_n_down"],
                                       universes[cfg["universe"]], cfg["fdr_value"])
            if res is None:
                rows.append({**cfg, "n_up": 0, "n_down": 0, "auroc": np.nan, "status": "no_signature"})
                continue
            up_set, down_set = res
            scores = singscore(cohort_expr, up_set, down_set).loc[common].to_numpy()
            if pd.isna(scores).any():
                rows.append({**cfg, "n_up": len(up_set), "n_down": len(down_set), "auroc": np.nan, "status": "nan_scores"})
                continue
            if len(np.unique(y)) < 2:
                rows.append({**cfg, "n_up": len(up_set), "n_down": len(down_set), "auroc": np.nan, "status": "single_class"})
                continue
            auc_val = roc_auc_score(y, scores)
            rows.append({**cfg, "n_up": len(up_set), "n_down": len(down_set), "auroc": float(auc_val), "status": "ok"})
        except Exception as e:
            rows.append({**cfg, "n_up": 0, "n_down": 0, "auroc": np.nan, "status": f"error:{str(e)[:80]}"})
        if (i + 1) % 30 == 0:
            print(f"  {i+1}/{len(grid)} configs done", flush=True)

    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False, float_format="%.10g")

    # Top-quartile robustness check (pre-reg §8b)
    valid = df[df["status"] == "ok"].sort_values("auroc", ascending=False)
    if len(valid):
        # Pre-reg primary config: q=0.30, top_n_up=top_n_down=30, asym=sym, oncokb_all, fdr=none
        primary_row = valid[(valid["q"] == 0.30) & (valid["top_n_up"] == 30) &
                            (valid["top_n_down"] == 30) & (valid["asym"] == "sym") &
                            (valid["universe"] == "oncokb_all") & (valid["fdr_gate"] == "none")]
        if len(primary_row):
            primary_auc = primary_row["auroc"].values[0]
            rank = (valid["auroc"] > primary_auc).sum() + 1
            quartile = rank / len(valid)
            passed = quartile <= 0.25
        else:
            primary_auc = np.nan; rank = -1; quartile = np.nan; passed = False
    else:
        primary_auc = np.nan; rank = -1; quartile = np.nan; passed = False

    summary = {
        "cohort": args.cohort,
        "drug": args.drug,
        "n_grid_configs": len(grid),
        "n_valid_configs": int((df["status"] == "ok").sum()),
        "primary_config_auroc": float(primary_auc) if pd.notna(primary_auc) else None,
        "primary_config_rank": int(rank) if rank > 0 else None,
        "primary_config_quartile": float(quartile) if pd.notna(quartile) else None,
        "passes_top_quartile_gate": bool(passed),
        "criterion": "primary config must rank in top quartile (pre-reg §8b)",
    }
    with open(args.output.with_suffix(".json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
