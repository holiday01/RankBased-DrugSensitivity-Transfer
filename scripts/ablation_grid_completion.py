"""Phase 1 §4.5 — Ablation grid completion (36 cells).

Grid:
  q in {0.20, 0.25, 0.30}
  top_N in {10, 20, 30, 50}
  gene_universe in {OncoKB-all, OncoKB-only, protein-coding}
  = 3 x 4 x 3 = 36 cells

Pre-reg primary: q=0.30, top_N=30, OncoKB-all.

DISCREPANCY NOTE: The Phase 1 spec says "on TCGA (NOT on external)" but
constraints also say "DO NOT touch TCGA-BRCA". Since the pre-existing
run_ablation_grid.py is structured to evaluate AUROC on an external
cohort (CCLE-derived signature -> external cohort pCR), we follow that
pattern but restrict to GSE25066 (largest cohort, 488 patients) as the
robustness-grid evaluation target. The Phase 1 spec's "on TCGA" likely
refers to training-set (CCLE) signature derivation, NOT TCGA-BRCA
evaluation. We DO NOT touch TCGA-BRCA (Phase 3) here.

Outputs:
  results/phase_1/ablation_grid.tsv
  results/phase_1/ablation_grid_summary.json
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO_ROOT = Path(".")
REPO = REPO_ROOT / "."
OUT_DIR = REPO / "." / "results" / "phase_1"
CCLE_DRUG_PATH = REPO_ROOT / "enrichment_test" / "ccle_drug_common_with_tcga.csv"
CCLE_EXPR_PATH = REPO_ROOT / "CCLE25Q3" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_MODEL_PATH = REPO_ROOT / "CCLE25Q3" / "Model.csv"
ONCOKB_PATH = REPO_ROOT / "OncoKB" / "cancerGeneList.tsv"

TARGET_COHORT = "GSE25066"
TARGET_DRUG = "DOXORUBICIN"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clean_gene(g):
    return re.sub(r"\s*\(.*?\)", "", str(g)).upper().strip()


def norm_drug(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def load_solid():
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
    expr = expr.loc[~expr.index.isin(heme_ids)]

    sym_col = "Hugo Symbol" if "Hugo Symbol" in oncokb.columns else "Gene Symbol"
    oncogenes_all = set(oncokb[sym_col].astype(str).str.upper())
    if "Is Oncogene" in oncokb.columns:
        oncogenes_only = set(oncokb[oncokb["Is Oncogene"] == "Yes"][sym_col].astype(str).str.upper())
    else:
        oncogenes_only = oncogenes_all

    ccle_drug_map = {}
    for col in ccle_drug.columns[1:]:
        key = norm_drug(col.split("(")[0])
        if key not in ccle_drug_map:
            ccle_drug_map[key] = col
    if TARGET_DRUG not in ccle_drug_map:
        raise SystemExit(f"{TARGET_DRUG} not found in CCLE")
    auc = ccle_drug.set_index(ccle_drug.columns[0])[ccle_drug_map[TARGET_DRUG]].astype(float)
    return expr, auc, oncogenes_all, oncogenes_only


def derive(expr, auc, q, n_up, n_down, universe_set):
    common = expr.index.intersection(auc.index)
    auc_sub = auc.loc[common].dropna()
    if len(auc_sub) < 100:
        return None
    expr_sub = expr.loc[auc_sub.index]
    if universe_set is not None:
        cols = [c for c in expr_sub.columns if c in universe_set]
        if len(cols) < 50:
            return None
        expr_sub = expr_sub[cols]
    low = auc_sub.quantile(q); high = auc_sub.quantile(1 - q)
    sens = auc_sub[auc_sub <= low].index
    res = auc_sub[auc_sub >= high].index
    if len(sens) < 30 or len(res) < 30:
        return None
    diff = (expr_sub.loc[sens].mean(0) - expr_sub.loc[res].mean(0)).dropna()
    up = sorted(diff.nlargest(n_up).index.tolist())
    down = sorted(diff.nsmallest(n_down).index.tolist())
    return set(up), set(down)


def singscore(expr_cohort, up_set, down_set):
    if expr_cohort.shape[1] < 5000:
        raise ValueError("need full transcriptome")
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
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {"script": str(me), "script_sha256": sha256_of(me),
           "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
           "target_cohort": TARGET_COHORT, "target_drug": TARGET_DRUG,
           "grid_size": 36,
           "primary_config": {"q": 0.30, "top_n": 30, "universe": "oncokb_all"},
           "discrepancy_note": (
               "Spec says 'on TCGA' but constraints forbid touching TCGA-BRCA. "
               "Following pre-existing run_ablation_grid.py: signature derived from CCLE solid-only, "
               "evaluated against largest external cohort GSE25066 (488 patients). "
               "TCGA-BRCA Phase 3 deferred."),
           }

    print("Loading CCLE solid-only ...", flush=True)
    expr_ccle, auc, oncokb_all, oncokb_only = load_solid()

    cohort_expr = pd.read_csv(REPO / "external_data" / TARGET_COHORT / f"{TARGET_COHORT}_expression.tsv",
                              sep="\t", index_col=0)
    cohort_expr.columns = [c.upper() for c in cohort_expr.columns]
    labels = pd.read_csv(REPO / "external_data" / TARGET_COHORT / f"{TARGET_COHORT}_pcr.tsv", sep="\t")
    common = sorted(set(cohort_expr.index.astype(str)) & set(labels["patient_id"].astype(str)))
    cohort_expr = cohort_expr.loc[common]
    y = labels.set_index("patient_id").loc[common]["pcr"].astype(int).values

    universes = {
        "oncokb_all": oncokb_all,
        "oncokb_only": oncokb_only,
        "protein_coding": None,
    }

    rows = []
    for q in [0.20, 0.25, 0.30]:
        for top_n in [10, 20, 30, 50]:
            for uni_name, uni_set in universes.items():
                cfg = {"q": q, "top_n": top_n, "universe": uni_name}
                try:
                    res = derive(expr_ccle, auc, q, top_n, top_n, uni_set)
                    if res is None:
                        rows.append({**cfg, "n_up": 0, "n_down": 0, "auroc": np.nan, "status": "no_sig"})
                        continue
                    up, down = res
                    s = singscore(cohort_expr, up, down).loc[common].values
                    if np.isnan(s).any() or len(np.unique(y)) < 2:
                        rows.append({**cfg, "n_up": len(up), "n_down": len(down),
                                     "auroc": np.nan, "status": "skip"})
                        continue
                    auc_val = float(roc_auc_score(y, s))
                    rows.append({**cfg, "n_up": len(up), "n_down": len(down),
                                 "auroc": auc_val, "status": "ok"})
                except Exception as e:
                    rows.append({**cfg, "n_up": 0, "n_down": 0, "auroc": np.nan,
                                 "status": f"err:{str(e)[:80]}"})
                print(f"  q={q} N={top_n} u={uni_name}: AUROC={rows[-1]['auroc']}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "ablation_grid.tsv", sep="\t", index=False, float_format="%.6g")

    valid = df[df["status"] == "ok"].sort_values("auroc", ascending=False).reset_index(drop=True)
    primary = valid[(valid["q"] == 0.30) & (valid["top_n"] == 30) & (valid["universe"] == "oncokb_all")]
    if len(primary):
        primary_auc = float(primary["auroc"].values[0])
        rank = int((valid["auroc"] > primary_auc).sum() + 1)
        quartile = rank / len(valid)
    else:
        primary_auc = float("nan"); rank = -1; quartile = float("nan")

    summary = {
        "n_grid_configs": len(rows),
        "n_valid": int((df["status"] == "ok").sum()),
        "primary_config_AUROC": primary_auc,
        "primary_config_rank_within_valid": rank,
        "primary_config_quartile": quartile,
        "passes_top_quartile_check": (quartile <= 0.25 if np.isfinite(quartile) else False),
        "best_AUROC": float(valid["auroc"].max()) if len(valid) else float("nan"),
        "best_config": (valid.iloc[0].to_dict() if len(valid) else None),
        "median_AUROC": float(valid["auroc"].median()) if len(valid) else float("nan"),
    }
    log["summary"] = summary
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(OUT_DIR / "ablation_grid_summary.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    print(f"A5 done. Primary rank={rank}/{len(valid)} quartile={quartile} primary_AUROC={primary_auc:.3f}")


if __name__ == "__main__":
    main()
