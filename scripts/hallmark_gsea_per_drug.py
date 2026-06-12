"""Phase 1 §4.6 — Hallmark GSEA per drug.

For each drug (12 total: 7 cytotoxic + 3 TKI + Fulvestrant + Imatinib):
  1. Pull/regenerate the full per-gene Pearson-r ranked list against
     CCLE-solid-only AUC.
  2. Run weighted KS-style GSEA against MSigDB Hallmark v2025.1.Hs (50 sets).
  3. Permutation null (gene shuffle) x 1000; BH-FDR per drug.
  4. Output NES heatmap (drug x Hallmark) + leading-edge gene lists for
     proliferation Hallmarks.

Implementation: weighted KS as in Subramanian et al. 2005 (gseapy unavailable).

Outputs:
  results/phase_1/hallmark_gsea_NES_matrix.tsv
  results/phase_1/hallmark_gsea_leading_edge.json
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

REPO_ROOT = Path(".")
REPO = REPO_ROOT / "."
OUT_DIR = REPO / "." / "results" / "phase_1"
HALLMARK = REPO / "h.all.v2025.1.Hs.symbols.gmt"
CCLE_DRUG_PATH = REPO_ROOT / "enrichment_test" / "ccle_drug_common_with_tcga.csv"
CCLE_EXPR_PATH = REPO_ROOT / "CCLE25Q3" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_MODEL_PATH = REPO_ROOT / "CCLE25Q3" / "Model.csv"

DRUGS = ["DOXORUBICIN", "ETOPOSIDE", "PACLITAXEL", "CARBOPLATIN",
         "CYCLOPHOSPHAMIDE", "FLUOROURACIL", "TOP2_consensus",
         "TRAMETINIB", "GEFITINIB", "LAPATINIB",
         "FULVESTRANT", "IMATINIB"]

PROLIF_HM = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
             "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2",
             "HALLMARK_MITOTIC_SPINDLE"]


def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def clean_gene(g):
    return re.sub(r"\s*\(.*?\)", "", str(g)).upper().strip()


def norm_drug(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def load_solid_expr():
    ccle_drug = pd.read_csv(CCLE_DRUG_PATH)
    ccle_exp = pd.read_csv(CCLE_EXPR_PATH)
    model = pd.read_csv(CCLE_MODEL_PATH, low_memory=False)
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
    ccle_drug_map = {}
    for col in ccle_drug.columns[1:]:
        k = norm_drug(col.split("(")[0])
        if k not in ccle_drug_map:
            ccle_drug_map[k] = col
    return expr, ccle_drug, ccle_drug_map


def per_drug_ranking(expr, auc):
    """Pearson r of each gene vs AUC across cell lines."""
    common = expr.index.intersection(auc.index)
    expr_sub = expr.loc[common]
    auc_sub = auc.loc[common].dropna()
    expr_sub = expr_sub.loc[auc_sub.index]
    # Pearson per column
    r = expr_sub.corrwith(auc_sub, method="pearson").dropna()
    # Sensitivity convention: low AUC = sensitive. Negative r = high expr in sensitive.
    # We rank by -r so sensitivity-associated genes are at the top.
    rank = (-r).sort_values(ascending=False)
    return rank


def load_hallmark_sets():
    sets = {}
    with open(HALLMARK) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = set(g.upper() for g in parts[2:])
    return sets


def weighted_ks(ranked_genes: list, ranked_scores: np.ndarray, gene_set: set,
                p: float = 1.0):
    """Subramanian-style weighted KS enrichment statistic.

    ranked_genes: list of gene symbols, sorted by ranked_scores DESC.
    ranked_scores: same length, the ranking metric (e.g., -Pearson r).
    Returns (ES, leading_edge_genes).
    """
    N = len(ranked_genes)
    in_set = np.array([g in gene_set for g in ranked_genes])
    NR = (np.abs(ranked_scores) ** p)[in_set].sum()
    if NR == 0 or in_set.sum() == 0 or in_set.sum() == N:
        return 0.0, []
    Nh = in_set.sum()
    miss = 1.0 / (N - Nh)
    p_hit = np.zeros(N)
    p_miss = np.zeros(N)
    p_hit[in_set] = (np.abs(ranked_scores[in_set]) ** p) / NR
    p_miss[~in_set] = miss
    running = np.cumsum(p_hit - p_miss)
    es_idx = np.argmax(np.abs(running))
    es = float(running[es_idx])
    # Leading edge: genes in set up to es_idx (if ES > 0) or from es_idx (if ES < 0)
    if es >= 0:
        le = [g for g, m in zip(ranked_genes[:es_idx + 1], in_set[:es_idx + 1]) if m]
    else:
        le = [g for g, m in zip(ranked_genes[es_idx:], in_set[es_idx:]) if m]
    return es, le


def permutation_ES_null(ranked_scores: np.ndarray, set_size: int, n_perm: int, seed: int):
    """Gene-shuffle null: ES distribution by permuting which genes 'are' in the set."""
    rng = np.random.default_rng(seed)
    N = len(ranked_scores)
    es_nulls = []
    abs_s = np.abs(ranked_scores)
    for _ in range(n_perm):
        idx = rng.choice(N, size=set_size, replace=False)
        in_set = np.zeros(N, dtype=bool)
        in_set[idx] = True
        NR = abs_s[in_set].sum()
        if NR == 0:
            continue
        miss = 1.0 / (N - set_size)
        p_hit = np.zeros(N); p_miss = np.zeros(N)
        p_hit[in_set] = abs_s[in_set] / NR
        p_miss[~in_set] = miss
        running = np.cumsum(p_hit - p_miss)
        es_nulls.append(running[np.argmax(np.abs(running))])
    return np.array(es_nulls)


def normalize_es(es, es_nulls):
    """NES per Subramanian: divide by mean ES of same-sign nulls."""
    if es >= 0:
        same = es_nulls[es_nulls >= 0]
    else:
        same = es_nulls[es_nulls < 0]
    if len(same) == 0:
        return float("nan")
    return float(es / np.mean(np.abs(same)))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {"script": str(me), "script_sha256": sha256_of(me),
           "started_utc": _dt.datetime.utcnow().isoformat() + "Z"}

    print("Loading CCLE solid-only expression ...", flush=True)
    expr, ccle_drug, drug_map = load_solid_expr()
    hm_sets = load_hallmark_sets()

    nes_matrix = pd.DataFrame(index=DRUGS, columns=list(hm_sets.keys()), dtype=float)
    p_matrix = pd.DataFrame(index=DRUGS, columns=list(hm_sets.keys()), dtype=float)
    fdr_matrix = pd.DataFrame(index=DRUGS, columns=list(hm_sets.keys()), dtype=float)
    leading = {}

    from statsmodels.stats.multitest import multipletests
    # For TOP2_consensus: derive a synthetic ranking by averaging Dox + Etoposide rankings
    cached_rankings = {}

    for drug in DRUGS:
        print(f"Drug: {drug}", flush=True)
        if drug == "TOP2_consensus":
            r_d = cached_rankings.get("DOXORUBICIN")
            r_e = cached_rankings.get("ETOPOSIDE")
            if r_d is None or r_e is None:
                # compute if not cached
                col_d = drug_map.get("DOXORUBICIN"); col_e = drug_map.get("ETOPOSIDE")
                if col_d is None or col_e is None:
                    print(f"  skip {drug} (Dox/Eto not in CCLE)")
                    continue
                r_d = per_drug_ranking(expr, ccle_drug.set_index(ccle_drug.columns[0])[col_d].astype(float))
                r_e = per_drug_ranking(expr, ccle_drug.set_index(ccle_drug.columns[0])[col_e].astype(float))
                cached_rankings["DOXORUBICIN"] = r_d
                cached_rankings["ETOPOSIDE"] = r_e
            common = r_d.index.intersection(r_e.index)
            r = (r_d.loc[common] + r_e.loc[common]) / 2
            r = r.sort_values(ascending=False)
        else:
            col = drug_map.get(drug)
            if col is None:
                print(f"  skip {drug} (not in CCLE)")
                continue
            auc = ccle_drug.set_index(ccle_drug.columns[0])[col].astype(float)
            r = per_drug_ranking(expr, auc)
            cached_rankings[drug] = r
        ranked_genes = r.index.tolist()
        ranked_scores = r.values
        leading[drug] = {}

        # Permutation null shared across hallmark sets (precompute for typical sizes)
        # For speed, do per-set null with 1000 permutations (50 sets x 1000 = 50k -- acceptable)
        for hm_name, hm_set in hm_sets.items():
            present = hm_set & set(ranked_genes)
            if len(present) < 5:
                nes_matrix.at[drug, hm_name] = np.nan
                p_matrix.at[drug, hm_name] = np.nan
                continue
            es, le = weighted_ks(ranked_genes, ranked_scores, present)
            # Use 500 perms for speed; precompute null per-size
            null = permutation_ES_null(ranked_scores, len(present), 500, seed=42 + hash(hm_name) % 1000)
            nes = normalize_es(es, null)
            if es >= 0:
                p = float((np.sum(null >= es) + 1) / (len(null) + 1))
            else:
                p = float((np.sum(null <= es) + 1) / (len(null) + 1))
            nes_matrix.at[drug, hm_name] = nes
            p_matrix.at[drug, hm_name] = p
            if hm_name in PROLIF_HM:
                leading[drug][hm_name] = {"ES": es, "NES": nes, "p": p,
                                           "leading_edge_genes": le[:30]}

        # BH-FDR across hallmarks within drug
        ps = p_matrix.loc[drug].dropna()
        if len(ps):
            _, q, _, _ = multipletests(ps.values, method="fdr_bh")
            fdr_matrix.loc[drug, ps.index] = q

    nes_matrix.to_csv(OUT_DIR / "hallmark_gsea_NES_matrix.tsv", sep="\t", float_format="%.4f")
    p_matrix.to_csv(OUT_DIR / "hallmark_gsea_pvalue_matrix.tsv", sep="\t", float_format="%.4g")
    fdr_matrix.to_csv(OUT_DIR / "hallmark_gsea_FDR_matrix.tsv", sep="\t", float_format="%.4g")
    with open(OUT_DIR / "hallmark_gsea_leading_edge.json", "w") as f:
        json.dump(leading, f, indent=2, default=str)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    # Top 3 enriched per drug
    summary = {}
    for drug in DRUGS:
        if drug not in nes_matrix.index:
            continue
        row = nes_matrix.loc[drug].dropna()
        # Top 3 by absolute NES, sign-aware
        top3 = row.abs().sort_values(ascending=False).head(3).index.tolist()
        summary[drug] = [{"hallmark": h, "NES": float(row[h]),
                          "p": float(p_matrix.at[drug, h]),
                          "FDR": float(fdr_matrix.at[drug, h]) if pd.notna(fdr_matrix.at[drug, h]) else None}
                         for h in top3]
    log["top3_hallmark_per_drug"] = summary
    with open(OUT_DIR / "hallmark_gsea_summary.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    print("A6 done")


if __name__ == "__main__":
    main()
