"""Stage 0: Freeze CCLE-derived drug signatures (Round 1+2+3, parallelized, 2026-05-28).

Round 1 fixes (P0.1, P2.2, P2.3):
  - Spearman sign convention corrected: negative rho -> sensitivity (matches manuscript L93)
  - groupby(axis=1) replaced with pandas-2.x compatible idiom
  - Deterministic tie-break: secondary sort by gene_symbol asc
  - Pinned float serialization (%.10g) so SHA256s are pandas-version stable
  - direction = "sensitivity" / "resistance" (was "up" / "down")
  - MANIFEST adds: python/pandas/numpy/scipy versions, platform, script_sha256,
    skipped_drugs, drug-column-collision warnings

Round 2 (P2.1): bootstrap-Jaccard stability + delta-dropoff per signature
Round 3 (P1.1): solid-tumor-only variant + TOP2-poison consensus (Doxorubicin + Etoposide)

Parallelism: multiprocessing.Pool with N_WORKERS workers (default = nproc-2);
  parallelizes at (pool_label, strategy, drug_norm) task granularity (186 tasks).
"""

from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import os
import platform as _platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = Path("/home/holiday01/2026_ISMB_code")
CCLE_DRUG_PATH = REPO_ROOT / "enrichment_test" / "ccle_drug_common_with_tcga.csv"
CCLE_EXPR_PATH = REPO_ROOT / "CCLE25Q3" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
CCLE_MODEL_PATH = REPO_ROOT / "CCLE25Q3" / "Model.csv"
ONCOKB_PATH = REPO_ROOT / "OncoKB" / "cancerGeneList.tsv"
TCGA_DRUG_PATH = REPO_ROOT / "enrichment_test" / "tcga_drug_common_with_ccle.csv"

OUT_DIR = REPO_ROOT / "revise_bioadv" / "frozen_signatures"
SCRIPT_PATH = Path(__file__).resolve()

# ---------------------------------------------------------------------------
# Strategy config
# ---------------------------------------------------------------------------
STRATEGIES = [
    {"name": "quantile_30_oncokb_all", "auc_q": 0.3, "use_spearman": False, "gene_filter": "oncokb_all"},
    {"name": "quantile_20_protein",    "auc_q": 0.2, "use_spearman": False, "gene_filter": "protein_coding"},
    {"name": "spearman_top20_oncokb",  "auc_q": 0.5, "use_spearman": True,  "gene_filter": "oncokb_only"},
]
STRATEGY_BY_NAME = {s["name"]: s for s in STRATEGIES}

MIN_CELL_LINES = 100
MIN_SAMPLES_PER_GROUP = 50
MIN_VALID_GENES = 5
QUANTILE_TOP_N = 30
SPEARMAN_TOP_N = 20

N_BOOTSTRAP = 100
BOOTSTRAP_SEED = 42

TOP2_POISON_DRUGS = ["DOXORUBICIN", "ETOPOSIDE"]

N_WORKERS = max(1, os.cpu_count() - 2)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def clean_gene(g):
    return re.sub(r"\s*\(.*?\)", "", str(g)).upper().strip()


def norm_drug(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def deterministic_top_n_by_diff(diff: pd.Series, n: int, largest: bool) -> list:
    df = pd.DataFrame({"delta": diff.values, "gene": diff.index})
    df = df.sort_values(["delta", "gene"], ascending=[not largest, True], kind="mergesort")
    return df.head(n)["gene"].tolist()


def deterministic_top_n_by_rho(rho_dict: dict, n: int, sign: str) -> list:
    if sign == "pos":
        candidates = [(g, c) for g, c in rho_dict.items() if c > 0]
    else:
        candidates = [(g, c) for g, c in rho_dict.items() if c < 0]
    return sorted(candidates, key=lambda x: (-abs(x[1]), x[0]))[:n]


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


# ---------------------------------------------------------------------------
# Load & preprocess (called once in main)
# ---------------------------------------------------------------------------
def build_oncokb_sets(oncokb: pd.DataFrame) -> tuple[set, set]:
    if "Hugo Symbol" in oncokb.columns:
        sym_col = "Hugo Symbol"
    elif "Gene Symbol" in oncokb.columns:
        sym_col = "Gene Symbol"
    else:
        raise ValueError("OncoKB Hugo Symbol column not found")
    oncogenes_all = set(oncokb[sym_col].astype(str).str.upper())
    if "Is Oncogene" in oncokb.columns:
        oncogenes_only = set(oncokb[oncokb["Is Oncogene"] == "Yes"][sym_col].astype(str).str.upper())
    elif "Gene Type" in oncokb.columns:
        oncogenes_only = set(oncokb[oncokb["Gene Type"].str.upper().isin(["ONCOGENE", "ONCO"])][sym_col].astype(str).str.upper())
    else:
        oncogenes_only = oncogenes_all
    return oncogenes_all, oncogenes_only


def build_expr_base(ccle_drug: pd.DataFrame, ccle_exp: pd.DataFrame) -> pd.DataFrame:
    common_models = ccle_drug.iloc[:, 0].astype(str)
    expr = ccle_exp[ccle_exp.ModelID.isin(common_models)].copy()
    expr = expr.set_index("ModelID")
    meta_cols = ["Unnamed: 0", "SequencingID", "IsDefaultEntryForModel",
                 "ModelConditionID", "IsDefaultEntryForMC"]
    expr = expr.drop(columns=[c for c in meta_cols if c in expr.columns], errors="ignore")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr[~expr.index.duplicated(keep="first")]
    expr.columns = [clean_gene(g) for g in expr.columns]
    expr = expr.T.groupby(level=0).mean().T
    return expr


def identify_heme_model_ids(model_df: pd.DataFrame) -> set:
    heme_lineage_keys = ["myeloid", "lymphoid", "plasma cell"]
    heme_disease_keys = ["leukemia", "lymphoma", "myeloma", "myelodysplastic"]
    lineage = model_df["OncotreeLineage"].astype(str).str.lower().fillna("")
    disease = model_df["OncotreePrimaryDisease"].astype(str).str.lower().fillna("")
    lin_mask = lineage.apply(lambda s: any(k in s for k in heme_lineage_keys))
    dis_mask = disease.apply(lambda s: any(k in s for k in heme_disease_keys))
    return set(model_df.loc[lin_mask | dis_mask, "ModelID"].astype(str))


def build_ccle_drug_map_with_collisions(ccle_drug: pd.DataFrame) -> tuple[dict, list]:
    collisions = []
    seen = {}
    for col in ccle_drug.columns[1:]:
        base = col.split("(")[0]
        key = norm_drug(base)
        if key in seen:
            collisions.append({"key": key, "kept": seen[key], "dropped": col})
        else:
            seen[key] = col
    return seen, collisions


# ---------------------------------------------------------------------------
# Signature build
# ---------------------------------------------------------------------------
def build_signature(drug_norm: str, strategy: dict, ccle_drug: pd.DataFrame,
                    ccle_drug_map: dict, expr_pool: pd.DataFrame,
                    auc_override: pd.Series | None = None) -> dict | None:
    if auc_override is None:
        if drug_norm not in ccle_drug_map:
            return {"_skip_reason": "drug_not_in_ccle"}
        drug_col = ccle_drug_map[drug_norm]
        auc = ccle_drug.set_index(ccle_drug.columns[0])[drug_col].astype(float)
    else:
        auc = auc_override
        drug_col = None

    common = expr_pool.index.intersection(auc.index)
    if len(common) < MIN_CELL_LINES:
        return {"_skip_reason": f"too_few_cell_lines({len(common)}<{MIN_CELL_LINES})"}

    auc_sub = auc.loc[common].dropna()
    expr_sub = expr_pool.loc[auc_sub.index]
    if len(auc_sub) < MIN_CELL_LINES:
        return {"_skip_reason": f"too_few_after_dropna({len(auc_sub)})"}

    if strategy["use_spearman"]:
        rho_series = expr_sub.corrwith(auc_sub, method="spearman").dropna()
        if len(rho_series) < MIN_VALID_GENES:
            return {"_skip_reason": f"too_few_rho({len(rho_series)})"}
        rho_dict = rho_series.to_dict()
        # P0.1 FIX: negative rho -> sensitivity, positive rho -> resistance
        sens = deterministic_top_n_by_rho(rho_dict, SPEARMAN_TOP_N, sign="neg")
        res = deterministic_top_n_by_rho(rho_dict, SPEARMAN_TOP_N, sign="pos")
        up_sensitive = [(g, c) for g, c in sens]
        up_resistant = [(g, c) for g, c in res]
        n_sens = n_res = int(len(auc_sub))
        score_col = "rho"
    else:
        q = strategy["auc_q"]
        low = auc_sub.quantile(q)
        high = auc_sub.quantile(1 - q)
        sens_idx = auc_sub[auc_sub <= low].index
        res_idx = auc_sub[auc_sub >= high].index
        if len(sens_idx) < MIN_SAMPLES_PER_GROUP or len(res_idx) < MIN_SAMPLES_PER_GROUP:
            return {"_skip_reason": f"too_few_per_group({len(sens_idx)},{len(res_idx)})"}
        sens_expr = expr_sub.loc[sens_idx]
        res_expr = expr_sub.loc[res_idx]
        diff = (sens_expr.mean(axis=0) - res_expr.mean(axis=0)).dropna()
        if len(diff) < MIN_VALID_GENES * 2:
            return {"_skip_reason": f"too_few_diff({len(diff)})"}
        up_s_genes = deterministic_top_n_by_diff(diff, QUANTILE_TOP_N, largest=True)
        up_r_genes = deterministic_top_n_by_diff(diff, QUANTILE_TOP_N, largest=False)
        up_sensitive = [(g, float(diff[g])) for g in up_s_genes]
        up_resistant = [(g, float(diff[g])) for g in up_r_genes]
        n_sens = int(len(sens_idx))
        n_res = int(len(res_idx))
        score_col = "delta"

    if len(up_sensitive) < 5 or len(up_resistant) < 5:
        return {"_skip_reason": "fewer_than_5_genes"}

    top_s = abs(up_sensitive[0][1])
    last_s = abs(up_sensitive[-1][1])
    dropoff_s = last_s / top_s if top_s > 0 else float("nan")
    top_r = abs(up_resistant[0][1])
    last_r = abs(up_resistant[-1][1])
    dropoff_r = last_r / top_r if top_r > 0 else float("nan")

    return {
        "up_sensitive": [(g.upper(), s) for g, s in up_sensitive],
        "up_resistant": [(g.upper(), s) for g, s in up_resistant],
        "n_sensitive_cells": n_sens,
        "n_resistant_cells": n_res,
        "n_cell_lines": int(len(common)),
        "drug_col": drug_col,
        "score_col": score_col,
        "dropoff_ratio_sensitivity": float(dropoff_s),
        "dropoff_ratio_resistance": float(dropoff_r),
    }


def stability_jaccard(drug_norm: str, strategy: dict, ccle_drug: pd.DataFrame,
                       ccle_drug_map: dict, expr_pool: pd.DataFrame,
                       orig_sens: set, orig_res: set,
                       n_bootstrap: int = N_BOOTSTRAP, seed: int = BOOTSTRAP_SEED) -> dict:
    if drug_norm not in ccle_drug_map:
        return {k: float("nan") for k in ["jaccard_sens_median", "jaccard_sens_iqr",
                                          "jaccard_res_median", "jaccard_res_iqr",
                                          "n_bootstrap_successful"]}
    drug_col = ccle_drug_map[drug_norm]
    auc_full = ccle_drug.set_index(ccle_drug.columns[0])[drug_col].astype(float)
    common = expr_pool.index.intersection(auc_full.index)
    auc_full = auc_full.loc[common].dropna()
    if len(auc_full) < MIN_CELL_LINES:
        return {k: float("nan") for k in ["jaccard_sens_median", "jaccard_sens_iqr",
                                          "jaccard_res_median", "jaccard_res_iqr",
                                          "n_bootstrap_successful"]}
    rng = np.random.default_rng(seed)
    cell_ids = auc_full.index.tolist()
    js_sens, js_res = [], []
    for _ in range(n_bootstrap):
        sample = rng.choice(cell_ids, size=len(cell_ids), replace=True)
        auc_boot = auc_full.loc[sample].reset_index(drop=True)
        expr_boot = expr_pool.loc[sample].reset_index(drop=True)
        auc_boot.index = expr_boot.index
        sig = build_signature(drug_norm, strategy, None, None, expr_boot, auc_override=auc_boot)
        if sig is None or "_skip_reason" in sig:
            continue
        boot_s = {g for g, _ in sig["up_sensitive"]}
        boot_r = {g for g, _ in sig["up_resistant"]}
        js_sens.append(jaccard(boot_s, orig_sens))
        js_res.append(jaccard(boot_r, orig_res))
    return {
        "jaccard_sens_median": float(np.median(js_sens)) if js_sens else float("nan"),
        "jaccard_sens_iqr": float(np.percentile(js_sens, 75) - np.percentile(js_sens, 25)) if js_sens else float("nan"),
        "jaccard_res_median": float(np.median(js_res)) if js_res else float("nan"),
        "jaccard_res_iqr": float(np.percentile(js_res, 75) - np.percentile(js_res, 25)) if js_res else float("nan"),
        "n_bootstrap_successful": int(len(js_sens)),
    }


def build_top2_poison_consensus(per_drug_sigs: dict, drugs: list = TOP2_POISON_DRUGS) -> dict | None:
    sigs = [per_drug_sigs.get(d) for d in drugs]
    sigs = [s for s in sigs if s is not None and "_skip_reason" not in s]
    if len(sigs) < 2:
        return None

    def rank_map(pairs):
        return {g.upper(): r for r, (g, _) in enumerate(pairs, start=1)}

    sens_ranks = [rank_map(s["up_sensitive"]) for s in sigs]
    res_ranks = [rank_map(s["up_resistant"]) for s in sigs]

    sens_genes, res_genes = {}, {}
    for rmap in sens_ranks:
        for g, r in rmap.items():
            sens_genes.setdefault(g, []).append(r)
    for rmap in res_ranks:
        for g, r in rmap.items():
            res_genes.setdefault(g, []).append(r)

    sens_combined = sorted(
        [(g, sum(rs) / len(rs)) for g, rs in sens_genes.items() if len(rs) == len(sigs)],
        key=lambda x: (x[1], x[0]),
    )
    res_combined = sorted(
        [(g, sum(rs) / len(rs)) for g, rs in res_genes.items() if len(rs) == len(sigs)],
        key=lambda x: (x[1], x[0]),
    )
    return {
        "up_sensitive": sens_combined[:QUANTILE_TOP_N],
        "up_resistant": res_combined[:QUANTILE_TOP_N],
        "n_input_drugs": len(sigs),
        "input_drugs": list(drugs[:len(sigs)]),
        "score_col": "mean_rank",
    }


def stability_consensus_jaccard(
    pool_label: str,
    strategy: dict,
    expr_pool: pd.DataFrame,
    orig_sens: set,
    orig_res: set,
    drugs: list = TOP2_POISON_DRUGS,
    n_bootstrap: int = N_BOOTSTRAP,
    seed: int = BOOTSTRAP_SEED,
) -> dict:
    """Bootstrap stability for the rank-averaged consensus.

    Restrict to cells with AUC for ALL input drugs (so the resample is
    well-defined for every drug), then for each resample rebuild each drug
    signature and the consensus, and compute Jaccard with the frozen original.
    """
    # Pre-compute per-drug AUC restricted to the expression pool
    auc_by_drug = {}
    for drug_norm in drugs:
        if drug_norm not in _CCLE_DRUG_MAP:
            continue
        drug_col = _CCLE_DRUG_MAP[drug_norm]
        auc = _CCLE_DRUG.set_index(_CCLE_DRUG.columns[0])[drug_col].astype(float)
        auc = auc.loc[expr_pool.index.intersection(auc.index)].dropna()
        auc_by_drug[drug_norm] = auc

    if len(auc_by_drug) < 2:
        return {k: float("nan") for k in ["jaccard_sens_median", "jaccard_sens_iqr",
                                          "jaccard_res_median", "jaccard_res_iqr",
                                          "n_bootstrap_successful"]}

    # Bootstrap pool = intersection of all drugs' AUC-available cells
    common_cells = set.intersection(*[set(s.index) for s in auc_by_drug.values()])
    common_cells = sorted(common_cells & set(expr_pool.index))
    if len(common_cells) < MIN_CELL_LINES:
        return {k: float("nan") for k in ["jaccard_sens_median", "jaccard_sens_iqr",
                                          "jaccard_res_median", "jaccard_res_iqr",
                                          "n_bootstrap_successful"]}

    expr_common = expr_pool.loc[common_cells]
    auc_common = {d: s.loc[common_cells] for d, s in auc_by_drug.items()}

    rng = np.random.default_rng(seed)
    js_sens, js_res = [], []
    for _ in range(n_bootstrap):
        sample = rng.choice(common_cells, size=len(common_cells), replace=True)
        expr_boot = expr_common.loc[sample].reset_index(drop=True)
        boot_index = expr_boot.index
        boot_sigs = {}
        for drug_norm, auc_full in auc_common.items():
            auc_boot = auc_full.loc[sample].reset_index(drop=True)
            auc_boot.index = boot_index
            sig = build_signature(drug_norm, strategy, None, None, expr_boot, auc_override=auc_boot)
            if sig is None or "_skip_reason" in sig:
                continue
            boot_sigs[drug_norm] = sig
        if len(boot_sigs) < 2:
            continue
        consensus_boot = build_top2_poison_consensus(boot_sigs, drugs)
        if consensus_boot is None:
            continue
        boot_s = {g for g, _ in consensus_boot["up_sensitive"]}
        boot_r = {g for g, _ in consensus_boot["up_resistant"]}
        js_sens.append(jaccard(boot_s, orig_sens))
        js_res.append(jaccard(boot_r, orig_res))
    return {
        "jaccard_sens_median": float(np.median(js_sens)) if js_sens else float("nan"),
        "jaccard_sens_iqr": float(np.percentile(js_sens, 75) - np.percentile(js_sens, 25)) if js_sens else float("nan"),
        "jaccard_res_median": float(np.median(js_res)) if js_res else float("nan"),
        "jaccard_res_iqr": float(np.percentile(js_res, 75) - np.percentile(js_res, 25)) if js_res else float("nan"),
        "n_bootstrap_successful": int(len(js_sens)),
    }


# ---------------------------------------------------------------------------
# TSV dump
# ---------------------------------------------------------------------------
def dump_signature_tsv(out_path: Path, sig: dict):
    score_col = sig["score_col"]
    rows = []
    for rank, (gene, score) in enumerate(sig["up_sensitive"], start=1):
        rows.append({"gene_symbol": gene, "direction": "sensitivity",
                     "rank_in_direction": rank, score_col: score})
    for rank, (gene, score) in enumerate(sig["up_resistant"], start=1):
        rows.append({"gene_symbol": gene, "direction": "resistance",
                     "rank_in_direction": rank, score_col: score})
    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False, float_format="%.10g", lineterminator="\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Worker globals (set in main, inherited by children via fork)
# ---------------------------------------------------------------------------
_CCLE_DRUG = None
_CCLE_DRUG_MAP = None
_EXPR_POOLS = None  # dict: {pool_label: {gene_filter: DataFrame}}


def worker_task(task):
    pool_label, strategy_name, drug_norm = task
    strategy = STRATEGY_BY_NAME[strategy_name]
    expr_pool = _EXPR_POOLS[pool_label][strategy["gene_filter"]]
    sig = build_signature(drug_norm, strategy, _CCLE_DRUG, _CCLE_DRUG_MAP, expr_pool)
    out = {"task": task, "sig": sig}
    if sig is None or "_skip_reason" in sig:
        return out
    out["stability"] = stability_jaccard(
        drug_norm, strategy, _CCLE_DRUG, _CCLE_DRUG_MAP, expr_pool,
        orig_sens={g for g, _ in sig["up_sensitive"]},
        orig_res={g for g, _ in sig["up_resistant"]},
    )
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    global _CCLE_DRUG, _CCLE_DRUG_MAP, _EXPR_POOLS

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    t0 = datetime.now()
    print(f"[{datetime.now()}] Loading data ...", flush=True)
    ccle_drug = pd.read_csv(CCLE_DRUG_PATH)
    ccle_exp = pd.read_csv(CCLE_EXPR_PATH)
    oncokb = pd.read_csv(ONCOKB_PATH, sep="\t")
    tcga_drug = pd.read_csv(TCGA_DRUG_PATH)
    model_df = pd.read_csv(CCLE_MODEL_PATH, low_memory=False)

    oncogenes_all, oncogenes_only = build_oncokb_sets(oncokb)
    print(f"[{datetime.now()}] Building expression matrices ...", flush=True)
    expr_base_pan = build_expr_base(ccle_drug, ccle_exp)
    heme_ids = identify_heme_model_ids(model_df)
    expr_base_solid = expr_base_pan.loc[~expr_base_pan.index.isin(heme_ids)]

    expr_pools = {
        "pan_cancer": {
            "oncokb_all": expr_base_pan[[g for g in expr_base_pan.columns if g in oncogenes_all]],
            "oncokb_only": expr_base_pan[[g for g in expr_base_pan.columns if g in oncogenes_only]],
            "protein_coding": expr_base_pan,
        },
        "solid_only": {
            "oncokb_all": expr_base_solid[[g for g in expr_base_solid.columns if g in oncogenes_all]],
            "oncokb_only": expr_base_solid[[g for g in expr_base_solid.columns if g in oncogenes_only]],
            "protein_coding": expr_base_solid,
        },
    }
    ccle_drug_map, collisions = build_ccle_drug_map_with_collisions(ccle_drug)
    if "drug_name" in tcga_drug.columns:
        tcga_drugs = tcga_drug["drug_name"].dropna().astype(str).unique()
    else:
        tcga_drugs = tcga_drug[tcga_drug.columns[0]].dropna().astype(str).unique()
    tcga_norm_map = {norm_drug(d): d for d in tcga_drugs}
    common_drugs_norm = sorted(set(tcga_norm_map.keys()) & set(ccle_drug_map.keys()))

    print(f"  pan-cancer cells: {expr_base_pan.shape[0]}, solid-only cells: {expr_base_solid.shape[0]}", flush=True)
    print(f"  common drugs: {len(common_drugs_norm)}", flush=True)
    print(f"  gene pool sizes: oncokb_all={expr_pools['pan_cancer']['oncokb_all'].shape[1]}, "
          f"oncokb_only={expr_pools['pan_cancer']['oncokb_only'].shape[1]}, "
          f"protein_coding={expr_pools['pan_cancer']['protein_coding'].shape[1]}", flush=True)

    _CCLE_DRUG = ccle_drug
    _CCLE_DRUG_MAP = ccle_drug_map
    _EXPR_POOLS = expr_pools

    # Build task list
    tasks = []
    for pool_label in ["pan_cancer", "solid_only"]:
        for strategy in STRATEGIES:
            for drug_norm in common_drugs_norm:
                tasks.append((pool_label, strategy["name"], drug_norm))
    print(f"\n[{datetime.now()}] Submitting {len(tasks)} tasks to {N_WORKERS} workers ...", flush=True)

    # Parallel execution. Use fork (default on Linux); children inherit _EXPR_POOLS etc.
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=N_WORKERS) as pool:
        results = []
        for i, r in enumerate(pool.imap_unordered(worker_task, tasks, chunksize=2)):
            results.append(r)
            if (i + 1) % 30 == 0:
                print(f"  [{datetime.now()}] {i + 1}/{len(tasks)} done", flush=True)
    print(f"[{datetime.now()}] All {len(tasks)} tasks complete.", flush=True)

    # Group + write
    pool_manifests = {"pan_cancer": {"outputs": {}, "skipped_drugs": []},
                      "solid_only": {"outputs": {}, "skipped_drugs": []}}
    per_drug_sigs_by_pool_strat = {p: {s["name"]: {} for s in STRATEGIES} for p in pool_manifests}

    for r in results:
        pool_label, strategy_name, drug_norm = r["task"]
        sig = r["sig"]
        if sig is None or "_skip_reason" in sig:
            pool_manifests[pool_label]["skipped_drugs"].append({
                "drug_norm": drug_norm, "strategy": strategy_name,
                "reason": sig["_skip_reason"] if sig else "build_returned_none",
            })
            continue
        per_drug_sigs_by_pool_strat[pool_label][strategy_name][drug_norm] = sig

        # Write TSV
        out_path = (OUT_DIR if pool_label == "pan_cancer" else OUT_DIR / "solid_only") / strategy_name / f"{drug_norm}.tsv"
        dump_signature_tsv(out_path, sig)
        rel = str(out_path.relative_to(REPO_ROOT))
        pool_manifests[pool_label]["outputs"].setdefault(strategy_name, []).append({
            "drug_norm": drug_norm,
            "drug_ccle_col": sig["drug_col"],
            "n_up": len(sig["up_sensitive"]),
            "n_down": len(sig["up_resistant"]),
            "n_sensitive_cells": sig["n_sensitive_cells"],
            "n_resistant_cells": sig["n_resistant_cells"],
            "n_cell_lines_used": sig["n_cell_lines"],
            "dropoff_ratio_sensitivity": sig["dropoff_ratio_sensitivity"],
            "dropoff_ratio_resistance": sig["dropoff_ratio_resistance"],
            "stability": r.get("stability"),
            "file": rel,
            "sha256": sha256_of(out_path),
        })

    # Consensus per (pool, strategy) + bootstrap stability for each
    print(f"[{datetime.now()}] Computing consensus + consensus stability ...", flush=True)
    consensus_tasks = []
    consensus_baselines = {}
    for pool_label, by_strat in per_drug_sigs_by_pool_strat.items():
        for strat_name, drug_to_sig in by_strat.items():
            if all(d in drug_to_sig for d in TOP2_POISON_DRUGS):
                consensus = build_top2_poison_consensus(drug_to_sig)
                if consensus is not None:
                    consensus_baselines[(pool_label, strat_name)] = consensus
                    consensus_tasks.append((pool_label, strat_name))

    # Compute consensus stability in parallel
    def _consensus_stability_task(task):
        pool_label, strat_name = task
        strategy = STRATEGY_BY_NAME[strat_name]
        expr_pool = _EXPR_POOLS[pool_label][strategy["gene_filter"]]
        consensus = consensus_baselines[(pool_label, strat_name)]
        stab = stability_consensus_jaccard(
            pool_label, strategy, expr_pool,
            orig_sens={g for g, _ in consensus["up_sensitive"]},
            orig_res={g for g, _ in consensus["up_resistant"]},
        )
        return (pool_label, strat_name, stab)

    # Inline serial (consensus tasks are only 6, each ~30s; not worth pickling consensus_baselines)
    consensus_stability_results = {}
    for task in consensus_tasks:
        pool_label, strat_name, stab = _consensus_stability_task(task)
        consensus_stability_results[(pool_label, strat_name)] = stab
        print(f"  [{pool_label}/{strat_name}] consensus stability: sens={stab['jaccard_sens_median']:.3f}, "
              f"res={stab['jaccard_res_median']:.3f}", flush=True)

    # Write consensus TSVs + manifest entries
    for (pool_label, strat_name), consensus in consensus_baselines.items():
        out_path = (OUT_DIR if pool_label == "pan_cancer" else OUT_DIR / "solid_only") / strat_name / "TOP2_POISON_CONSENSUS.tsv"
        dump_signature_tsv(out_path, consensus)
        pool_manifests[pool_label]["outputs"][strat_name].append({
            "drug_norm": "TOP2_POISON_CONSENSUS",
            "drug_ccle_col": None,
            "n_up": len(consensus["up_sensitive"]),
            "n_down": len(consensus["up_resistant"]),
            "n_input_drugs": consensus["n_input_drugs"],
            "input_drugs": consensus["input_drugs"],
            "score_col": consensus["score_col"],
            "stability": consensus_stability_results.get((pool_label, strat_name)),
            "file": str(out_path.relative_to(REPO_ROOT)),
            "sha256": sha256_of(out_path),
        })

    # Sort each strategy's drugs alphabetically for deterministic manifest
    for pm in pool_manifests.values():
        for k in pm["outputs"]:
            pm["outputs"][k].sort(key=lambda x: x["drug_norm"])
        pm["skipped_drugs"].sort(key=lambda x: (x["strategy"], x["drug_norm"]))

    # Capture BLAS / threading env
    import io
    np_config_buf = io.StringIO()
    try:
        np.show_config(mode="dicts") if "mode" in np.show_config.__code__.co_varnames else np.show_config()
    except Exception:
        pass
    blas_env = {k: os.environ.get(k, "(unset)") for k in
                ["OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"]}

    # Auxiliary artifacts: every script + every planning .md
    AUX_DIR = REPO_ROOT / "revise_bioadv"
    aux_scripts = sorted((AUX_DIR / "scripts").glob("*.py"))
    aux_mds = sorted([p for p in AUX_DIR.glob("*.md")])
    auxiliary = []
    for p in aux_scripts + aux_mds:
        auxiliary.append({
            "path": str(p.relative_to(REPO_ROOT)),
            "sha256": sha256_of(p),
            "size_bytes": p.stat().st_size,
        })

    manifest = {
        "stage": "0_freeze_signatures_round123_parallel",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "wall_time_seconds": (datetime.now() - t0).total_seconds(),
        "software": {
            "python_version": sys.version,
            "platform": _platform.platform(),
            "n_workers": N_WORKERS,
            "pandas": pd.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "blas_env": blas_env,
            "mp_start_method": mp.get_start_method(allow_none=True) or "fork",
            "script_path": str(SCRIPT_PATH),
            "script_sha256": sha256_of(SCRIPT_PATH),
        },
        "auxiliary_artifacts": auxiliary,
        "source_files": {
            "ccle_drug": {"path": str(CCLE_DRUG_PATH), "sha256": sha256_of(CCLE_DRUG_PATH)},
            "ccle_expression": {"path": str(CCLE_EXPR_PATH), "sha256": sha256_of(CCLE_EXPR_PATH)},
            "ccle_model": {"path": str(CCLE_MODEL_PATH), "sha256": sha256_of(CCLE_MODEL_PATH)},
            "oncokb": {"path": str(ONCOKB_PATH), "sha256": sha256_of(ONCOKB_PATH)},
            "tcga_drug": {"path": str(TCGA_DRUG_PATH), "sha256": sha256_of(TCGA_DRUG_PATH)},
        },
        "params": {
            "strategies": STRATEGIES,
            "MIN_CELL_LINES": MIN_CELL_LINES, "MIN_SAMPLES_PER_GROUP": MIN_SAMPLES_PER_GROUP,
            "MIN_VALID_GENES": MIN_VALID_GENES, "QUANTILE_TOP_N": QUANTILE_TOP_N,
            "SPEARMAN_TOP_N": SPEARMAN_TOP_N, "N_BOOTSTRAP": N_BOOTSTRAP,
            "BOOTSTRAP_SEED": BOOTSTRAP_SEED, "TOP2_POISON_DRUGS": TOP2_POISON_DRUGS,
            "round_1_fixes": ["spearman_sign_corrected_neg_rho_to_sensitivity",
                              "deterministic_tie_break_by_gene_symbol",
                              "float_format_pinned_%.10g",
                              "direction_labels_sensitivity_resistance",
                              "pandas2x_compatible_groupby"],
            "round_2": ["bootstrap_jaccard_stability_100x", "delta_dropoff_ratio"],
            "round_3": ["solid_only_pool_excluding_heme", "TOP2_poison_consensus(Doxorubicin+Etoposide)"],
            "parallel_backend": "multiprocessing.Pool(fork)",
        },
        "drug_column_collisions": collisions,
        "drugs": [{"drug_norm": d, "tcga_name": tcga_norm_map.get(d)} for d in common_drugs_norm],
        "pan_cancer": pool_manifests["pan_cancer"],
        "solid_only": pool_manifests["solid_only"],
    }
    with open(OUT_DIR / "MANIFEST.json", "w") as f:
        json.dump(manifest, f, indent=2)

    n_pan = sum(len(v) for v in pool_manifests["pan_cancer"]["outputs"].values())
    n_solid = sum(len(v) for v in pool_manifests["solid_only"]["outputs"].values())
    print(f"\n[{datetime.now()}] DONE in {(datetime.now()-t0).total_seconds():.1f}s.")
    print(f"  Pan-cancer artifacts: {n_pan}; Solid-only artifacts: {n_solid}")
    print(f"  Pan-cancer skipped: {len(pool_manifests['pan_cancer']['skipped_drugs'])}; "
          f"Solid-only skipped: {len(pool_manifests['solid_only']['skipped_drugs'])}")
    print(f"  MANIFEST: {OUT_DIR / 'MANIFEST.json'}")


if __name__ == "__main__":
    main()
