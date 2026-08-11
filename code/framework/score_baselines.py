"""B4 trivial-baseline + negative-control scoring on external cohorts.

Computes:
  - Trivial baselines (A4 + R1-M4):
      * MKI67 single-gene score (positive-control proliferation marker)
      * Hallmark E2F_TARGETS singscore
      * Hallmark G2M_CHECKPOINT singscore
      * Random gene-set baseline (size-matched, seed-pinned)
  - Negative-control signatures (A4):
      * DASATINIB / GEFITINIB / ERLOTINIB (solid_only/q30 frozen TSVs)
      * Expected AUROC ≈ 0.5 in GSE16446 (wrong MoA for anthracycline cohort)
  - Compares against the primary framework signature (DOXORUBICIN solid_only/q30)
      via paired DeLong test for the R1-M4 success criterion (Δ ≥ +0.03,
      paired p < 0.0083 in ≥ 2 of 3 cohorts).

Usage:
  python score_baselines.py \
      --cohort GSE16446 \
      --expression external_data/GSE16446/GSE16446_expression.tsv \
      --labels external_data/GSE16446/GSE16446_pcr.tsv \
      --hallmark h.all.v2025.1.Hs.symbols.gmt \
      --output external_data/GSE16446/baseline_results.tsv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score


REPO_ROOT = Path("/home/holiday01/2026_ISMB_code")
HALLMARK_PATH = REPO_ROOT / "enrichment_test" / "h.all.v2025.1.Hs.symbols.gmt"
FROZEN_DIR = REPO_ROOT / "revise_bioadv" / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all"

NEGATIVE_CONTROL_DRUGS = ["DASATINIB", "GEFITINIB", "ERLOTINIB"]
PRIMARY_DRUG = "DOXORUBICIN"

PROLIFERATION_HALLMARKS = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def read_gmt(path: Path, pathways: list = None) -> dict:
    sets = {}
    with open(path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            if pathways is not None and name not in pathways:
                continue
            sets[name] = set(g.upper() for g in parts[2:] if g.strip())
    return sets


def load_signature_tsv(path: Path) -> tuple[set, set]:
    df = pd.read_csv(path, sep="\t")
    sens = set(df[df["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(df[df["direction"] == "resistance"]["gene_symbol"].str.upper())
    return sens, res


def compute_singscore_bidirectional(expr: pd.DataFrame, up: set, down: set) -> pd.Series:
    """Full-transcriptome rank-based singscore with up + down."""
    if expr.shape[1] < 5000:
        raise ValueError(f"Expression matrix has only {expr.shape[1]} genes — refusing to score (full-transcriptome required).")
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    up_p = sorted(up & set(expr.columns))
    down_p = sorted(down & set(expr.columns))
    if len(up_p) < 3:
        raise ValueError(f"Too few up genes: {len(up_p)}")

    def norm_score(m, n):
        mn = (1 + n) / (2 * n_total)
        mx = 1 - mn
        return 2 * (m - mn) / (mx - mn) - 1

    up_s = norm_score(ranks[up_p].mean(axis=1), len(up_p))
    if down_p:
        down_s = norm_score((1 - ranks[down_p]).mean(axis=1), len(down_p))
        raw = up_s - down_s
    else:
        raw = up_s
    return raw - raw.median()


def compute_singscore_unidirectional(expr: pd.DataFrame, gene_set: set) -> pd.Series:
    """Up-only singscore (for Hallmark sets and single-gene MKI67)."""
    if expr.shape[1] < 5000:
        raise ValueError(f"Expression matrix has only {expr.shape[1]} genes — refusing to score.")
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    p = sorted(gene_set & set(expr.columns))
    if len(p) < 3:
        return pd.Series(index=expr.index, dtype=float)
    mn = (1 + len(p)) / (2 * n_total)
    mx = 1 - mn
    raw = 2 * (ranks[p].mean(axis=1) - mn) / (mx - mn) - 1
    return raw - raw.median()


def paired_delong(y_true: np.ndarray, scores_a: np.ndarray, scores_b: np.ndarray) -> tuple[float, float, float]:
    """Paired DeLong test for comparing two AUROCs on the SAME samples.
    Returns (auroc_a, auroc_b, p_value).
    Reference: DeLong et al. 1988; Sun & Xu 2014 (fast O(N log N) implementation).
    """
    y = np.asarray(y_true).astype(int)
    a = np.asarray(scores_a)
    b = np.asarray(scores_b)
    pos = y == 1
    neg = y == 0
    n_pos = pos.sum()
    n_neg = neg.sum()
    if n_pos < 2 or n_neg < 2:
        return float("nan"), float("nan"), float("nan")

    auc_a = roc_auc_score(y, a)
    auc_b = roc_auc_score(y, b)

    def midrank(x):
        order = np.argsort(x)
        ranks = np.empty(len(x))
        i = 0
        while i < len(x):
            j = i
            while j < len(x) - 1 and x[order[j + 1]] == x[order[i]]:
                j += 1
            avg_rank = 0.5 * (i + j) + 1
            for k in range(i, j + 1):
                ranks[order[k]] = avg_rank
            i = j + 1
        return ranks

    # V_a^{10}, V_a^{01}, V_b^{10}, V_b^{01}
    tx_a = midrank(a[pos]); ty_a = midrank(a[neg]); tz_a = midrank(a)
    tx_b = midrank(b[pos]); ty_b = midrank(b[neg]); tz_b = midrank(b)

    V_a_10 = (tz_a[pos] - tx_a) / n_neg
    V_a_01 = 1.0 - (tz_a[neg] - ty_a) / n_pos
    V_b_10 = (tz_b[pos] - tx_b) / n_neg
    V_b_01 = 1.0 - (tz_b[neg] - ty_b) / n_pos

    var_a = np.var(V_a_10, ddof=1) / n_pos + np.var(V_a_01, ddof=1) / n_neg
    var_b = np.var(V_b_10, ddof=1) / n_pos + np.var(V_b_01, ddof=1) / n_neg
    cov = np.cov(V_a_10, V_b_10, ddof=1)[0, 1] / n_pos + \
          np.cov(V_a_01, V_b_01, ddof=1)[0, 1] / n_neg
    var_diff = var_a + var_b - 2 * cov
    if var_diff <= 0:
        return auc_a, auc_b, float("nan")
    z = (auc_a - auc_b) / np.sqrt(var_diff)
    p = 2 * (1 - norm.cdf(abs(z)))
    return float(auc_a), float(auc_b), float(p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--expression", required=True, type=Path)
    ap.add_argument("--labels", required=True, type=Path)
    ap.add_argument("--hallmark", default=HALLMARK_PATH, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument("--n-random-baselines", type=int, default=10)
    args = ap.parse_args()

    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    labels_df = pd.read_csv(args.labels, sep="\t")
    common = sorted(set(expr.index.astype(str)) & set(labels_df["patient_id"].astype(str)))
    if len(common) < 10:
        raise SystemExit(f"Too few common patients: {len(common)}")
    expr = expr.loc[common]
    labels = labels_df.set_index("patient_id").loc[common]["pcr"].astype(int).to_numpy()

    # Hallmark gene sets
    hallmark = read_gmt(args.hallmark, pathways=PROLIFERATION_HALLMARKS)

    rows = []

    # Primary framework signature
    sens, res = load_signature_tsv(FROZEN_DIR / f"{PRIMARY_DRUG}.tsv")
    primary_scores = compute_singscore_bidirectional(expr, sens, res).loc[common].to_numpy()
    auc_primary = roc_auc_score(labels, primary_scores)
    rows.append({"label": f"PRIMARY_{PRIMARY_DRUG}", "type": "framework", "auroc": auc_primary,
                 "delta_vs_primary": 0.0, "paired_p": np.nan,
                 "signature_sha256": sha256_of(FROZEN_DIR / f"{PRIMARY_DRUG}.tsv")})

    # Trivial baselines: MKI67 single-gene + Hallmark E2F + Hallmark G2M + random
    if "MKI67" in expr.columns:
        mki67_scores = expr["MKI67"].rank(pct=True).to_numpy()
        auc, auc_b, p = paired_delong(labels, primary_scores, mki67_scores)
        rows.append({"label": "MKI67_single_gene", "type": "trivial_baseline",
                     "auroc": auc_b, "delta_vs_primary": auc - auc_b, "paired_p": p,
                     "signature_sha256": "na"})
    for hm_name, hm_set in hallmark.items():
        try:
            sc = compute_singscore_unidirectional(expr, hm_set).loc[common].to_numpy()
            auc, auc_b, p = paired_delong(labels, primary_scores, sc)
            rows.append({"label": hm_name, "type": "trivial_baseline",
                         "auroc": auc_b, "delta_vs_primary": auc - auc_b, "paired_p": p,
                         "signature_sha256": "na"})
        except Exception as e:
            rows.append({"label": hm_name, "type": "trivial_baseline",
                         "auroc": np.nan, "delta_vs_primary": np.nan, "paired_p": np.nan,
                         "signature_sha256": "na", "error": str(e)})

    # Random baselines (size-matched to primary signature)
    rng = np.random.default_rng(args.random_seed)
    universe = list(expr.columns)
    n_up, n_down = len(sens), len(res)
    random_aucs = []
    for i in range(args.n_random_baselines):
        sampled = rng.choice(universe, size=n_up + n_down, replace=False)
        up_r = set(sampled[:n_up]); down_r = set(sampled[n_up:])
        try:
            sc = compute_singscore_bidirectional(expr, up_r, down_r).loc[common].to_numpy()
            random_aucs.append(roc_auc_score(labels, sc))
        except Exception:
            continue
    if random_aucs:
        med_random_auc = float(np.median(random_aucs))
        rows.append({"label": f"random_baseline_median_of_{args.n_random_baselines}",
                     "type": "trivial_baseline", "auroc": med_random_auc,
                     "delta_vs_primary": auc_primary - med_random_auc, "paired_p": np.nan,
                     "signature_sha256": "na",
                     "details": f"random_aucs_range=[{min(random_aucs):.3f}, {max(random_aucs):.3f}]"})

    # Negative-control signatures (wrong-MoA TKIs on anthracycline cohort)
    for drug in NEGATIVE_CONTROL_DRUGS:
        sig_path = FROZEN_DIR / f"{drug}.tsv"
        if not sig_path.exists():
            rows.append({"label": f"NEG_CTRL_{drug}", "type": "negative_control",
                         "auroc": np.nan, "error": "signature_not_found",
                         "signature_sha256": "na"})
            continue
        sens_n, res_n = load_signature_tsv(sig_path)
        try:
            sc = compute_singscore_bidirectional(expr, sens_n, res_n).loc[common].to_numpy()
            auc_n = roc_auc_score(labels, sc)
            rows.append({"label": f"NEG_CTRL_{drug}", "type": "negative_control",
                         "auroc": auc_n, "delta_vs_primary": auc_primary - auc_n, "paired_p": np.nan,
                         "signature_sha256": sha256_of(sig_path)})
        except Exception as e:
            rows.append({"label": f"NEG_CTRL_{drug}", "type": "negative_control",
                         "auroc": np.nan, "error": str(e), "signature_sha256": "na"})

    df = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False, float_format="%.10g")

    # Pre-reg success-criterion evaluation
    baselines = df[df["type"] == "trivial_baseline"]
    if len(baselines) and "auroc" in baselines.columns:
        best_baseline_auc = baselines["auroc"].max()
        best_baseline_label = baselines.loc[baselines["auroc"].idxmax(), "label"]
        delta = auc_primary - best_baseline_auc
        best_baseline_row = baselines.loc[baselines["auroc"].idxmax()]
        paired_p = best_baseline_row.get("paired_p", np.nan)
        # Per pre-reg §8: paired DeLong p < 0.0083 AND Δ ≥ 0.03
        passed = (pd.notna(delta) and delta >= 0.03 and pd.notna(paired_p) and paired_p < 0.0083)
        summary = {
            "cohort": args.cohort,
            "primary_auroc": auc_primary,
            "best_baseline_label": best_baseline_label,
            "best_baseline_auroc": float(best_baseline_auc),
            "delta_vs_primary": float(delta) if pd.notna(delta) else None,
            "paired_p_vs_primary": float(paired_p) if pd.notna(paired_p) else None,
            "pre_reg_passed": bool(passed),
            "criterion": "delta >= 0.03 AND paired_DeLong p < 0.0083 (pre-reg §8)",
        }
    else:
        summary = {"cohort": args.cohort, "primary_auroc": auc_primary, "error": "no_baselines_run"}
    with open(args.output.with_suffix(".json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
