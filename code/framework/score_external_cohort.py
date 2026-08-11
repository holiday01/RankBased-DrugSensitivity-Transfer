"""Stage 3: Score an external cohort using a frozen CCLE-derived signature.

Applies a frozen signature TSV (from Stage 0) to an external cohort's expression
matrix via full-transcriptome rank-based singscore. Outputs per-patient scores
+ AUROC vs binary endpoint + bootstrap CI + permutation null p-values.

CRITICAL implementation invariant:
  Rank is computed over the FULL gene set of the external cohort, NOT
  subset to the signature genes. This honors the singscore foundational
  assumption (see manuscript Methods Eq. 1) and the audit finding that
  `task_common.score_signature_samples` has a subset-rank bug.

Usage:
  python score_external_cohort.py \
      --signature frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv \
      --expression <cohort_expression.tsv>   \  # rows = patient_id, cols = HGNC_symbol
      --labels <cohort_pcr_labels.tsv>       \  # cols = patient_id, pcr (0/1)
      --output <result.tsv>                  \
      --stratum-col <subtype>                \  # optional, for stratified AUROC

Output columns:
  patient_id, signature_score, ...

Plus a JSON sidecar with: signature_sha256, signature_genes_up_found_n,
  signature_genes_down_found_n, total_genes_ranked_n, n_patients_scored,
  auroc, auroc_ci_low, auroc_ci_high, perm_p_label_shuffle,
  perm_p_gene_set_shuffle.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_signature(path: Path) -> tuple[set, set, str, dict]:
    """Returns (sensitivity_genes, resistance_genes, score_col, info)."""
    df = pd.read_csv(path, sep="\t")
    required_cols = {"gene_symbol", "direction"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Signature {path} missing required columns")
    sens = set(df[df["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(df[df["direction"] == "resistance"]["gene_symbol"].str.upper())
    score_col = next((c for c in df.columns if c not in {"gene_symbol", "direction", "rank_in_direction"}), None)
    info = {
        "signature_path": str(path),
        "signature_sha256": sha256_of(path),
        "n_sensitivity_genes": len(sens),
        "n_resistance_genes": len(res),
        "score_col_in_signature": score_col,
    }
    return sens, res, score_col, info


def compute_singscore(
    expr: pd.DataFrame,
    up_genes: set,
    down_genes: set,
) -> tuple[pd.Series, dict]:
    """Full-transcriptome rank-based singscore (matches manuscript Methods Eq. 1).

    expr: DataFrame, index = patient_id, columns = HGNC_symbol (UPPERCASE).
          MUST be the full cohort transcriptome (not subsetted to signature).
    """
    # Normalize gene-symbol case
    expr_cols_upper = {c: c.upper() for c in expr.columns}
    expr = expr.rename(columns=expr_cols_upper)

    n_genes_total = expr.shape[1]
    if n_genes_total < 5000:
        # Defensive: refuse to score on subsetted matrix (Stage 0 audit issue)
        raise ValueError(
            f"Expression matrix has only {n_genes_total} genes — looks like a subset. "
            "Stage 3 requires full-transcriptome (≥5000 genes) for valid rank."
        )

    # Per-sample percentile ranks (axis=1 = within each patient)
    rank_matrix = expr.rank(axis=1, method="average", pct=True)

    up_present = sorted(up_genes & set(expr.columns))
    down_present = sorted(down_genes & set(expr.columns))

    if len(up_present) < 3 or len(down_present) < 3:
        raise ValueError(
            f"Too few signature genes present: up={len(up_present)}, down={len(down_present)}. "
            "Verify probe-to-symbol mapping + HGNC alias resolution."
        )

    # Theoretical normalization range (matches manuscript Eq. 1)
    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    up_score = normalize(rank_matrix[up_present].mean(axis=1), len(up_present))
    down_score = normalize((1 - rank_matrix[down_present]).mean(axis=1), len(down_present))
    raw = up_score - down_score
    centered = raw - raw.median()

    diag = {
        "n_genes_in_expression": n_genes_total,
        "n_up_genes_in_signature": len(up_genes),
        "n_up_genes_present_in_cohort": len(up_present),
        "n_down_genes_in_signature": len(down_genes),
        "n_down_genes_present_in_cohort": len(down_present),
        "up_genes_missing": sorted(up_genes - set(expr.columns)),
        "down_genes_missing": sorted(down_genes - set(expr.columns)),
    }
    return centered, diag


def auroc_with_bootstrap_ci(scores: np.ndarray, labels: np.ndarray, n_boot: int = 1000, seed: int = 42):
    if len(np.unique(labels)) < 2:
        return float("nan"), float("nan"), float("nan")
    obs = roc_auc_score(labels, scores)
    rng = np.random.default_rng(seed)
    boot = []
    n = len(labels)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(labels[idx])) < 2:
            continue
        boot.append(roc_auc_score(labels[idx], scores[idx]))
    if len(boot) < 10:
        return float(obs), float("nan"), float("nan")
    return float(obs), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def permutation_p_within_stratum(scores: np.ndarray, labels: np.ndarray, strata: np.ndarray | None,
                                  n_perm: int = 1000, seed: int = 1337):
    """Pre-reg null #1: shuffle labels WITHIN strata only.

    REFUSES to fall back to global shuffle if strata is None. This is the audit-
    flagged failure mode that produced the original 0.752 base-rate artifact in
    BIOADV-2026-206.
    """
    if strata is None or len(np.unique(strata)) < 1:
        raise ValueError(
            "permutation_p_within_stratum REQUIRES a stratum variable; refusing to "
            "fall back to global shuffle (would replicate the base-rate artifact null). "
            "For single-stratum cohorts, pass a constant stratum array to make the "
            "degeneration explicit, OR use a within-cohort biological covariate "
            "(PAM50 basal class for GSE16446)."
        )
    obs = roc_auc_score(labels, scores) if len(np.unique(labels)) >= 2 else float("nan")
    if not np.isfinite(obs):
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    perm_aucs = []
    for _ in range(n_perm):
        perm_labels = labels.copy()
        for s in np.unique(strata):
            idx = np.where(strata == s)[0]
            perm_labels[idx] = rng.permutation(labels[idx])
        if len(np.unique(perm_labels)) >= 2:
            perm_aucs.append(roc_auc_score(perm_labels, scores))
    if not perm_aucs:
        return obs, float("nan")
    p = (np.sum(np.array(perm_aucs) >= obs) + 1) / (len(perm_aucs) + 1)
    return obs, float(p)


def permutation_p_gene_set(expr: pd.DataFrame, labels: np.ndarray,
                            up_n: int, down_n: int, gene_universe: list,
                            obs_auroc: float, n_perm: int = 1000, seed: int = 1338):
    """Pre-reg null #2: replace signature with random same-size gene set."""
    rng = np.random.default_rng(seed)
    perm_aucs = []
    n_genes_total = expr.shape[1]
    rank_matrix = expr.rank(axis=1, method="average", pct=True)
    expr_genes = set(expr.columns)
    universe_filtered = [g for g in gene_universe if g in expr_genes]
    if len(universe_filtered) < up_n + down_n + 10:
        return float("nan")
    for _ in range(n_perm):
        sampled = rng.choice(universe_filtered, size=up_n + down_n, replace=False)
        up_rand = set(sampled[:up_n])
        down_rand = set(sampled[up_n:])
        try:
            up_score = rank_matrix[list(up_rand)].mean(axis=1)
            down_score = (1 - rank_matrix[list(down_rand)]).mean(axis=1)
            score = (up_score - down_score)
            score = score - score.median()
            if len(np.unique(labels)) >= 2:
                perm_aucs.append(roc_auc_score(labels, score.to_numpy()))
        except Exception:
            continue
    if not perm_aucs:
        return float("nan")
    return (np.sum(np.array(perm_aucs) >= obs_auroc) + 1) / (len(perm_aucs) + 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signature", required=True, type=Path)
    ap.add_argument("--expression", required=True, type=Path,
                    help="TSV: rows=patient_id, cols=HGNC_symbol, values=expression. Full transcriptome required.")
    ap.add_argument("--labels", required=True, type=Path,
                    help="TSV with cols: patient_id, pcr (0/1). Optional stratum col.")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--stratum-col", default=None)
    ap.add_argument("--gene-universe", default=None, type=Path,
                    help="Optional: TSV with one column of gene_symbols for permutation null #2.")
    ap.add_argument("--n-bootstrap", type=int, default=1000)
    ap.add_argument("--n-permutation", type=int, default=1000)
    ap.add_argument("--skip-perm1", action="store_true",
                    help="Skip within-stratum label permutation (Null #1). Use when cohort has no valid stratification variable (e.g., GSE16446 PAM50 fallback per pre-reg §6.2).")
    ap.add_argument("--seed", type=int, default=42,
                    help="Bootstrap seed (pre-reg: 42). Permutation null #1 seed=1337, #2 seed=1338.")
    args = ap.parse_args()

    sens, res, score_col, info = load_signature(args.signature)
    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    labels_df = pd.read_csv(args.labels, sep="\t")
    if "patient_id" not in labels_df.columns or "pcr" not in labels_df.columns:
        raise ValueError("labels file must have columns: patient_id, pcr")

    # Align
    common = sorted(set(expr.index.astype(str)) & set(labels_df["patient_id"].astype(str)))
    if len(common) < 10:
        raise ValueError(f"Too few common patients: {len(common)}")
    expr_aligned = expr.loc[common]
    labels_aligned = labels_df.set_index("patient_id").loc[common]

    scores, diag = compute_singscore(expr_aligned, sens, res)
    info.update(diag)
    info["n_patients_scored"] = len(common)

    y = labels_aligned["pcr"].astype(int).to_numpy()
    s = scores.loc[common].to_numpy()

    strata = labels_aligned[args.stratum_col].to_numpy() if args.stratum_col else None

    auroc, ci_lo, ci_hi = auroc_with_bootstrap_ci(s, y, args.n_bootstrap, args.seed)
    if args.skip_perm1:
        perm_p_label = None
    else:
        _, perm_p_label = permutation_p_within_stratum(s, y, strata, args.n_permutation, seed=1337)

    if args.gene_universe and args.gene_universe.exists():
        gu = pd.read_csv(args.gene_universe, sep="\t")
        gene_universe_list = gu.iloc[:, 0].astype(str).str.upper().tolist()
        perm_p_gs = permutation_p_gene_set(
            expr_aligned, y, len(sens), len(res), gene_universe_list, auroc,
            args.n_permutation, args.seed,
        )
    else:
        perm_p_gs = None

    info["auroc"] = auroc
    info["auroc_ci_low"] = ci_lo
    info["auroc_ci_high"] = ci_hi
    info["perm_p_label_within_stratum"] = perm_p_label
    info["perm_p_gene_set_shuffle"] = perm_p_gs
    info["seed"] = args.seed

    out_scores = pd.DataFrame({"patient_id": common, "signature_score": s, "pcr": y})
    if args.stratum_col:
        out_scores[args.stratum_col] = strata
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_scores.to_csv(args.output, sep="\t", index=False, float_format="%.10g")
    with open(args.output.with_suffix(".json"), "w") as f:
        json.dump(info, f, indent=2)
    pp_label_str = f"{perm_p_label:.4f}" if perm_p_label is not None else "skipped"
    pp_gs_str = f"{perm_p_gs:.4f}" if perm_p_gs is not None and not pd.isna(perm_p_gs) else "na"
    print(f"AUROC: {auroc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]; "
          f"perm_p_label={pp_label_str}; perm_p_gene_set={pp_gs_str}; n={len(common)}")


if __name__ == "__main__":
    main()
