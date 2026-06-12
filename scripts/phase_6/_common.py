"""Phase 6 shared helpers: singscore (full-transcriptome rank), alias mapping,
axis decomposition residual, bootstrap AUROC, 3x3 verdict, sha256.

Identical mathematics to Phase 0 / Phase 2 / Phase 4 to keep the
cross-signature test apples-to-apples.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_alias_map(path: Path) -> dict[str, str]:
    """Return alias_symbol(upper) -> canonical_symbol(upper)."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols = {c.lower(): c for c in df.columns}
    m: dict[str, str] = {}
    if "alias_symbol" in cols and "canonical_symbol" in cols:
        for _, row in df.iterrows():
            alias = str(row[cols["alias_symbol"]]).upper().strip()
            canon = str(row[cols["canonical_symbol"]]).upper().strip()
            if alias and canon and alias != canon:
                m.setdefault(alias, canon)
    return m


def resolve_genes(wanted: Iterable[str], cohort_genes: set[str], alias_map: dict[str, str]):
    cohort_upper = {g.upper() for g in cohort_genes}
    present, mapping, missing = [], {}, []
    for g in wanted:
        gU = g.upper()
        if gU in cohort_upper:
            present.append(gU)
            mapping[gU] = gU
            continue
        if gU in alias_map and alias_map[gU] in cohort_upper:
            present.append(alias_map[gU])
            mapping[gU] = alias_map[gU]
            continue
        hit = None
        for alias, approved in alias_map.items():
            if approved == gU and alias in cohort_upper:
                hit = alias
                break
        if hit:
            present.append(hit)
            mapping[gU] = hit
            continue
        missing.append(gU)
    # de-duplicate, preserving order
    seen = set()
    present = [g for g in present if not (g in seen or seen.add(g))]
    return present, mapping, missing


def compute_singscore(expr: pd.DataFrame, up_genes: list[str], down_genes: list[str]) -> pd.Series:
    """Full-transcriptome-rank singscore (same as Phase 0 / Phase 3)."""
    n_genes_total = expr.shape[1]
    if n_genes_total < 5000:
        raise ValueError(f"Expression matrix has only {n_genes_total} genes (<5000).")
    rank_matrix = expr.rank(axis=1, method="average", pct=True)

    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    # Match Phase 3 mammaprint_score.py exactly: raw = up_score - down_score
    # where down_score is normalized over (1 - rank). This is the same formula
    # used in Phase 0 score_external_cohort and Phase 3, preserved here so
    # cross-signature comparison is apples-to-apples within this analysis.
    if up_genes:
        up_score = normalize(rank_matrix[up_genes].mean(axis=1), len(up_genes))
    else:
        up_score = pd.Series(0.0, index=expr.index)
    if down_genes:
        down_score = normalize((1 - rank_matrix[down_genes]).mean(axis=1), len(down_genes))
    else:
        down_score = pd.Series(0.0, index=expr.index)
    raw = up_score - down_score
    return raw - raw.median()


def auroc_with_boot_ci(scores, labels, n_boot=1000, seed=42):
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels)
    if len(np.unique(labels)) < 2:
        return float("nan"), float("nan"), float("nan"), []
    obs = roc_auc_score(labels, scores)
    rng = np.random.default_rng(seed)
    n = len(labels)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(labels[idx])) < 2:
            continue
        boot.append(roc_auc_score(labels[idx], scores[idx]))
    if len(boot) < 10:
        return float(obs), float("nan"), float("nan"), boot
    return float(obs), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)), boot


def fit_residual(df: pd.DataFrame, sig_col: str, axis_cols: list[str], label_col: str):
    needed = [sig_col, label_col] + axis_cols
    sub = df.dropna(subset=needed).copy()
    if len(sub) < 20:
        return None, None
    X = sub[axis_cols].values
    y = sub[sig_col].values
    lr = LinearRegression()
    lr.fit(X, y)
    sub["residual"] = y - lr.predict(X)
    coef = {c: float(v) for c, v in zip(axis_cols, lr.coef_)}
    coef["intercept"] = float(lr.intercept_)
    return sub, coef


def verdict_3x3(point: float, ci_lo: float):
    if not (math.isfinite(point) and math.isfinite(ci_lo)):
        return "INDETERMINATE_NAN"
    # consider absolute-distance-from-0.5 so that an inverse-signal residual
    # still counts as "signal survives" (FALSIFIED).
    info_point = abs(point - 0.5)
    info_ci_lo = abs(ci_lo - 0.5)  # not used for the standard rule;
    # Apply the pre-specified rule on the raw (signed) AUROC:
    if point <= 0.55 and ci_lo <= 0.55:
        return "H-b CONFIRMED (no residual signal after deconfounding)"
    if point > 0.60 and ci_lo > 0.55:
        return "H-b FALSIFIED (residual signal survives)"
    return "GRAY ZONE"
