"""Phase 4 / Critical Gaps — Task B.

Combined-cohort meta-statistics for 2-gene oracle (MKI67 - ESR1 signed z) vs framework Dox signature.

Cohorts: GSE16446 (burned), GSE25066 (burned), GSE22226 (burned), GSE32646 (held-out — Phase 2)

Computes:
  1. Per-cohort paired DeLong test (proper, Sun & Xu 2014). Re-runs from raw expr/score
     to replace Phase 1 bootstrap permutation p-values.
  2. Holm-adjusted p across all 4 cohorts
  3. Stouffer's combined p (one-sided, weighting by sqrt(n))
  4. Fixed-effects meta AUROC-delta (inverse-variance weighted) + 95% CI
  5. Random-effects meta (DerSimonian-Laird) as sensitivity

Outputs:
  results/phase_4/critical_gaps/meta_stats.json
  results/phase_4/critical_gaps/meta_stats_summary.tsv
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

REPO = Path(".")
PHASE2 = REPO / "." / "results" / "phase_2"
OUT_DIR = REPO / "." / "results" / "phase_4" / "critical_gaps"

BURNED = ["GSE16446", "GSE25066", "GSE22226"]


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# -------- DeLong (Sun & Xu 2014 fast) --------
def _midrank(x: np.ndarray) -> np.ndarray:
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def delong_var_cov(scores_list, y):
    pos = (y == 1)
    m = int(pos.sum())
    n = int((~pos).sum())
    k = len(scores_list)
    V10 = np.zeros((k, m), dtype=float)
    V01 = np.zeros((k, n), dtype=float)
    aucs = np.zeros(k, dtype=float)
    for r, sc in enumerate(scores_list):
        X = sc[pos].astype(float)
        Y = sc[~pos].astype(float)
        TX = _midrank(X)
        TY = _midrank(Y)
        TZ = _midrank(np.concatenate([X, Y]))
        TZX = TZ[:m]
        TZY = TZ[m:]
        aucs[r] = (TZX.sum() / m - (m + 1) / 2.0) / n
        V10[r, :] = 1.0 - (TZX - TX) / n
        V01[r, :] = (TZY - TY) / m
    if k == 1:
        S10 = np.array([[np.var(V10[0], ddof=1)]])
        S01 = np.array([[np.var(V01[0], ddof=1)]])
    else:
        S10 = np.cov(V10, ddof=1)
        S01 = np.cov(V01, ddof=1)
    S = S10 / m + S01 / n
    return aucs, S


def delong_paired(y, s_a, s_b):
    aucs, S = delong_var_cov([s_a, s_b], y)
    delta = aucs[0] - aucs[1]
    var_delta = float(S[0, 0] + S[1, 1] - 2.0 * S[0, 1])
    if var_delta <= 0 or not np.isfinite(var_delta):
        return float(aucs[0]), float(aucs[1]), float(delta), float("nan"), float("nan"), float("nan")
    z = delta / np.sqrt(var_delta)
    p_two = 2.0 * (1.0 - norm.cdf(abs(z)))
    return float(aucs[0]), float(aucs[1]), float(delta), float(var_delta), float(z), float(p_two)


def load_burned_cohort(cohort: str):
    """Load expr, pcr, framework score for a Phase 1 burned cohort."""
    base = REPO / "external_data" / cohort
    expr = pd.read_csv(base / f"{cohort}_expression.tsv", sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    pcr = pd.read_csv(base / f"{cohort}_pcr.tsv", sep="\t")
    fw = pd.read_csv(base / "score_DOXORUBICIN.tsv", sep="\t")
    fw_col = "signature_score"
    return expr, pcr, fw, fw_col, {
        "expr_sha256": sha256_of(base / f"{cohort}_expression.tsv"),
        "pcr_sha256": sha256_of(base / f"{cohort}_pcr.tsv"),
        "framework_score_sha256": sha256_of(base / "score_DOXORUBICIN.tsv"),
    }


def load_gse32646():
    expr = pd.read_csv(PHASE2 / "GSE32646_expression.tsv", sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    pcr = pd.read_csv(PHASE2 / "GSE32646_pcr.tsv", sep="\t")
    fw = pd.read_csv(PHASE2 / "score_DOXORUBICIN_GSE32646.tsv", sep="\t")
    fw_col = "sig_dox"
    return expr, pcr, fw, fw_col, {
        "expr_sha256": sha256_of(PHASE2 / "GSE32646_expression.tsv"),
        "pcr_sha256": sha256_of(PHASE2 / "GSE32646_pcr.tsv"),
        "framework_score_sha256": sha256_of(PHASE2 / "score_DOXORUBICIN_GSE32646.tsv"),
    }


def compute_cohort(cohort: str, expr, pcr, fw, fw_col):
    common = (set(expr.index.astype(str)) &
              set(pcr["patient_id"].astype(str)) &
              set(fw["patient_id"].astype(str)))
    common = sorted(common)
    expr_c = expr.loc[common]
    pcr_c = pcr.set_index("patient_id").loc[common]["pcr"].astype(int).values
    fw_c = fw.set_index("patient_id").loc[common][fw_col].astype(float).values

    if "MKI67" not in expr_c.columns or "ESR1" not in expr_c.columns:
        raise RuntimeError(f"{cohort}: MKI67/ESR1 missing")

    mki = expr_c["MKI67"].astype(float)
    esr = expr_c["ESR1"].astype(float)
    zmki = (mki - mki.mean()) / (mki.std(ddof=0) + 1e-12)
    zesr = (esr - esr.mean()) / (esr.std(ddof=0) + 1e-12)
    oracle = (zmki - zesr).values  # signed sum: high MKI67 - high ESR1

    mask = np.isfinite(oracle) & np.isfinite(fw_c) & np.isfinite(pcr_c)
    y = pcr_c[mask]
    o = oracle[mask]
    f = fw_c[mask]

    auc_o, auc_f, delta, var_delta, z, p_two = delong_paired(y, o, f)
    n = int(mask.sum())
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())

    # one-sided p (2-gene > framework)
    if np.isfinite(z):
        p_one_sided = 1.0 - norm.cdf(z)
    else:
        p_one_sided = float("nan")

    return {
        "cohort": cohort,
        "n": n,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "framework_AUROC": auc_f,
        "two_gene_AUROC": auc_o,
        "delta": delta,
        "var_delta": var_delta,
        "se_delta": float(np.sqrt(var_delta)) if (np.isfinite(var_delta) and var_delta > 0) else float("nan"),
        "z": z,
        "DeLong_p_two_sided": p_two,
        "DeLong_p_one_sided_2g_gt_fw": p_one_sided,
    }


def holm(pvals):
    """Holm-Bonferroni adjusted p-values."""
    p = np.array(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n, dtype=float)
    running_max = 0.0
    for rank, idx in enumerate(order):
        v = (n - rank) * p[idx]
        running_max = max(running_max, v)
        adj[idx] = min(running_max, 1.0)
    return adj.tolist()


def stouffer(z_list, weights=None):
    z = np.array(z_list, dtype=float)
    if weights is None:
        w = np.ones_like(z)
    else:
        w = np.array(weights, dtype=float)
    z_c = float(np.sum(w * z) / np.sqrt(np.sum(w ** 2)))
    p_one = float(1.0 - norm.cdf(z_c))
    p_two = float(2.0 * (1.0 - norm.cdf(abs(z_c))))
    return z_c, p_one, p_two


def fixed_effects_meta(deltas, vars_):
    d = np.array(deltas, dtype=float)
    v = np.array(vars_, dtype=float)
    w = 1.0 / v
    mu = float(np.sum(w * d) / np.sum(w))
    se = float(np.sqrt(1.0 / np.sum(w)))
    lo, hi = mu - 1.96 * se, mu + 1.96 * se
    z = mu / se
    p_two = float(2.0 * (1.0 - norm.cdf(abs(z))))
    # Heterogeneity (Q, I^2)
    Q = float(np.sum(w * (d - mu) ** 2))
    df = len(d) - 1
    if df > 0:
        I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    else:
        I2 = float("nan")
    return {"delta": mu, "se": se, "ci_95": [lo, hi], "z": z, "p_two_sided": p_two,
            "Q": Q, "df": df, "I2": I2}


def random_effects_meta(deltas, vars_):
    """DerSimonian-Laird random-effects meta."""
    d = np.array(deltas, dtype=float)
    v = np.array(vars_, dtype=float)
    k = len(d)
    w_fe = 1.0 / v
    mu_fe = float(np.sum(w_fe * d) / np.sum(w_fe))
    Q = float(np.sum(w_fe * (d - mu_fe) ** 2))
    df = k - 1
    C = float(np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe))
    if df > 0 and C > 0:
        tau2 = max(0.0, (Q - df) / C)
    else:
        tau2 = 0.0
    w_re = 1.0 / (v + tau2)
    mu_re = float(np.sum(w_re * d) / np.sum(w_re))
    se_re = float(np.sqrt(1.0 / np.sum(w_re)))
    lo, hi = mu_re - 1.96 * se_re, mu_re + 1.96 * se_re
    z = mu_re / se_re
    p_two = float(2.0 * (1.0 - norm.cdf(abs(z))))
    return {"delta": mu_re, "se": se_re, "ci_95": [lo, hi], "z": z, "p_two_sided": p_two,
            "tau2": tau2, "Q": Q, "df": df}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {
        "script": str(me),
        "script_sha256": sha256_of(me),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "convention": {
            "oracle_formula": "z(MKI67) - z(ESR1)",
            "delong": "Sun & Xu 2014 paired (replaces Phase 1 bootstrap permutation)",
            "stouffer": "one-sided (2-gene > framework), weights = sqrt(n)",
            "holm": "Holm-Bonferroni across 4 cohorts",
            "fixed_effects": "inverse-variance weighted, var from DeLong covariance",
            "random_effects": "DerSimonian-Laird",
        },
        "cohort_sources": {},
    }

    per_cohort = []
    # --- 3 burned cohorts ---
    for c in BURNED:
        expr, pcr, fw, fw_col, sha = load_burned_cohort(c)
        log["cohort_sources"][c] = sha
        per_cohort.append(compute_cohort(c, expr, pcr, fw, fw_col))

    # --- GSE32646 ---
    expr, pcr, fw, fw_col, sha = load_gse32646()
    log["cohort_sources"]["GSE32646"] = sha
    per_cohort.append(compute_cohort("GSE32646", expr, pcr, fw, fw_col))

    # --- Holm across 4 ---
    pvals = [r["DeLong_p_two_sided"] for r in per_cohort]
    holm_adj = holm(pvals)
    for r, p_h in zip(per_cohort, holm_adj):
        r["Holm_p"] = p_h

    # --- Stouffer one-sided, weights sqrt(n) ---
    zs = [r["z"] for r in per_cohort if np.isfinite(r["z"])]
    ns = [r["n"] for r in per_cohort if np.isfinite(r["z"])]
    weights = np.sqrt(ns)
    z_stouffer, p_one_sided, p_two_sided = stouffer(zs, weights)

    # --- Fixed + random effects meta ---
    deltas = [r["delta"] for r in per_cohort if np.isfinite(r["var_delta"]) and r["var_delta"] > 0]
    vars_ = [r["var_delta"] for r in per_cohort if np.isfinite(r["var_delta"]) and r["var_delta"] > 0]
    fe = fixed_effects_meta(deltas, vars_)
    re_ = random_effects_meta(deltas, vars_)

    summary = {
        "per_cohort": per_cohort,
        "stouffer_combined": {
            "z": z_stouffer,
            "p_one_sided_2g_gt_fw": p_one_sided,
            "p_two_sided": p_two_sided,
            "weights": "sqrt(n)",
            "k_cohorts": int(len(zs)),
        },
        "fixed_effects_meta": fe,
        "random_effects_meta": re_,
    }

    log["summary"] = summary
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(OUT_DIR / "meta_stats.json", "w") as fh:
        json.dump(log, fh, indent=2)

    # --- Build TSV ---
    rows = []
    for r in per_cohort:
        rows.append({
            "row_type": "cohort",
            "cohort": r["cohort"],
            "n": r["n"],
            "framework_AUROC": r["framework_AUROC"],
            "two_gene_AUROC": r["two_gene_AUROC"],
            "delta": r["delta"],
            "DeLong_p": r["DeLong_p_two_sided"],
            "Holm_p": r["Holm_p"],
            "n_used": r["n"],
        })
    n_total = sum(r["n"] for r in per_cohort)
    rows.append({
        "row_type": "meta_fixed_effects",
        "cohort": "META_FE",
        "n": n_total,
        "framework_AUROC": np.nan,
        "two_gene_AUROC": np.nan,
        "delta": fe["delta"],
        "DeLong_p": fe["p_two_sided"],
        "Holm_p": np.nan,
        "n_used": n_total,
    })
    rows.append({
        "row_type": "meta_random_effects",
        "cohort": "META_RE",
        "n": n_total,
        "framework_AUROC": np.nan,
        "two_gene_AUROC": np.nan,
        "delta": re_["delta"],
        "DeLong_p": re_["p_two_sided"],
        "Holm_p": np.nan,
        "n_used": n_total,
    })

    tsv = pd.DataFrame(rows)
    tsv.to_csv(OUT_DIR / "meta_stats_summary.tsv", sep="\t", index=False)

    # Headline
    print("[TASK B] Per-cohort:")
    for r in per_cohort:
        print(f"  {r['cohort']:>10s} n={r['n']:4d}  fw={r['framework_AUROC']:.3f}  2g={r['two_gene_AUROC']:.3f}  "
              f"delta={r['delta']:+.3f}  DeLong_p={r['DeLong_p_two_sided']:.4g}  Holm={r['Holm_p']:.4g}")
    print(f"\n  Stouffer (one-sided, sqrt(n)): z={z_stouffer:.3f}  p={p_one_sided:.4g}")
    print(f"  Fixed-effects meta: delta={fe['delta']:+.4f} [{fe['ci_95'][0]:+.4f}, {fe['ci_95'][1]:+.4f}]  "
          f"p={fe['p_two_sided']:.4g}  I2={fe['I2']:.2f}")
    print(f"  Random-effects meta: delta={re_['delta']:+.4f} [{re_['ci_95'][0]:+.4f}, {re_['ci_95'][1]:+.4f}]  "
          f"tau2={re_['tau2']:.4f}")


if __name__ == "__main__":
    main()
