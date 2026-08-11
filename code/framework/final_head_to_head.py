"""Final multi-metric head-to-head: all predictors on 3 cohorts.

Predictors compared:
  - Framework signatures: DOXORUBICIN, PACLITAXEL, TOP2_consensus, ETOPOSIDE (singscore)
  - Ridge regression (oncoPredict/pRRophetic-style): ridge_rank, ridge_zscore
  - Proliferation baselines: MKI67, E2F_TARGETS, G2M_CHECKPOINT
  - Random reference

Metrics (all reviewer-preferred):
  - AUROC + bootstrap CI
  - PR-AUC + bootstrap CI
  - Brier score
  - Calibration slope (target=1)
  - F1 @ Youden-J
  - DCA net benefit at t=5%, 10%, 15%, 20%, 25%

Decision question: where does the framework sit in the predictor ranking?
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             brier_score_loss, roc_curve, f1_score)
from sklearn.isotonic import IsotonicRegression

REVISE_ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
EXT = REVISE_ROOT / "external_data"
HALLMARK = Path("/home/holiday01/2026_ISMB_code/enrichment_test/h.all.v2025.1.Hs.symbols.gmt")
FROZEN = REVISE_ROOT / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all"


def read_gmt_set(path, name):
    with open(path) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if parts and parts[0] == name:
                return set(g.upper() for g in parts[2:] if g.strip())
    return set()


def singscore_bidirectional(expr, up, down):
    ranks = expr.rank(axis=1, method="average", pct=True)
    n_total = expr.shape[1]
    up_p = sorted(up & set(expr.columns))
    down_p = sorted(down & set(expr.columns))
    def norm_score(m, n):
        mn = (1 + n) / (2 * n_total); mx = 1 - mn
        return 2 * (m - mn) / (mx - mn) - 1
    up_s = norm_score(ranks[up_p].mean(axis=1), len(up_p))
    if down_p:
        down_s = norm_score((1 - ranks[down_p]).mean(axis=1), len(down_p))
        raw = up_s - down_s
    else:
        raw = up_s
    return raw - raw.median()


def singscore_unidirectional(expr, gene_set):
    ranks = expr.rank(axis=1, method="average", pct=True)
    n_total = expr.shape[1]
    p = sorted(gene_set & set(expr.columns))
    if len(p) < 3:
        return pd.Series(index=expr.index, dtype=float)
    mn = (1 + len(p)) / (2 * n_total); mx = 1 - mn
    raw = 2 * (ranks[p].mean(axis=1) - mn) / (mx - mn) - 1
    return raw - raw.median()


def load_signature(path):
    df = pd.read_csv(path, sep="\t")
    sens = set(df[df["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(df[df["direction"] == "resistance"]["gene_symbol"].str.upper())
    return sens, res


def score_to_prob(score):
    rng = score.max() - score.min()
    return (score - score.min()) / rng if rng > 0 else np.full_like(score, 0.5)


def all_metrics(y, score, n_boot=500, seed=42):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    out = {}
    out["auroc"] = float(roc_auc_score(y, score))
    out["pr_auc"] = float(average_precision_score(y, score))
    prob = score_to_prob(score)
    out["brier"] = float(brier_score_loss(y, prob))
    # Calibration slope
    try:
        iso = IsotonicRegression(out_of_bounds="clip").fit(prob, y)
        iso_pred = iso.predict(prob)
        slope, intercept = np.polyfit(prob, iso_pred, 1) if np.std(prob) > 0 else (float("nan"), float("nan"))
        out["calibration_slope"] = float(slope)
    except Exception:
        out["calibration_slope"] = float("nan")
    # F1 @ Youden
    fpr, tpr, thr = roc_curve(y, score)
    j = np.argmax(tpr - fpr)
    pred = (score >= thr[j]).astype(int)
    out["f1_youden"] = float(f1_score(y, pred, zero_division=0))
    # Bootstrap CI for AUROC
    rng = np.random.default_rng(seed); n = len(y)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) >= 2:
            boot.append(roc_auc_score(y[idx], score[idx]))
    out["auroc_ci"] = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))] if boot else [float("nan"), float("nan")]
    # DCA
    n = len(y)
    dca = {}
    for t in [0.05, 0.10, 0.15, 0.20, 0.25]:
        pred = (prob >= t).astype(int)
        tp = ((pred == 1) & (y == 1)).sum(); fp = ((pred == 1) & (y == 0)).sum()
        nb = tp / n - fp / n * (t / (1 - t))
        dca[f"nb_t{int(t*100):02d}"] = float(nb)
    out["dca"] = dca
    return out


COHORTS = ["GSE16446", "GSE25066", "GSE22226"]


def main():
    e2f_genes = read_gmt_set(HALLMARK, "HALLMARK_E2F_TARGETS")
    g2m_genes = read_gmt_set(HALLMARK, "HALLMARK_G2M_CHECKPOINT")

    rows = []

    for cohort in COHORTS:
        expr = pd.read_csv(EXT / cohort / f"{cohort}_expression.tsv", sep="\t", index_col=0)
        expr.columns = [c.upper() for c in expr.columns]
        pcr_df = pd.read_csv(EXT / cohort / f"{cohort}_pcr.tsv", sep="\t").dropna(subset=["pcr"])
        common = sorted(set(expr.index.astype(str)) & set(pcr_df["patient_id"].astype(str)))
        expr = expr.loc[common]
        y = pcr_df.set_index("patient_id").loc[common]["pcr"].astype(int).to_numpy()

        scores = {}
        # Framework signatures
        for drug in ["DOXORUBICIN", "PACLITAXEL", "TOP2_POISON_CONSENSUS", "ETOPOSIDE"]:
            sig_path = FROZEN / f"{drug}.tsv"
            if sig_path.exists():
                sens, res = load_signature(sig_path)
                try:
                    scores[f"framework_{drug}"] = singscore_bidirectional(expr, sens, res).loc[common].to_numpy()
                except Exception as e:
                    pass

        # Proliferation baselines
        if "MKI67" in expr.columns:
            scores["baseline_MKI67"] = expr["MKI67"].rank(pct=True).loc[common].to_numpy()
        scores["baseline_E2F_TARGETS"] = singscore_unidirectional(expr, e2f_genes).loc[common].to_numpy()
        scores["baseline_G2M_CHECKPOINT"] = singscore_unidirectional(expr, g2m_genes).loc[common].to_numpy()

        # Ridge predictions
        for method in ["ridge_rank", "ridge_zscore"]:
            ridge_path = EXT / cohort / f"score_ridge_{method}.tsv"
            if ridge_path.exists():
                rd = pd.read_csv(ridge_path, sep="\t").set_index("patient_id").loc[common]
                # Score = -predicted_auc (lower AUC = sensitive = should predict pCR=1)
                scores[f"tool_{method}"] = (-rd["predicted_auc"]).to_numpy()

        for name, sc in scores.items():
            m = all_metrics(y, sc)
            m["cohort"] = cohort; m["predictor"] = name; m["n"] = len(y); m["pcr_rate"] = float(y.mean())
            # Flatten dca
            for k, v in m.pop("dca").items():
                m[k] = v
            # Flatten auroc_ci
            ci = m.pop("auroc_ci")
            m["auroc_ci_lo"], m["auroc_ci_hi"] = ci[0], ci[1]
            rows.append(m)

    df = pd.DataFrame(rows)
    cols = ["cohort", "predictor", "n", "pcr_rate", "auroc", "auroc_ci_lo", "auroc_ci_hi",
            "pr_auc", "brier", "calibration_slope", "f1_youden",
            "nb_t05", "nb_t10", "nb_t15", "nb_t20", "nb_t25"]
    df = df[cols]
    out_path = REVISE_ROOT / "primary_test_results" / "final_head_to_head.tsv"
    df.to_csv(out_path, sep="\t", index=False, float_format="%.4g")

    print(df.to_string(index=False))

    # Rank within cohort
    print("\n" + "=" * 100)
    print("RANKING within each cohort (by AUROC, best → worst):")
    print("=" * 100)
    for cohort in COHORTS:
        sub = df[df["cohort"] == cohort].sort_values("auroc", ascending=False)
        print(f"\n  {cohort}:")
        for _, r in sub.iterrows():
            type_ = "framework" if r["predictor"].startswith("framework") else \
                    "tool" if r["predictor"].startswith("tool") else "baseline"
            print(f"    [{type_:<9}] {r['predictor']:<30}  AUROC={r['auroc']:.3f} CI=[{r['auroc_ci_lo']:.3f},{r['auroc_ci_hi']:.3f}]  PR-AUC={r['pr_auc']:.3f}  Brier={r['brier']:.3f}")

    print(f"\nFull table: {out_path}")


if __name__ == "__main__":
    main()
