"""Multi-metric fair comparison: framework signatures vs proliferation baselines.

Addresses R1-M5 + R2-2: AUROC alone is not the right metric for clinical
biomarker evaluation. Pre-reg §6.1 + §6.4 require PR-AUC, Brier, calibration,
and DCA. The §8 baseline comparison only ran AUROC, so this script fills the
gap.

Metrics computed (per cohort × per predictor):
  - AUROC + 95% bootstrap CI (already done; re-computed here for completeness)
  - PR-AUC (average precision) + 95% bootstrap CI
  - Brier score
  - Calibration slope (linear fit to isotonic predictions)
  - F1 + precision + recall at Youden-J threshold
  - DCA: net benefit at threshold range 5%-40% (Vickers)

Paired comparisons (framework vs each baseline, same patients):
  - Paired bootstrap Δ-PR-AUC + 95% CI
  - Paired bootstrap Δ-Brier + 95% CI
  - DCA net-benefit difference at clinically relevant thresholds (5%, 10%, 15%, 20%)

Decision rule: framework "wins" on a metric in a cohort if
  - Δ_metric > 0 AND 95% CI excludes 0
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             brier_score_loss, roc_curve, f1_score,
                             precision_score, recall_score)
from sklearn.isotonic import IsotonicRegression

REVISE_ROOT = Path(".")
EXT = REVISE_ROOT / "external_data"
HALLMARK = Path("./enrichment_test/h.all.v2025.1.Hs.symbols.gmt")
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


def score_to_prob(score):
    """Min-max scaling to [0,1] for probabilistic metrics (Brier, DCA)."""
    rng = score.max() - score.min()
    if rng > 0:
        return (score - score.min()) / rng
    return np.full_like(score, 0.5)


def all_metrics(y, score, n_boot=1000, seed=42):
    """Returns dict of metric → value (+ CI for select metrics)."""
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    out = {}
    out["auroc"] = float(roc_auc_score(y, score))
    out["pr_auc"] = float(average_precision_score(y, score))
    prob = score_to_prob(score)
    out["brier"] = float(brier_score_loss(y, prob))

    # Calibration slope (linear fit to isotonic predictions)
    try:
        iso = IsotonicRegression(out_of_bounds="clip").fit(prob, y)
        iso_pred = iso.predict(prob)
        if np.std(prob) > 0:
            slope, intercept = np.polyfit(prob, iso_pred, 1)
            out["calibration_slope"] = float(slope)
            out["calibration_intercept"] = float(intercept)
        else:
            out["calibration_slope"] = float("nan")
            out["calibration_intercept"] = float("nan")
    except Exception:
        out["calibration_slope"] = float("nan")
        out["calibration_intercept"] = float("nan")

    # F1 @ Youden-J
    try:
        fpr, tpr, thr = roc_curve(y, score)
        j_idx = np.argmax(tpr - fpr)
        t = float(thr[j_idx])
        pred = (score >= t).astype(int)
        out["f1_youden"] = float(f1_score(y, pred, zero_division=0))
        out["precision_youden"] = float(precision_score(y, pred, zero_division=0))
        out["recall_youden"] = float(recall_score(y, pred, zero_division=0))
    except Exception:
        out["f1_youden"] = out["precision_youden"] = out["recall_youden"] = float("nan")

    # Bootstrap CI for AUROC + PR-AUC
    rng = np.random.default_rng(seed)
    n = len(y)
    auroc_boot, prauc_boot = [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        auroc_boot.append(roc_auc_score(y[idx], score[idx]))
        prauc_boot.append(average_precision_score(y[idx], score[idx]))
    out["auroc_ci"] = [float(np.percentile(auroc_boot, 2.5)), float(np.percentile(auroc_boot, 97.5))] if auroc_boot else [float("nan"), float("nan")]
    out["pr_auc_ci"] = [float(np.percentile(prauc_boot, 2.5)), float(np.percentile(prauc_boot, 97.5))] if prauc_boot else [float("nan"), float("nan")]

    # DCA at key thresholds
    thresholds = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]
    dca = {}
    for t in thresholds:
        pred = (prob >= t).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        nb_model = tp / n - fp / n * (t / (1 - t))
        tp_all = int((y == 1).sum())
        fp_all = int((y == 0).sum())
        nb_all = tp_all / n - fp_all / n * (t / (1 - t))
        dca[f"nb_model_t{int(t*100):02d}"] = float(nb_model)
        dca[f"nb_treat_all_t{int(t*100):02d}"] = float(nb_all)
        dca[f"nb_diff_vs_treat_all_t{int(t*100):02d}"] = float(nb_model - nb_all)
    out["dca"] = dca

    return out


def paired_bootstrap_delta(y, score_a, score_b, metric_fn, n_boot=1000, seed=42):
    """Paired bootstrap Δ(metric_a - metric_b). Returns observed, ci_low, ci_high."""
    y = np.asarray(y).astype(int)
    a = np.asarray(score_a, dtype=float)
    b = np.asarray(score_b, dtype=float)
    obs = metric_fn(y, a) - metric_fn(y, b)
    rng = np.random.default_rng(seed)
    n = len(y)
    deltas = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            d = metric_fn(y[idx], a[idx]) - metric_fn(y[idx], b[idx])
            if np.isfinite(d):
                deltas.append(d)
        except Exception:
            continue
    if len(deltas) < 10:
        return float(obs), float("nan"), float("nan")
    return float(obs), float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))


COHORTS = {
    "GSE16446": dict(expr=EXT/"GSE16446/GSE16446_expression.tsv", pcr=EXT/"GSE16446/GSE16446_pcr.tsv"),
    "GSE25066": dict(expr=EXT/"GSE25066/GSE25066_expression.tsv", pcr=EXT/"GSE25066/GSE25066_pcr.tsv"),
    "GSE22226": dict(expr=EXT/"GSE22226/GSE22226_expression.tsv", pcr=EXT/"GSE22226/GSE22226_pcr.tsv"),
}

DRUG_SIGS = {
    "DOXORUBICIN":    FROZEN / "DOXORUBICIN.tsv",
    "PACLITAXEL":     FROZEN / "PACLITAXEL.tsv",
    "TOP2_consensus": FROZEN / "TOP2_POISON_CONSENSUS.tsv",
    "ETOPOSIDE":      FROZEN / "ETOPOSIDE.tsv",
}


def load_signature(path):
    df = pd.read_csv(path, sep="\t")
    sens = set(df[df["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(df[df["direction"] == "resistance"]["gene_symbol"].str.upper())
    return sens, res


def main():
    e2f_genes = read_gmt_set(HALLMARK, "HALLMARK_E2F_TARGETS")
    g2m_genes = read_gmt_set(HALLMARK, "HALLMARK_G2M_CHECKPOINT")

    all_results = {}
    paired_results = []

    for cohort_name, paths in COHORTS.items():
        print(f"\n=== {cohort_name} ===")
        expr = pd.read_csv(paths["expr"], sep="\t", index_col=0)
        expr.columns = [c.upper() for c in expr.columns]
        pcr_df = pd.read_csv(paths["pcr"], sep="\t").dropna(subset=["pcr"])
        common = sorted(set(expr.index.astype(str)) & set(pcr_df["patient_id"].astype(str)))
        expr = expr.loc[common]
        labels_df = pcr_df.set_index("patient_id").loc[common]
        y = labels_df["pcr"].astype(int).to_numpy()
        print(f"  n={len(y)}, pcr_rate={y.mean():.3f}")

        # Scores: framework signatures + baselines
        scores = {}
        for drug_name, sig_path in DRUG_SIGS.items():
            sens, res = load_signature(sig_path)
            try:
                scores[drug_name] = singscore_bidirectional(expr, sens, res).loc[common].to_numpy()
            except Exception as e:
                print(f"  {drug_name}: ERROR {e}")
        scores["MKI67"] = expr["MKI67"].rank(pct=True).loc[common].to_numpy() if "MKI67" in expr.columns else None
        scores["E2F_TARGETS"] = singscore_unidirectional(expr, e2f_genes).loc[common].to_numpy() if e2f_genes else None
        scores["G2M_CHECKPOINT"] = singscore_unidirectional(expr, g2m_genes).loc[common].to_numpy() if g2m_genes else None
        scores = {k: v for k, v in scores.items() if v is not None}

        # Compute all metrics
        cohort_results = {}
        for name, score in scores.items():
            try:
                m = all_metrics(y, score)
                cohort_results[name] = m
                print(f"  {name:<18}: AUROC={m['auroc']:.3f} PR-AUC={m['pr_auc']:.3f} Brier={m['brier']:.3f} F1={m['f1_youden']:.3f}")
            except Exception as e:
                print(f"  {name}: ERROR {e}")
        all_results[cohort_name] = dict(n=int(len(y)), pcr_rate=float(y.mean()), per_predictor=cohort_results)

        # Identify best proliferation baseline for this cohort (by AUROC)
        baseline_options = {b: cohort_results[b]["auroc"] for b in ["MKI67", "E2F_TARGETS", "G2M_CHECKPOINT"] if b in cohort_results}
        best_baseline = max(baseline_options, key=baseline_options.get)
        print(f"  Best baseline = {best_baseline}")

        # Paired comparisons: each framework signature vs best baseline
        for drug_name in DRUG_SIGS.keys():
            if drug_name not in scores or best_baseline not in scores:
                continue
            a = scores[drug_name]
            b = scores[best_baseline]
            d_auroc, lo_a, hi_a = paired_bootstrap_delta(y, a, b, roc_auc_score)
            d_prauc, lo_p, hi_p = paired_bootstrap_delta(y, a, b, average_precision_score)
            def neg_brier(yy, ss):
                return -brier_score_loss(yy, score_to_prob(ss))
            d_brier_neg, lo_b, hi_b = paired_bootstrap_delta(y, a, b, neg_brier)

            # DCA differences at key thresholds (model net benefit, framework - baseline)
            prob_a = score_to_prob(a); prob_b = score_to_prob(b); n = len(y)
            dca_diffs = {}
            for t in [0.05, 0.10, 0.15, 0.20, 0.25]:
                pred_a = (prob_a >= t).astype(int); pred_b = (prob_b >= t).astype(int)
                nb_a = ((pred_a == 1) & (y == 1)).sum() / n - ((pred_a == 1) & (y == 0)).sum() / n * (t / (1 - t))
                nb_b = ((pred_b == 1) & (y == 1)).sum() / n - ((pred_b == 1) & (y == 0)).sum() / n * (t / (1 - t))
                dca_diffs[f"t{int(t*100):02d}"] = float(nb_a - nb_b)

            entry = dict(
                cohort=cohort_name, framework_signature=drug_name, baseline=best_baseline,
                framework_auroc=float(roc_auc_score(y, a)),
                baseline_auroc=float(roc_auc_score(y, b)),
                delta_auroc=d_auroc, delta_auroc_ci=[lo_a, hi_a],
                framework_pr_auc=float(average_precision_score(y, a)),
                baseline_pr_auc=float(average_precision_score(y, b)),
                delta_pr_auc=d_prauc, delta_pr_auc_ci=[lo_p, hi_p],
                framework_brier=float(brier_score_loss(y, score_to_prob(a))),
                baseline_brier=float(brier_score_loss(y, score_to_prob(b))),
                delta_neg_brier=d_brier_neg, delta_neg_brier_ci=[lo_b, hi_b],
                dca_nb_diff_at_threshold=dca_diffs,
                pcr_rate=float(y.mean()), n=int(len(y)),
            )
            # Decision: framework wins on this metric if Δ > 0 and CI excludes 0
            entry["wins_auroc"] = bool(d_auroc > 0 and lo_a > 0)
            entry["wins_pr_auc"] = bool(d_prauc > 0 and lo_p > 0)
            entry["wins_brier"] = bool(d_brier_neg > 0 and lo_b > 0)  # negative Brier (lower=better), so positive Δ means framework better
            paired_results.append(entry)

    # Save
    out_dir = REVISE_ROOT / "primary_test_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "multi_metric_comparison_full.json", "w") as f:
        json.dump(all_results, f, indent=2)
    with open(out_dir / "multi_metric_paired.json", "w") as f:
        json.dump(paired_results, f, indent=2)
    pd.DataFrame(paired_results).to_csv(out_dir / "multi_metric_paired.tsv", sep="\t", index=False, float_format="%.4g")

    # Headline summary
    print("\n" + "=" * 90)
    print("HEADLINE: framework signatures vs best proliferation baseline, paired bootstrap deltas")
    print("=" * 90)
    print(f"{'Cohort':<10} {'Drug':<16} {'Baseline':<16} {'ΔAUROC':>10} {'ΔPR-AUC':>10} {'ΔBrier(neg)':>12}  Wins")
    for r in paired_results:
        wins = []
        if r["wins_auroc"]: wins.append("AUROC")
        if r["wins_pr_auc"]: wins.append("PR-AUC")
        if r["wins_brier"]: wins.append("Brier")
        wins_str = "+".join(wins) if wins else "—"
        print(f"{r['cohort']:<10} {r['framework_signature']:<16} {r['baseline']:<16} "
              f"{r['delta_auroc']:>+10.3f} {r['delta_pr_auc']:>+10.3f} {r['delta_neg_brier']:>+12.3f}  {wins_str}")

    # DCA threshold table
    print("\n" + "=" * 90)
    print("DCA: Δ(framework − baseline) net benefit at clinical thresholds")
    print("=" * 90)
    print(f"{'Cohort':<10} {'Drug':<16} {'t=05%':>10} {'t=10%':>10} {'t=15%':>10} {'t=20%':>10} {'t=25%':>10}")
    for r in paired_results:
        d = r["dca_nb_diff_at_threshold"]
        print(f"{r['cohort']:<10} {r['framework_signature']:<16} "
              f"{d['t05']:>+10.4f} {d['t10']:>+10.4f} {d['t15']:>+10.4f} {d['t20']:>+10.4f} {d['t25']:>+10.4f}")

    print(f"\nResults written to {out_dir}/multi_metric_*.{{json,tsv}}")


if __name__ == "__main__":
    main()
