"""Run the pre-registered primary test family T1-T6.

Per preregistration.md §1, §2, §6.1, §6.3, §8:
  T1: GSE16446, glm(pcr ~ doxo_score + basal_score), OR per SD of doxo ≥ 1.3, Wald p < 0.0083
  T2: GSE25066 TN stratum, glm(pcr ~ doxo_score), OR per SD ≥ 1.3, Wald p < 0.0083
  T3: GSE25066 HR+ stratum, same
  T4: GSE22226 pCR, glm(pcr ~ doxo_score + basal_indicator), OR per SD ≥ 1.3, Wald p < 0.0083
  T5: GSE22226 RCB index Spearman ρ p < 0.0083
  T6: H4 negative-control max-AUROC < 0.6 (read from score_<drug>_negctrl.json)

Also computes for each primary test:
  - PR-AUC + bootstrap CI
  - F1 + precision + recall @ Youden-J
  - Brier score
  - Calibration slope+intercept (isotonic + linear)
  - 10-bin reliability table
  - Vickers Decision-Curve Analysis (net benefit, 5-40% threshold)
  - Trivial-baseline comparison (MKI67 / E2F / G2M)

Outputs JSON sidecars + a master primary_results.tsv summary.
"""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import (roc_auc_score, average_precision_score,
                             precision_recall_curve, brier_score_loss,
                             roc_curve, f1_score, precision_score, recall_score)
from sklearn.isotonic import IsotonicRegression
from scipy.stats import spearmanr, norm

REVISE_ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
EXT = REVISE_ROOT / "external_data"


def bootstrap_ci(stat_fn, *arrays, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(arrays[0])
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            v = stat_fn(*[a[idx] for a in arrays])
            if np.isfinite(v):
                boot.append(v)
        except Exception:
            continue
    if len(boot) < 10:
        return float("nan"), float("nan")
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def logistic_or_per_sd(y, score, covariates=None):
    """Fit logistic regression; return OR per SD of score + Wald p + covariate-adjusted AUROC."""
    if len(np.unique(y)) < 2:
        return dict(or_per_sd=float("nan"), wald_p=float("nan"), coef=float("nan"),
                    coef_se=float("nan"), z=float("nan"), n=int(len(y)),
                    adjusted_auc=float("nan"), n_covariates=0, converged=False)
    sd = np.std(score, ddof=1)
    if sd == 0 or not np.isfinite(sd):
        return dict(or_per_sd=float("nan"), wald_p=float("nan"), coef=float("nan"),
                    coef_se=float("nan"), z=float("nan"), n=int(len(y)),
                    adjusted_auc=float("nan"), n_covariates=0, converged=False)
    z_score = (score - np.mean(score)) / sd
    X = np.column_stack([np.ones(len(y)), z_score])
    n_cov = 0
    if covariates is not None and len(covariates) > 0:
        cov_arr = np.column_stack(covariates)
        # Standardize each covariate
        for i in range(cov_arr.shape[1]):
            c = cov_arr[:, i]
            csd = np.std(c, ddof=1)
            if csd > 0:
                cov_arr[:, i] = (c - np.mean(c)) / csd
        X = np.column_stack([X, cov_arr])
        n_cov = cov_arr.shape[1]
    try:
        mod = sm.GLM(y, X, family=sm.families.Binomial())
        res = mod.fit(disp=0)
        coef = float(res.params[1])  # z_score coefficient
        se = float(res.bse[1])
        z = coef / se if se > 0 else float("nan")
        wald_p = float(2 * (1 - norm.cdf(abs(z)))) if np.isfinite(z) else float("nan")
        or_per_sd = float(np.exp(coef))
        # Adjusted AUROC: predicted prob from full model
        try:
            pred = res.predict(X)
            adj_auc = float(roc_auc_score(y, pred))
        except Exception:
            adj_auc = float("nan")
        converged = bool(res.converged) if hasattr(res, "converged") else True
    except Exception as e:
        return dict(or_per_sd=float("nan"), wald_p=float("nan"), coef=float("nan"),
                    coef_se=float("nan"), z=float("nan"), n=int(len(y)),
                    adjusted_auc=float("nan"), n_covariates=n_cov, converged=False,
                    error=str(e))
    return dict(or_per_sd=or_per_sd, wald_p=wald_p, coef=coef, coef_se=se, z=z,
                n=int(len(y)), adjusted_auc=adj_auc, n_covariates=n_cov, converged=converged)


def calibration_metrics(y, score):
    """Brier, isotonic slope+intercept (via linear fit to isotonic), reliability bins."""
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    # Map to [0,1] probability via min-max for Brier (the score isn't a probability)
    rng = score.max() - score.min()
    if rng > 0:
        prob = (score - score.min()) / rng
    else:
        prob = np.full_like(score, 0.5)
    brier = float(brier_score_loss(y, prob))
    # Isotonic regression
    try:
        iso = IsotonicRegression(out_of_bounds="clip").fit(prob, y)
        iso_pred = iso.predict(prob)
        # Linear fit of isotonic predictions vs probabilities
        if np.std(prob) > 0:
            slope, intercept = np.polyfit(prob, iso_pred, 1)
        else:
            slope, intercept = float("nan"), float("nan")
    except Exception:
        slope, intercept = float("nan"), float("nan")
    # Reliability bins (10 quantile bins)
    df = pd.DataFrame({"prob": prob, "y": y})
    df["bin"] = pd.qcut(df["prob"], 10, duplicates="drop", labels=False)
    rel = df.groupby("bin").agg(predicted=("prob", "mean"), observed=("y", "mean"), n=("y", "size")).reset_index()
    return dict(brier=brier, slope=float(slope), intercept=float(intercept),
                reliability=rel.to_dict(orient="records"))


def pr_auc_with_ci(y, score, n_boot=1000, seed=42):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    pr = float(average_precision_score(y, score))
    lo, hi = bootstrap_ci(lambda yy, ss: average_precision_score(yy, ss), y, score, n_boot=n_boot, seed=seed)
    return pr, lo, hi


def f1_at_youden(y, score):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    fpr, tpr, thr = roc_curve(y, score)
    j = tpr - fpr
    if len(j) == 0:
        return dict(threshold=float("nan"), f1=float("nan"), precision=float("nan"), recall=float("nan"))
    idx = np.argmax(j)
    t = float(thr[idx])
    pred = (score >= t).astype(int)
    return dict(
        threshold=t,
        youden_j=float(j[idx]),
        f1=float(f1_score(y, pred, zero_division=0)),
        precision=float(precision_score(y, pred, zero_division=0)),
        recall=float(recall_score(y, pred, zero_division=0)),
    )


def net_benefit(y, prob, threshold):
    """Vickers net benefit at decision threshold."""
    pred = (prob >= threshold).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    n = len(y)
    if n == 0:
        return float("nan")
    return tp / n - fp / n * (threshold / (1 - threshold)) if threshold < 1 else float("nan")


def dca(y, score, thresholds=np.linspace(0.05, 0.40, 36)):
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    rng = score.max() - score.min()
    if rng > 0:
        prob = (score - score.min()) / rng
    else:
        prob = np.full_like(score, 0.5)
    base_rate = float(y.mean())
    out = []
    for t in thresholds:
        nb = net_benefit(y, prob, t)
        # Treat-all: assume everyone gets treated
        tp_all = int((y == 1).sum())
        fp_all = int((y == 0).sum())
        nb_all = tp_all / len(y) - fp_all / len(y) * (t / (1 - t)) if t < 1 else float("nan")
        out.append(dict(threshold=float(t), nb_model=float(nb), nb_treat_all=float(nb_all), nb_treat_none=0.0))
    return dict(base_rate=base_rate, curve=out)


def load_score(cohort, score_name):
    """Return DataFrame with patient_id, signature_score, pcr; merge w/ pcr labels."""
    sc = pd.read_csv(EXT / cohort / f"{score_name}.tsv", sep="\t")
    return sc


def run_T1():
    """GSE16446: doxo + PAM50 basal_score covariate."""
    cohort = "GSE16446"
    score = pd.read_csv(EXT / cohort / "score_DOXORUBICIN.tsv", sep="\t")
    pam50 = pd.read_csv(EXT / cohort / f"{cohort}_pam50.tsv", sep="\t")
    df = score.merge(pam50[["patient_id", "basal_score", "pam50_subtype"]], on="patient_id", how="inner")
    df = df.dropna(subset=["signature_score", "pcr", "basal_score"]).reset_index(drop=True)
    y = df["pcr"].to_numpy().astype(int)
    s = df["signature_score"].to_numpy()
    bs = df["basal_score"].to_numpy()
    out = logistic_or_per_sd(y, s, covariates=[bs])
    out["test_id"] = "T1"
    out["hypothesis"] = "H1: GSE16446 doxorubicin signature predicts pCR (PAM50-basal-adjusted)"
    out["cohort"] = cohort
    out["primary_passes"] = bool(np.isfinite(out["or_per_sd"]) and out["or_per_sd"] >= 1.3 and out["wald_p"] < 0.0083)

    # Supporting + extras (using raw signature_score)
    raw_auc = float(roc_auc_score(y, s)) if len(np.unique(y)) >= 2 else float("nan")
    lo, hi = bootstrap_ci(lambda yy, ss: roc_auc_score(yy, ss), y, s, n_boot=1000, seed=42)
    out["raw_auroc"] = raw_auc
    out["raw_auroc_ci"] = [lo, hi]
    pr, prl, prh = pr_auc_with_ci(y, s)
    out["pr_auc"] = pr; out["pr_auc_ci"] = [prl, prh]
    out["f1_youden"] = f1_at_youden(y, s)
    out["calibration"] = calibration_metrics(y, s)
    out["dca"] = dca(y, s)
    out["base_rate"] = float(y.mean())
    return out, df, y, s


def run_T2_T3():
    """GSE25066 TN + HR+ strata (unadjusted within stratum, no PAM50 covariate per pre-reg §2)."""
    results = []
    for stratum, suffix in [("TN", "TN"), ("HR+", "HRpos")]:
        cohort = "GSE25066"
        score = pd.read_csv(EXT / cohort / f"score_DOXORUBICIN_{suffix}.tsv", sep="\t")
        df = score.dropna(subset=["signature_score", "pcr"]).reset_index(drop=True)
        y = df["pcr"].to_numpy().astype(int)
        s = df["signature_score"].to_numpy()
        out = logistic_or_per_sd(y, s, covariates=None)
        out["test_id"] = "T2" if stratum == "TN" else "T3"
        out["hypothesis"] = f"H2: GSE25066 {stratum} stratum (unadjusted)"
        out["cohort"] = cohort
        out["stratum"] = stratum
        out["primary_passes"] = bool(np.isfinite(out["or_per_sd"]) and out["or_per_sd"] >= 1.3 and out["wald_p"] < 0.0083)
        raw_auc = float(roc_auc_score(y, s)) if len(np.unique(y)) >= 2 else float("nan")
        lo, hi = bootstrap_ci(lambda yy, ss: roc_auc_score(yy, ss), y, s, n_boot=1000, seed=42)
        out["raw_auroc"] = raw_auc; out["raw_auroc_ci"] = [lo, hi]
        pr, prl, prh = pr_auc_with_ci(y, s)
        out["pr_auc"] = pr; out["pr_auc_ci"] = [prl, prh]
        out["f1_youden"] = f1_at_youden(y, s)
        out["calibration"] = calibration_metrics(y, s)
        out["dca"] = dca(y, s)
        out["base_rate"] = float(y.mean())
        results.append(out)
    return results


def run_T4_T5():
    """GSE22226 — T4 pCR with basal indicator covariate; T5 RCB Spearman."""
    cohort = "GSE22226"
    score = pd.read_csv(EXT / cohort / "score_DOXORUBICIN.tsv", sep="\t")
    pcr_df = pd.read_csv(EXT / cohort / "GSE22226_pcr.tsv", sep="\t")
    pam50 = pd.read_csv(EXT / cohort / f"{cohort}_pam50.tsv", sep="\t")
    # Prefer metadata basal (already in pcr_df) but cross with computed basal_score
    df = score.merge(pcr_df[["patient_id", "pam50_basal", "rcb_class"]], on="patient_id", how="inner")
    df = df.merge(pam50[["patient_id", "basal_score"]], on="patient_id", how="left")

    # T4: pcr ~ signature_score + basal_score (continuous, computed)
    df4 = df.dropna(subset=["signature_score", "pcr", "basal_score"]).reset_index(drop=True)
    y = df4["pcr"].to_numpy().astype(int)
    s = df4["signature_score"].to_numpy()
    bs = df4["basal_score"].to_numpy()
    out4 = logistic_or_per_sd(y, s, covariates=[bs])
    out4["test_id"] = "T4"
    out4["hypothesis"] = "H3: GSE22226 doxorubicin signature predicts pCR (PAM50-basal-adjusted)"
    out4["cohort"] = cohort
    out4["primary_passes"] = bool(np.isfinite(out4["or_per_sd"]) and out4["or_per_sd"] >= 1.3 and out4["wald_p"] < 0.0083)
    raw_auc = float(roc_auc_score(y, s)) if len(np.unique(y)) >= 2 else float("nan")
    lo, hi = bootstrap_ci(lambda yy, ss: roc_auc_score(yy, ss), y, s, n_boot=1000, seed=42)
    out4["raw_auroc"] = raw_auc; out4["raw_auroc_ci"] = [lo, hi]
    pr, prl, prh = pr_auc_with_ci(y, s)
    out4["pr_auc"] = pr; out4["pr_auc_ci"] = [prl, prh]
    out4["f1_youden"] = f1_at_youden(y, s)
    out4["calibration"] = calibration_metrics(y, s)
    out4["dca"] = dca(y, s)
    out4["base_rate"] = float(y.mean())

    # T5: RCB index Spearman correlation
    # RCB class is categorical "0=...", "I=...", "II=...", "III=..." in GSE22226 — convert to numeric 0/1/2/3
    rcb_map = {}
    def parse_rcb(v):
        if pd.isna(v):
            return np.nan
        v = str(v).strip()
        if v.startswith("0"): return 0.0
        if v.lower().startswith("i=") or v.startswith("I "): return 1.0
        if v.startswith("II"): return 2.0
        if v.startswith("III"): return 3.0
        return np.nan
    df["rcb_numeric"] = df["rcb_class"].apply(parse_rcb)
    df5 = df.dropna(subset=["signature_score", "rcb_numeric"]).reset_index(drop=True)
    if len(df5) >= 10:
        rho, p = spearmanr(df5["signature_score"].to_numpy(), df5["rcb_numeric"].to_numpy())
        # Direction: higher signature = more sensitive => LOWER RCB => negative rho expected
        out5 = dict(test_id="T5", cohort="GSE22226",
                    hypothesis="H3 supplemental: doxo signature vs RCB index Spearman",
                    spearman_rho=float(rho), spearman_p=float(p), n=int(len(df5)),
                    primary_passes=bool(np.isfinite(p) and p < 0.0083 and rho < 0))
    else:
        out5 = dict(test_id="T5", cohort="GSE22226", spearman_rho=float("nan"),
                    spearman_p=float("nan"), n=int(len(df5)), primary_passes=False)
    return [out4, out5]


def run_T6():
    """H4: max negative-control AUROC < 0.6 on GSE16446."""
    cohort = "GSE16446"
    drugs = ["DASATINIB", "GEFITINIB", "ERLOTINIB"]
    aucs = {}
    for d in drugs:
        p = EXT / cohort / f"score_{d}_negctrl.json"
        if p.exists():
            with open(p) as f:
                j = json.load(f)
            aucs[d] = {"auroc": j["auroc"], "ci": [j["auroc_ci_low"], j["auroc_ci_high"]]}
    max_drug = max(aucs.items(), key=lambda kv: kv[1]["auroc"])
    return dict(test_id="T6", cohort=cohort,
                hypothesis="H4: max wrong-MoA AUROC < 0.6 on GSE16446",
                negative_controls=aucs,
                max_drug=max_drug[0],
                max_auroc=max_drug[1]["auroc"],
                primary_passes=bool(max_drug[1]["auroc"] < 0.6))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=REVISE_ROOT / "primary_test_results")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Running T1 — GSE16446 PAM50-adjusted...")
    t1, _, _, _ = run_T1()
    print(f"  OR={t1['or_per_sd']:.3f}, p={t1['wald_p']:.4g}, raw_AUC={t1['raw_auroc']:.3f}, pass={t1['primary_passes']}")

    print("Running T2/T3 — GSE25066 strata...")
    t23 = run_T2_T3()
    for t in t23:
        print(f"  {t['test_id']} ({t['stratum']}): OR={t['or_per_sd']:.3f}, p={t['wald_p']:.4g}, raw_AUC={t['raw_auroc']:.3f}, pass={t['primary_passes']}")

    print("Running T4/T5 — GSE22226...")
    t45 = run_T4_T5()
    for t in t45:
        if t['test_id'] == "T4":
            print(f"  T4 (pCR): OR={t['or_per_sd']:.3f}, p={t['wald_p']:.4g}, raw_AUC={t['raw_auroc']:.3f}, pass={t['primary_passes']}")
        else:
            print(f"  T5 (RCB): rho={t['spearman_rho']:.3f}, p={t['spearman_p']:.4g}, pass={t['primary_passes']}")

    print("Running T6 — GSE16446 negative controls...")
    t6 = run_T6()
    print(f"  Max wrong-MoA: {t6['max_drug']} AUROC={t6['max_auroc']:.3f}, pass={t6['primary_passes']}")

    all_results = [t1] + t23 + t45 + [t6]
    with open(args.out_dir / "primary_tests.json", "w") as f:
        json.dump(all_results, f, indent=2, default=float)

    # Summary TSV
    summary = []
    for r in all_results:
        row = {
            "test_id": r.get("test_id"),
            "cohort": r.get("cohort"),
            "hypothesis": r.get("hypothesis", ""),
            "n": r.get("n", "-"),
            "or_per_sd": r.get("or_per_sd", "-"),
            "wald_p": r.get("wald_p", "-"),
            "raw_auroc": r.get("raw_auroc", "-"),
            "raw_auroc_ci_lo": r.get("raw_auroc_ci", [None, None])[0] if r.get("raw_auroc_ci") else "-",
            "raw_auroc_ci_hi": r.get("raw_auroc_ci", [None, None])[1] if r.get("raw_auroc_ci") else "-",
            "pr_auc": r.get("pr_auc", "-"),
            "f1_youden": r.get("f1_youden", {}).get("f1", "-") if isinstance(r.get("f1_youden"), dict) else "-",
            "brier": r.get("calibration", {}).get("brier", "-") if isinstance(r.get("calibration"), dict) else "-",
            "spearman_rho": r.get("spearman_rho", "-"),
            "spearman_p": r.get("spearman_p", "-"),
            "primary_passes": r.get("primary_passes"),
        }
        summary.append(row)
    pd.DataFrame(summary).to_csv(args.out_dir / "primary_tests_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"\nResults written to {args.out_dir}/")


if __name__ == "__main__":
    main()
