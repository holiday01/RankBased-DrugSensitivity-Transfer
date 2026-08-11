"""Exploratory: PACLITAXEL signature OR-based test on GSE22226 (AC-T cohort).

EXPLORATORY status (per preregistration.md §6.3, §10):
  This is NOT in the locked PRIMARY test family (T1-T6). Reported as
  secondary mechanism evidence supporting the H3b exploratory hypothesis
  that the taxane signature outperforms the doxorubicin signature on
  AC-T cohorts.

Same machinery as primary_tests.py T4:
  glm(pcr ~ paclitaxel_score + pam50_basal_score, family=binomial)
  OR per SD + Wald p + adjusted AUROC + bootstrap CI + PR-AUC + calibration + DCA

Also runs Doxorubicin head-to-head on the same patients (paired bootstrap ΔAUROC).
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from scipy.stats import norm, spearmanr

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
    if len(np.unique(y)) < 2:
        return dict(or_per_sd=float("nan"), wald_p=float("nan"), coef=float("nan"),
                    coef_se=float("nan"), z=float("nan"), n=int(len(y)),
                    adjusted_auc=float("nan"), n_covariates=0)
    sd = np.std(score, ddof=1)
    if sd == 0 or not np.isfinite(sd):
        return dict(or_per_sd=float("nan"), wald_p=float("nan"), n=int(len(y)),
                    adjusted_auc=float("nan"), n_covariates=0)
    z_score = (score - np.mean(score)) / sd
    X = np.column_stack([np.ones(len(y)), z_score])
    n_cov = 0
    if covariates is not None and len(covariates) > 0:
        cov_arr = np.column_stack(covariates).astype(float)
        for i in range(cov_arr.shape[1]):
            c = cov_arr[:, i]
            csd = np.std(c, ddof=1)
            if csd > 0:
                cov_arr[:, i] = (c - np.mean(c)) / csd
        X = np.column_stack([X, cov_arr])
        n_cov = cov_arr.shape[1]
    mod = sm.GLM(y, X, family=sm.families.Binomial())
    res = mod.fit(disp=0)
    coef = float(res.params[1])
    se = float(res.bse[1])
    z = coef / se if se > 0 else float("nan")
    wald_p = float(2 * (1 - norm.cdf(abs(z)))) if np.isfinite(z) else float("nan")
    or_per_sd = float(np.exp(coef))
    pred = res.predict(X)
    adj_auc = float(roc_auc_score(y, pred))
    return dict(or_per_sd=or_per_sd, wald_p=wald_p, coef=coef, coef_se=se, z=z,
                n=int(len(y)), adjusted_auc=adj_auc, n_covariates=n_cov,
                converged=bool(res.converged) if hasattr(res, "converged") else True,
                covariate_coefs=[float(c) for c in res.params[2:]],
                covariate_pvals=[float(2 * (1 - norm.cdf(abs(p / se_)))) for p, se_ in zip(res.params[2:], res.bse[2:])])


def main():
    out_dir = REVISE_ROOT / "primary_test_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort = "GSE22226"

    # Load Paclitaxel score
    pac = pd.read_csv(EXT / cohort / "score_PACLITAXEL.tsv", sep="\t")
    pac = pac.rename(columns={"signature_score": "paclitaxel_score"})

    # Load Doxorubicin score (paired comparison)
    dox = pd.read_csv(EXT / cohort / "score_DOXORUBICIN.tsv", sep="\t")
    dox = dox.rename(columns={"signature_score": "doxorubicin_score"})

    # Load PAM50 basal score
    pam = pd.read_csv(EXT / cohort / f"{cohort}_pam50.tsv", sep="\t")[["patient_id", "basal_score"]]

    # Merge
    df = pac.merge(dox[["patient_id", "doxorubicin_score"]], on="patient_id", how="inner")
    df = df.merge(pam, on="patient_id", how="inner").dropna(
        subset=["paclitaxel_score", "doxorubicin_score", "basal_score", "pcr"]
    ).reset_index(drop=True)

    y = df["pcr"].astype(int).to_numpy()
    pac_s = df["paclitaxel_score"].to_numpy()
    dox_s = df["doxorubicin_score"].to_numpy()
    bs = df["basal_score"].to_numpy()

    print(f"n={len(df)}, pcr_rate={y.mean():.3f}, n_responders={int(y.sum())}")

    # PRIMARY-style test for Paclitaxel: glm(pcr ~ pac_score + basal_score)
    print("\n=== EXPLORATORY: PACLITAXEL OR test (PAM50-basal-adjusted) ===")
    pac_or = logistic_or_per_sd(y, pac_s, covariates=[bs])
    print(f"  OR per SD: {pac_or['or_per_sd']:.3f}")
    print(f"  Wald p:    {pac_or['wald_p']:.4g}")
    print(f"  Adjusted AUROC: {pac_or['adjusted_auc']:.3f}")
    print(f"  Pass PRIMARY threshold (OR ≥ 1.3 AND p < 0.0083)? {pac_or['or_per_sd'] >= 1.3 and pac_or['wald_p'] < 0.0083}")
    print(f"  basal_score covariate coef: {pac_or['covariate_coefs'][0]:.3f} (p={pac_or['covariate_pvals'][0]:.4g})")

    # Unadjusted Paclitaxel for comparison
    pac_or_unadj = logistic_or_per_sd(y, pac_s, covariates=None)
    print(f"\n=== EXPLORATORY: PACLITAXEL OR test (UNADJUSTED) ===")
    print(f"  OR per SD: {pac_or_unadj['or_per_sd']:.3f}, Wald p: {pac_or_unadj['wald_p']:.4g}")
    print(f"  AUROC: {pac_or_unadj['adjusted_auc']:.3f}")

    # Doxorubicin head-to-head (already in T4 but recompute on EXACTLY the same patients)
    print("\n=== HEAD-TO-HEAD: same n={n} patients ===".format(n=len(df)))
    dox_or = logistic_or_per_sd(y, dox_s, covariates=[bs])
    pac_auc = float(roc_auc_score(y, pac_s))
    dox_auc = float(roc_auc_score(y, dox_s))
    pac_ci = bootstrap_ci(lambda yy, ss: roc_auc_score(yy, ss), y, pac_s, n_boot=1000, seed=42)
    dox_ci = bootstrap_ci(lambda yy, ss: roc_auc_score(yy, ss), y, dox_s, n_boot=1000, seed=42)
    print(f"  Paclitaxel raw AUROC: {pac_auc:.3f} [{pac_ci[0]:.3f}, {pac_ci[1]:.3f}]")
    print(f"  Doxorubicin raw AUROC: {dox_auc:.3f} [{dox_ci[0]:.3f}, {dox_ci[1]:.3f}]")

    # Paired bootstrap ΔAUROC (H3b)
    rng = np.random.default_rng(42)
    deltas = []
    for _ in range(1000):
        idx = rng.integers(0, len(y), size=len(y))
        try:
            d = roc_auc_score(y[idx], pac_s[idx]) - roc_auc_score(y[idx], dox_s[idx])
            if np.isfinite(d):
                deltas.append(d)
        except Exception:
            continue
    delta_observed = pac_auc - dox_auc
    delta_ci = [float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))]
    print(f"  Paired Δ-AUROC (Pac − Dox): observed {delta_observed:+.3f}, 95% CI [{delta_ci[0]:+.3f}, {delta_ci[1]:+.3f}]")
    print(f"  Pre-reg H3b threshold (Δ ≥ +0.05 expected if taxane signature is the better proxy): {'CLEARED' if delta_observed >= 0.05 else 'NOT CLEARED'} (CI lower bound {'> 0' if delta_ci[0] > 0 else '<= 0'})")

    # Logistic regression: pcr ~ pac + dox + basal (which signature has independent contribution?)
    print("\n=== JOINT MODEL: pcr ~ Pac_z + Dox_z + Basal_z ===")
    pac_z = (pac_s - pac_s.mean()) / pac_s.std(ddof=1)
    dox_z = (dox_s - dox_s.mean()) / dox_s.std(ddof=1)
    bs_z = (bs - bs.mean()) / bs.std(ddof=1)
    X = np.column_stack([np.ones(len(y)), pac_z, dox_z, bs_z])
    mod = sm.GLM(y, X, family=sm.families.Binomial())
    jr = mod.fit(disp=0)
    print(f"  Intercept:   {jr.params[0]:+.3f}")
    print(f"  Pac coef:    {jr.params[1]:+.3f} (OR/SD={np.exp(jr.params[1]):.3f}, p={jr.pvalues[1]:.4g})")
    print(f"  Dox coef:    {jr.params[2]:+.3f} (OR/SD={np.exp(jr.params[2]):.3f}, p={jr.pvalues[2]:.4g})")
    print(f"  Basal coef:  {jr.params[3]:+.3f} (OR/SD={np.exp(jr.params[3]):.3f}, p={jr.pvalues[3]:.4g})")
    print(f"  Model joint AUROC: {roc_auc_score(y, jr.predict(X)):.3f}")

    # Save
    out = {
        "test_id": "EXPLORATORY_paclitaxel_GSE22226",
        "preregistration_status": "exploratory (not in alpha=0.0083 primary family per preregistration §6.3, §10)",
        "cohort": cohort,
        "regimen": "AC-T (Doxorubicin+Cyclophosphamide → Paclitaxel)",
        "n_patients": int(len(df)),
        "pcr_rate": float(y.mean()),
        "paclitaxel_adjusted": {
            "or_per_sd": pac_or["or_per_sd"],
            "wald_p": pac_or["wald_p"],
            "adjusted_auroc": pac_or["adjusted_auc"],
            "covariate_basal_coef": pac_or["covariate_coefs"][0],
            "covariate_basal_p": pac_or["covariate_pvals"][0],
            "passes_primary_threshold": bool(pac_or["or_per_sd"] >= 1.3 and pac_or["wald_p"] < 0.0083),
        },
        "paclitaxel_unadjusted": {
            "or_per_sd": pac_or_unadj["or_per_sd"],
            "wald_p": pac_or_unadj["wald_p"],
            "auroc": pac_or_unadj["adjusted_auc"],
        },
        "paclitaxel_raw_auroc": pac_auc,
        "paclitaxel_raw_auroc_ci": list(pac_ci),
        "doxorubicin_raw_auroc": dox_auc,
        "doxorubicin_raw_auroc_ci": list(dox_ci),
        "paired_delta_auroc_pac_minus_dox": float(delta_observed),
        "paired_delta_auroc_ci_95": delta_ci,
        "h3b_threshold_cleared": bool(delta_observed >= 0.05),
        "h3b_ci_above_zero": bool(delta_ci[0] > 0),
        "joint_model": {
            "pac_coef": float(jr.params[1]),
            "pac_or_per_sd": float(np.exp(jr.params[1])),
            "pac_p": float(jr.pvalues[1]),
            "dox_coef": float(jr.params[2]),
            "dox_or_per_sd": float(np.exp(jr.params[2])),
            "dox_p": float(jr.pvalues[2]),
            "basal_coef": float(jr.params[3]),
            "basal_p": float(jr.pvalues[3]),
            "joint_auroc": float(roc_auc_score(y, jr.predict(X))),
        },
    }
    with open(out_dir / "exploratory_paclitaxel_GSE22226.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults written to {out_dir}/exploratory_paclitaxel_GSE22226.json")


if __name__ == "__main__":
    main()
