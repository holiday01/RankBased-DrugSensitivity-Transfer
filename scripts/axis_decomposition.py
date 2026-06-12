"""Phase 1 §4.1 — Axis decomposition of the frozen Doxorubicin signature.

For each cohort:
  1. Fit OLS:  sig_dox ~ b0 + b1*P + b2*H + b3*B + b4*H_PAM50basal
  2. Compute R^2, partial-R^2 per predictor, VIF
  3. Ridge sensitivity: 10-fold CV alpha (seed 42); coefficient sign stability
  4. Cluster bootstrap 95% CI on R^2 (1000 reps, seed 42)
  5. Residual signature = OLS residual (in-cohort)
  6. Residual AUROC vs pCR + bootstrap 95% CI
  7. Synthetic positive control: sig_synth = alpha*P + beta*H + 0.3*TOP2A_z + 0.5*N(0,1)
     — residual AUROC of sig_synth should exceed ~0.6 (sanity)
  8. Pre-specified falsifier: residual AUROC 95% CI lower bound > 0.55 in >=2/3 cohorts
     => H-b REJECTED

Outputs:
  results/phase_1/axis_decomposition.json
  results/phase_1/axis_decomposition_summary.tsv
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

REPO = Path(".")
PHASE0 = REPO / "." / "results" / "phase_0"
OUT_DIR = REPO / "." / "results" / "phase_1"
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_dox_score(cohort: str) -> pd.DataFrame:
    return pd.read_csv(REPO / "external_data" / cohort / "score_DOXORUBICIN.tsv", sep="\t")


def load_axes(cohort: str) -> pd.DataFrame:
    return pd.read_csv(PHASE0 / f"axes_{cohort}.tsv", sep="\t")


def load_expression_topa(cohort: str) -> pd.Series:
    expr = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_expression.tsv",
                       sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    if "TOP2A" not in expr.columns:
        return pd.Series(dtype=float)
    z = expr["TOP2A"]
    z = (z - z.mean()) / (z.std(ddof=0) + 1e-12)
    return z


def partial_r2(model, X, y, predictors):
    """Partial-R^2 = (R^2_full - R^2_reduced) / (1 - R^2_reduced)."""
    full_r2 = model.rsquared
    out = {}
    for p in predictors:
        keep = [c for c in X.columns if c != p]
        Xr = sm.add_constant(X[keep])
        m_r = sm.OLS(y, Xr).fit()
        r2_r = m_r.rsquared
        pr2 = (full_r2 - r2_r) / (1.0 - r2_r) if (1.0 - r2_r) > 1e-12 else float("nan")
        out[p] = float(pr2)
    return out


def vif(X):
    """Variance inflation factor for each column."""
    from statsmodels.stats.outliers_influence import variance_inflation_factor as _vif
    Xc = sm.add_constant(X).values
    cols = ["const"] + list(X.columns)
    out = {}
    for i, c in enumerate(cols):
        if c == "const":
            continue
        try:
            v = float(_vif(Xc, i))
        except Exception:
            v = float("nan")
        out[c] = v
    return out


def boot_auroc(y, s, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    obs = float(roc_auc_score(y, s)) if len(np.unique(y)) >= 2 else float("nan")
    if not np.isfinite(obs):
        return obs, float("nan"), float("nan")
    boot = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        boot.append(roc_auc_score(y[idx], s[idx]))
    if len(boot) < 10:
        return obs, float("nan"), float("nan")
    return obs, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def cluster_bootstrap_r2(X, y, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(y)
    Xc = sm.add_constant(X)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            m = sm.OLS(y.iloc[idx].values, Xc.iloc[idx].values).fit()
            boot.append(m.rsquared)
        except Exception:
            continue
    if len(boot) < 10:
        return float("nan"), float("nan")
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def ridge_sensitivity(X, y, k=10, seed=42):
    rng = np.random.default_rng(seed)
    kf = KFold(n_splits=k, shuffle=True, random_state=seed)
    coefs = []
    alphas = []
    for tr, _ in kf.split(X):
        Xt = X.iloc[tr].values
        yt = y.iloc[tr].values
        m = RidgeCV(alphas=np.logspace(-3, 3, 25)).fit(Xt, yt)
        coefs.append(m.coef_)
        alphas.append(m.alpha_)
    coefs = np.array(coefs)
    sign = np.sign(coefs)
    # sign stability: fraction of folds where sign matches modal sign (per predictor)
    stability = {}
    for j, col in enumerate(X.columns):
        modal = float(np.sign(np.median(coefs[:, j])))
        if modal == 0:
            modal = 1.0
        stability[col] = float(np.mean(sign[:, j] == modal))
    return {
        "median_coefs": {c: float(np.median(coefs[:, j])) for j, c in enumerate(X.columns)},
        "sign_stability_fraction": stability,
        "median_alpha": float(np.median(alphas)),
    }


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {
        "script": str(me),
        "script_sha256": sha256_of(me),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "cohorts": {},
        "pre_reg_falsifier": (
            "residual AUROC 95% CI lower bound > 0.55 in >= 2/3 cohorts => H-b REJECTED"
        ),
    }
    summary_rows = []
    pass_count = 0

    for cohort in COHORTS:
        axes = load_axes(cohort)
        dox = load_dox_score(cohort)
        # join sig_dox + axes
        df = dox.merge(axes, left_on="patient_id", right_on="sample_id", how="inner")
        # Predictors: P, H, B, H_PAM50basal
        predictors = ["P", "H", "B", "H_PAM50basal"]
        df = df.dropna(subset=["signature_score"] + predictors)
        X = df[predictors].astype(float)
        y_sig = df["signature_score"].astype(float)

        # OLS
        m = sm.OLS(y_sig.values, sm.add_constant(X).values).fit()
        coefs = dict(zip(["const"] + predictors, [float(x) for x in m.params]))
        pvals = dict(zip(["const"] + predictors, [float(x) for x in m.pvalues]))
        r2 = float(m.rsquared)
        ci_lo, ci_hi = cluster_bootstrap_r2(X, y_sig, 1000, 42)
        pr2 = partial_r2(m, X, y_sig, predictors)
        vifs = vif(X)
        ridge = ridge_sensitivity(X, y_sig, k=10, seed=42)

        # Residual signature
        resid = y_sig.values - m.fittedvalues
        df_res = df.copy()
        df_res["sig_dox_residual"] = resid

        # Residual vs pCR
        # pCR may be missing (e.g., GSE16446 axes file shows blanks for some)
        valid = df_res.dropna(subset=["pCR"]).copy()
        valid["pCR"] = valid["pCR"].astype(int)
        if len(np.unique(valid["pCR"])) < 2 or len(valid) < 10:
            res_auc, res_lo, res_hi = float("nan"), float("nan"), float("nan")
        else:
            res_auc, res_lo, res_hi = boot_auroc(valid["pCR"].values,
                                                  valid["sig_dox_residual"].values, 1000, 42)

        # Synthetic positive control
        top2a = load_expression_topa(cohort)
        rng = np.random.default_rng(2026)
        # Build synth on the same patient set as df
        idx_int = df_res["patient_id"].astype(str).tolist()
        top2a_aligned = top2a.reindex(idx_int).astype(float)
        # If TOP2A missing, fill with 0
        top2a_aligned = top2a_aligned.fillna(0.0).values
        noise = rng.normal(0, 1, size=len(df_res))
        sig_synth = (0.7 * df_res["P"].values + 0.5 * df_res["H"].values
                     + 0.3 * top2a_aligned + 0.5 * noise)
        # Regress synth on same predictors, get residuals
        m_synth = sm.OLS(sig_synth, sm.add_constant(X).values).fit()
        synth_resid = sig_synth - m_synth.fittedvalues
        df_res["sig_synth_residual"] = synth_resid
        valid_s = df_res.dropna(subset=["pCR"]).copy()
        valid_s["pCR"] = valid_s["pCR"].astype(int)
        if len(np.unique(valid_s["pCR"])) < 2:
            synth_auc, synth_lo, synth_hi = float("nan"), float("nan"), float("nan")
        else:
            synth_auc, synth_lo, synth_hi = boot_auroc(valid_s["pCR"].values,
                                                       valid_s["sig_synth_residual"].values, 1000, 42)

        # Falsifier check
        residual_pass = (np.isfinite(res_lo) and res_lo > 0.55)
        if residual_pass:
            pass_count += 1

        cohort_block = {
            "n_samples_in_regression": int(len(df_res)),
            "n_samples_with_pCR": int(len(valid)),
            "n_pCR_positive": int(valid["pCR"].sum()) if len(valid) else 0,
            "OLS": {
                "R2": r2,
                "R2_boot_ci_95": [ci_lo, ci_hi],
                "coefficients": coefs,
                "p_values": pvals,
                "partial_R2": pr2,
                "VIF": vifs,
            },
            "ridge_sensitivity": ridge,
            "residual_pcr_AUROC": res_auc,
            "residual_pcr_AUROC_ci_95": [res_lo, res_hi],
            "residual_falsifier_pass": bool(residual_pass),
            "synthetic_positive_control": {
                "construction": "0.7*P + 0.5*H + 0.3*TOP2A_z + 0.5*N(0,1)",
                "residual_AUROC": synth_auc,
                "residual_AUROC_ci_95": [synth_lo, synth_hi],
                "sanity_check_pass": bool(np.isfinite(synth_auc) and synth_auc > 0.55),
            },
        }
        log["cohorts"][cohort] = cohort_block
        summary_rows.append({
            "cohort": cohort,
            "R2_full_model": r2,
            "R2_ci_low": ci_lo, "R2_ci_high": ci_hi,
            "partial_R2_P": pr2["P"],
            "partial_R2_H": pr2["H"],
            "partial_R2_B": pr2["B"],
            "partial_R2_H_PAM50basal": pr2["H_PAM50basal"],
            "VIF_P": vifs["P"], "VIF_H": vifs["H"], "VIF_B": vifs["B"],
            "VIF_H_PAM50basal": vifs["H_PAM50basal"],
            "residual_AUROC": res_auc,
            "residual_AUROC_ci_low": res_lo,
            "residual_AUROC_ci_high": res_hi,
            "synth_residual_AUROC": synth_auc,
            "residual_falsifier_pass": int(residual_pass),
        })

    # Aggregate falsifier verdict
    log["falsifier_pass_count"] = int(pass_count)
    log["H_b_REJECTED"] = bool(pass_count >= 2)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(OUT_DIR / "axis_decomposition.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    pd.DataFrame(summary_rows).to_csv(
        OUT_DIR / "axis_decomposition_summary.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"A1 done. H_b rejected = {log['H_b_REJECTED']} (pass={pass_count}/3)")
    for c in COHORTS:
        b = log["cohorts"][c]
        print(f"  {c}: R2={b['OLS']['R2']:.3f} | residual AUROC={b['residual_pcr_AUROC']:.3f}"
              f" CI=[{b['residual_pcr_AUROC_ci_95'][0]:.3f},{b['residual_pcr_AUROC_ci_95'][1]:.3f}]"
              f" | synth AUROC={b['synthetic_positive_control']['residual_AUROC']:.3f}")


if __name__ == "__main__":
    main()
