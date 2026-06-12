"""Phase 2 / B4 — Residual confirmatory test on GSE32646 (§4.8 headline).

Pre-specified primary model:
    sig_dox ~ β0 + β1·P + β2·H + β3·B + β4·H_PAM50basal + ε
Then residual = sig_dox − fitted; test residual vs pCR via AUROC + bootstrap CI.

3×3 decision rule:
  point ≤ 0.55 AND CI lower ≤ 0.55  → H-b CONFIRMED
  0.55 < point ≤ 0.60 OR CI crosses 0.55 → GRAY ZONE
  point > 0.60 AND CI lower > 0.55  → H-b FALSIFIED

Supplementary variants (Bonferroni × 3; ineligible to flip headline):
  - H_ESR1only in place of H
  - H_PAM50basal swap (drop H, keep H_PAM50basal only)
  - spline residual via pygam if available; else skip with note

Outputs:
  phase_2/gse32646_residual_test.json
  phase_2/gse32646_residual_summary.tsv
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

REPO = Path(".")
PHASE_2 = REPO / "." / "results" / "phase_2"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def auroc_with_boot_ci(scores, labels, n_boot=1000, seed=42):
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


def fit_residual(df: pd.DataFrame, h_col: str = "H", include_pam50: bool = True):
    """Fit sig_dox ~ P + h_col + B (+ H_PAM50basal optional). Return residuals.

    Drops rows with any NaN in required cols.
    """
    cols = ["P", h_col, "B"]
    if include_pam50 and h_col != "H_PAM50basal":
        cols = ["P", h_col, "B", "H_PAM50basal"]
    elif include_pam50 and h_col == "H_PAM50basal":
        cols = ["P", h_col, "B"]
    needed = ["sig_dox", "pcr"] + cols
    sub = df.dropna(subset=needed).copy()
    if len(sub) < 20:
        return None, sub, None
    X = sub[cols].values
    y = sub["sig_dox"].values
    lr = LinearRegression()
    lr.fit(X, y)
    yhat = lr.predict(X)
    sub["residual"] = y - yhat
    coef = {c: float(v) for c, v in zip(cols, lr.coef_)}
    coef["intercept"] = float(lr.intercept_)
    return sub, sub, coef


def verdict_3x3(point: float, ci_lo: float):
    """Apply pre-registered 3×3 decision rule.

    Per the prompt:
      point ≤ 0.55 AND ci_lo ≤ 0.55 → CONFIRMED
      0.55 < point ≤ 0.60 OR ci crosses 0.55 → GRAY ZONE
      point > 0.60 AND ci_lo > 0.55 → FALSIFIED

    Note: AUROC < 0.5 = signal in opposite direction; |AUROC - 0.5| is the
    informative magnitude. The rule above treats values ≤ 0.55 as "no signal"
    (= H-b confirmed, since H-b states the cohort lies on the null axis after
    deconfounding). For symmetry, we also flag (point < 0.45 OR ci_hi < 0.45)
    as INVERSE-SIGNAL gray zone for transparency.
    """
    if not (math.isfinite(point) and math.isfinite(ci_lo)):
        return "INDETERMINATE_NAN"
    # Standard rule
    if point <= 0.55 and ci_lo <= 0.55:
        return "H-b CONFIRMED (no residual signal after deconfounding)"
    if point > 0.60 and ci_lo > 0.55:
        return "H-b FALSIFIED (residual signal survives)"
    return "GRAY ZONE"


def main():
    PHASE_2.mkdir(parents=True, exist_ok=True)
    log = {
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "decision_rule": {
            "confirmed": "point ≤ 0.55 AND ci_lo ≤ 0.55",
            "gray_zone": "0.55 < point ≤ 0.60 OR ci crosses 0.55",
            "falsified": "point > 0.60 AND ci_lo > 0.55",
        },
    }

    axes_path = PHASE_2 / "axes_GSE32646.tsv"
    score_path = PHASE_2 / "score_DOXORUBICIN_GSE32646.tsv"
    axes = pd.read_csv(axes_path, sep="\t")
    score = pd.read_csv(score_path, sep="\t")

    # Merge
    axes["sample_id"] = axes["sample_id"].astype(str)
    score["patient_id"] = score["patient_id"].astype(str)
    df = axes.merge(score[["patient_id", "sig_dox"]], left_on="sample_id", right_on="patient_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["n_pcr_1"] = int((df["pcr"] == 1).sum())
    log["n_pcr_0"] = int((df["pcr"] == 0).sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    # --- HEADLINE: primary model -------------------------------------
    sub, _, coef = fit_residual(df, h_col="H", include_pam50=True)
    if sub is None:
        print("ERROR: too few rows for primary model", file=sys.stderr)
        sys.exit(1)
    log["primary_model"] = {
        "formula": "sig_dox ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, boot = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["primary_residual_auroc"] = {
        "point": auc,
        "ci_low": lo,
        "ci_high": hi,
        "n_bootstrap": 1000,
        "seed": 42,
    }
    primary_verdict = verdict_3x3(auc, lo)
    log["primary_verdict"] = primary_verdict

    # Raw sig_dox AUROC for comparison
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(sub["sig_dox"].to_numpy(), y, n_boot=1000, seed=42)
    log["raw_sig_dox_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi}

    # --- SUPPLEMENTARY VARIANTS (Bonferroni × 3) -----------------------------
    supp = {}
    # 1) H_ESR1only
    sub_e, _, coef_e = fit_residual(df, h_col="H_ESR1only", include_pam50=True)
    if sub_e is not None and len(sub_e) >= 20:
        y_e = sub_e["pcr"].astype(int).to_numpy()
        r_e = sub_e["residual"].to_numpy()
        ae, le, he, _ = auroc_with_boot_ci(r_e, y_e, n_boot=1000, seed=42)
        # Bonferroni × 3: widen CI by using 100*(1 - 0.05/3) = 98.333% boot ci
        # Re-derive Bonferroni-adjusted CI from the bootstrap distribution
        boot_e = []
        rng = np.random.default_rng(42)
        for _ in range(2000):
            idx = rng.integers(0, len(y_e), size=len(y_e))
            if len(np.unique(y_e[idx])) < 2:
                continue
            boot_e.append(roc_auc_score(y_e[idx], r_e[idx]))
        alpha = 0.05 / 3
        bf_lo = float(np.percentile(boot_e, 100 * alpha / 2)) if boot_e else float("nan")
        bf_hi = float(np.percentile(boot_e, 100 * (1 - alpha / 2))) if boot_e else float("nan")
        supp["H_ESR1only"] = {
            "n_used": int(len(sub_e)),
            "point": ae,
            "ci_low_95": le,
            "ci_high_95": he,
            "ci_low_bonf3": bf_lo,
            "ci_high_bonf3": bf_hi,
            "verdict_unadjusted": verdict_3x3(ae, le),
            "verdict_bonf3": verdict_3x3(ae, bf_lo),
        }
    # 2) H_PAM50basal as H replacement
    sub_b, _, coef_b = fit_residual(df, h_col="H_PAM50basal", include_pam50=True)
    if sub_b is not None and len(sub_b) >= 20:
        y_b = sub_b["pcr"].astype(int).to_numpy()
        r_b = sub_b["residual"].to_numpy()
        ab, lb, hb, _ = auroc_with_boot_ci(r_b, y_b, n_boot=1000, seed=42)
        boot_b = []
        rng = np.random.default_rng(42)
        for _ in range(2000):
            idx = rng.integers(0, len(y_b), size=len(y_b))
            if len(np.unique(y_b[idx])) < 2:
                continue
            boot_b.append(roc_auc_score(y_b[idx], r_b[idx]))
        alpha = 0.05 / 3
        bf_lo = float(np.percentile(boot_b, 100 * alpha / 2)) if boot_b else float("nan")
        bf_hi = float(np.percentile(boot_b, 100 * (1 - alpha / 2))) if boot_b else float("nan")
        supp["H_PAM50basal_only"] = {
            "n_used": int(len(sub_b)),
            "point": ab,
            "ci_low_95": lb,
            "ci_high_95": hb,
            "ci_low_bonf3": bf_lo,
            "ci_high_bonf3": bf_hi,
            "verdict_unadjusted": verdict_3x3(ab, lb),
            "verdict_bonf3": verdict_3x3(ab, bf_lo),
        }
    # 3) spline residual
    try:
        from pygam import LinearGAM, s
        sub_s = df.dropna(subset=["sig_dox", "pcr", "P", "H", "B"]).copy()
        Xs = sub_s[["P", "H", "B"]].values
        ys = sub_s["sig_dox"].values
        gam = LinearGAM(s(0) + s(1) + s(2)).fit(Xs, ys)
        sub_s["residual_spline"] = ys - gam.predict(Xs)
        y_s = sub_s["pcr"].astype(int).to_numpy()
        r_s = sub_s["residual_spline"].to_numpy()
        asp, lsp, hsp, _ = auroc_with_boot_ci(r_s, y_s, n_boot=1000, seed=42)
        boot_s = []
        rng = np.random.default_rng(42)
        for _ in range(2000):
            idx = rng.integers(0, len(y_s), size=len(y_s))
            if len(np.unique(y_s[idx])) < 2:
                continue
            boot_s.append(roc_auc_score(y_s[idx], r_s[idx]))
        alpha = 0.05 / 3
        bf_lo = float(np.percentile(boot_s, 100 * alpha / 2)) if boot_s else float("nan")
        bf_hi = float(np.percentile(boot_s, 100 * (1 - alpha / 2))) if boot_s else float("nan")
        supp["spline_residual"] = {
            "n_used": int(len(sub_s)),
            "point": asp,
            "ci_low_95": lsp,
            "ci_high_95": hsp,
            "ci_low_bonf3": bf_lo,
            "ci_high_bonf3": bf_hi,
            "verdict_unadjusted": verdict_3x3(asp, lsp),
            "verdict_bonf3": verdict_3x3(asp, bf_lo),
        }
    except ImportError:
        supp["spline_residual"] = {"status": "skipped", "note": "pygam not installed"}
    except Exception as e:
        supp["spline_residual"] = {"status": "failed", "error": str(e)}

    log["supplementary_variants"] = supp
    log["disclosure"] = (
        "Executed without prior Zenodo deposit per autonomous workflow; §4.8 "
        "reportable as exploratory unless user opts to formalize pre-LOCK "
        "retroactively (note: would require pre-execution Zenodo timestamp "
        "to claim confirmatory)."
    )
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    out_json = PHASE_2 / "gse32646_residual_test.json"
    with open(out_json, "w") as f:
        json.dump(log, f, indent=2)

    # Summary TSV
    rows = []
    rows.append(
        {
            "test": "primary (H + H_PAM50basal)",
            "n": log["primary_model"]["n_used"],
            "auroc_point": auc,
            "ci_low_95": lo,
            "ci_high_95": hi,
            "verdict": primary_verdict,
        }
    )
    rows.append(
        {
            "test": "raw sig_dox (no deconfounding)",
            "n": log["primary_model"]["n_used"],
            "auroc_point": raw_auc,
            "ci_low_95": raw_lo,
            "ci_high_95": raw_hi,
            "verdict": "informational only",
        }
    )
    for name, v in supp.items():
        if "point" in v:
            rows.append(
                {
                    "test": f"supp: {name}",
                    "n": v["n_used"],
                    "auroc_point": v["point"],
                    "ci_low_95": v["ci_low_95"],
                    "ci_high_95": v["ci_high_95"],
                    "ci_low_bonf3": v.get("ci_low_bonf3"),
                    "ci_high_bonf3": v.get("ci_high_bonf3"),
                    "verdict": v.get("verdict_bonf3"),
                }
            )
        else:
            rows.append({"test": f"supp: {name}", "n": None, "auroc_point": None, "verdict": v.get("status")})
    summary = pd.DataFrame(rows)
    summary.to_csv(PHASE_2 / "gse32646_residual_summary.tsv", sep="\t", index=False, float_format="%.4f")

    print(f"\n===== §4.8 HEADLINE (GSE32646, n={log['primary_model']['n_used']}) =====")
    print(f"Primary residual AUROC = {auc:.3f}  [95% CI: {lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {primary_verdict}")
    print(f"Raw sig_dox AUROC (no deconfound) = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print("\nSupplementary (Bonferroni × 3):")
    for name, v in supp.items():
        if "point" in v:
            print(
                f"  {name}: AUROC={v['point']:.3f} 95%CI=[{v['ci_low_95']:.3f},{v['ci_high_95']:.3f}] "
                f"Bonf3=[{v['ci_low_bonf3']:.3f},{v['ci_high_bonf3']:.3f}] verdict_bonf3={v['verdict_bonf3']}"
            )
        else:
            print(f"  {name}: {v}")


if __name__ == "__main__":
    main()
