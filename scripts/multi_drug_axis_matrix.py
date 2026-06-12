"""Phase 1 §4.2 — Multi-drug axis correlation matrix.

For 11 signatures x 3 cohorts:
  Cytotoxics (7): Doxorubicin, Etoposide, Paclitaxel, Carboplatin, Cyclophosphamide,
                  Fluorouracil, TOP2_consensus
  TKIs (3): Trametinib, Gefitinib, Lapatinib
  Positive H control: Fulvestrant
  Negative null: Imatinib

1) Per (drug, cohort): r(sig, P), r(sig, H), r(sig, B) — Pearson + Spearman + boot CI
2) 11x11 inter-signature correlation matrix per cohort (Pearson)
3) Headline test:
   Delta = mean_cytotoxic |r(sig,P)| - mean_(TKI + Imatinib) |r(sig,P)|
4) Pre-registered threshold checks per drug class.

Outputs:
  results/phase_1/multi_drug_axis.json
  results/phase_1/multi_drug_axis_table.tsv
  results/phase_1/correlation_matrix_<cohort>.tsv  (3 files)
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REPO = Path(".")
PHASE0 = REPO / "." / "results" / "phase_0"
OUT_DIR = REPO / "." / "results" / "phase_1"
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]

CYTOTOXICS = ["DOXORUBICIN", "ETOPOSIDE", "PACLITAXEL", "CARBOPLATIN",
              "CYCLOPHOSPHAMIDE", "FLUOROURACIL", "TOP2_consensus"]
TKIS = ["TRAMETINIB", "GEFITINIB", "LAPATINIB"]
H_CTRL = ["FULVESTRANT"]
NULLS = ["IMATINIB"]
DRUGS = CYTOTOXICS + TKIS + H_CTRL + NULLS  # 7 + 3 + 1 + 1 = 12 entries (TOP2_consensus is 7th cytotoxic, total 11 if no overlap with Dox; we report all 12)


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_score_file(drug: str, cohort: str) -> Path:
    """Search for drug score file in phase_0 then external_data."""
    candidates = [
        PHASE0 / f"score_{drug}_{cohort}.tsv",
        REPO / "external_data" / cohort / f"score_{drug}.tsv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def boot_r(x, y, n_boot=1000, seed=42, method="pearson"):
    """Bootstrap CI for Pearson or Spearman r."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 5:
        return float("nan"), float("nan"), float("nan")
    if method == "pearson":
        r_obs = pearsonr(x, y).statistic
    else:
        r_obs = spearmanr(x, y).statistic
    boot = []
    n = len(x)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            if method == "pearson":
                boot.append(pearsonr(x[idx], y[idx]).statistic)
            else:
                boot.append(spearmanr(x[idx], y[idx]).statistic)
        except Exception:
            continue
    if len(boot) < 10:
        return float(r_obs), float("nan"), float("nan")
    return float(r_obs), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {
        "script": str(me),
        "script_sha256": sha256_of(me),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "drugs_cytotoxic": CYTOTOXICS,
        "drugs_tki": TKIS,
        "drugs_h_control": H_CTRL,
        "drugs_null": NULLS,
        "cohorts": {},
    }

    table_rows = []
    delta_per_cohort = {}

    for cohort in COHORTS:
        axes = pd.read_csv(PHASE0 / f"axes_{cohort}.tsv", sep="\t")
        # Build per-cohort score panel
        score_panel = {}
        missing = []
        for d in DRUGS:
            p = find_score_file(d, cohort)
            if p is None:
                missing.append(d)
                continue
            df = pd.read_csv(p, sep="\t")
            score_panel[d] = df.set_index("patient_id")["signature_score"]

        if missing:
            log.setdefault("warnings", []).append(
                f"{cohort}: missing scores for {missing}")
        # Restrict to common patients across all available drugs + axes
        common_pts = set(axes["sample_id"].astype(str))
        for d, s in score_panel.items():
            common_pts &= set(s.index.astype(str))
        common_pts = sorted(common_pts)
        axes_idx = axes.set_index("sample_id").loc[common_pts]
        P = axes_idx["P"].astype(float).values
        H = axes_idx["H"].astype(float).values
        B = axes_idx["B"].astype(float).values

        # Per-drug correlations
        cohort_block = {"n_patients": len(common_pts), "per_drug": {}}
        for d in DRUGS:
            if d not in score_panel:
                cohort_block["per_drug"][d] = {"missing": True}
                continue
            sig = score_panel[d].loc[common_pts].astype(float).values
            r_p, lo_p, hi_p = boot_r(sig, P, 1000, 42, "pearson")
            r_h, lo_h, hi_h = boot_r(sig, H, 1000, 42, "pearson")
            r_b, lo_b, hi_b = boot_r(sig, B, 1000, 42, "pearson")
            rs_p, _, _ = boot_r(sig, P, 200, 42, "spearman")
            rs_h, _, _ = boot_r(sig, H, 200, 42, "spearman")
            rs_b, _, _ = boot_r(sig, B, 200, 42, "spearman")
            cohort_block["per_drug"][d] = {
                "pearson_r_P": r_p, "pearson_r_P_ci": [lo_p, hi_p],
                "pearson_r_H": r_h, "pearson_r_H_ci": [lo_h, hi_h],
                "pearson_r_B": r_b, "pearson_r_B_ci": [lo_b, hi_b],
                "spearman_r_P": rs_p, "spearman_r_H": rs_h, "spearman_r_B": rs_b,
            }
            klass = ("cytotoxic" if d in CYTOTOXICS else
                     "TKI" if d in TKIS else
                     "H_control" if d in H_CTRL else "null")
            table_rows.append({
                "cohort": cohort, "drug": d, "class": klass,
                "r_P": r_p, "r_P_ci_lo": lo_p, "r_P_ci_hi": hi_p,
                "r_H": r_h, "r_H_ci_lo": lo_h, "r_H_ci_hi": hi_h,
                "r_B": r_b, "r_B_ci_lo": lo_b, "r_B_ci_hi": hi_b,
                "rho_P": rs_p, "rho_H": rs_h, "rho_B": rs_b,
                "abs_r_P": abs(r_p),
            })

        # Inter-signature matrix
        present = [d for d in DRUGS if d in score_panel]
        mat_df = pd.DataFrame({d: score_panel[d].loc[common_pts].astype(float).values for d in present})
        corr = mat_df.corr(method="pearson")
        corr.to_csv(OUT_DIR / f"correlation_matrix_{cohort}.tsv", sep="\t", float_format="%.6g")

        # Delta cytotoxic - TKI+null
        rP_cyt = [abs(cohort_block["per_drug"][d]["pearson_r_P"])
                  for d in CYTOTOXICS if d in score_panel
                  and np.isfinite(cohort_block["per_drug"][d]["pearson_r_P"])]
        rP_oth = [abs(cohort_block["per_drug"][d]["pearson_r_P"])
                  for d in TKIS + NULLS if d in score_panel
                  and np.isfinite(cohort_block["per_drug"][d]["pearson_r_P"])]
        delta = float(np.mean(rP_cyt) - np.mean(rP_oth)) if rP_cyt and rP_oth else float("nan")
        cohort_block["headline_delta_cytotoxic_minus_TKInull"] = delta
        cohort_block["mean_abs_rP_cytotoxic"] = float(np.mean(rP_cyt)) if rP_cyt else float("nan")
        cohort_block["mean_abs_rP_TKI_null"] = float(np.mean(rP_oth)) if rP_oth else float("nan")
        delta_per_cohort[cohort] = delta

        # Pre-reg threshold checks
        checks = {}
        # All 6 cytotoxics |r(sig,P)| > 0.6 (exclude TOP2_consensus from this count; 6 single drugs)
        single_cyt = ["DOXORUBICIN", "ETOPOSIDE", "PACLITAXEL", "CARBOPLATIN",
                      "CYCLOPHOSPHAMIDE", "FLUOROURACIL"]
        checks["cytotoxics_all_rP_gt_06"] = {
            d: (abs(cohort_block["per_drug"][d]["pearson_r_P"]) > 0.6
                if d in score_panel and np.isfinite(cohort_block["per_drug"][d]["pearson_r_P"]) else None)
            for d in single_cyt
        }
        if "TRAMETINIB" in score_panel:
            checks["trametinib_abs_rP_lt_05"] = bool(
                abs(cohort_block["per_drug"]["TRAMETINIB"]["pearson_r_P"]) < 0.5)
        if "LAPATINIB" in score_panel:
            checks["lapatinib_rB_gt_03"] = bool(
                cohort_block["per_drug"]["LAPATINIB"]["pearson_r_B"] > 0.3)
        if "FULVESTRANT" in score_panel:
            checks["fulvestrant_rH_gt_03"] = bool(
                cohort_block["per_drug"]["FULVESTRANT"]["pearson_r_H"] > 0.3)
        if "IMATINIB" in score_panel:
            v = cohort_block["per_drug"]["IMATINIB"]
            checks["imatinib_all_abs_r_lt_03"] = bool(
                abs(v["pearson_r_P"]) < 0.3 and abs(v["pearson_r_H"]) < 0.3 and abs(v["pearson_r_B"]) < 0.3)
        cohort_block["preregistered_checks"] = checks
        log["cohorts"][cohort] = cohort_block

    # Cross-cohort weighted mean (weight by n_patients)
    weights = {c: log["cohorts"][c]["n_patients"] for c in COHORTS}
    total_w = sum(weights.values())
    cross_delta = sum(delta_per_cohort[c] * weights[c] for c in COHORTS
                      if np.isfinite(delta_per_cohort[c])) / total_w if total_w > 0 else float("nan")
    log["cross_cohort_weighted_delta"] = cross_delta
    log["delta_per_cohort"] = delta_per_cohort

    # BH-corrected paired comparison (cytotoxic vs TKI+null abs_rP arrays per cohort)
    # Use a paired Wilcoxon-like quick check (per-cohort)
    from scipy.stats import mannwhitneyu
    pvals = []
    for cohort in COHORTS:
        cb = log["cohorts"][cohort]
        rP_cyt = [abs(cb["per_drug"][d]["pearson_r_P"]) for d in CYTOTOXICS
                  if d in cb["per_drug"] and "pearson_r_P" in cb["per_drug"][d]
                  and np.isfinite(cb["per_drug"][d]["pearson_r_P"])]
        rP_oth = [abs(cb["per_drug"][d]["pearson_r_P"]) for d in TKIS + NULLS
                  if d in cb["per_drug"] and "pearson_r_P" in cb["per_drug"][d]
                  and np.isfinite(cb["per_drug"][d]["pearson_r_P"])]
        if len(rP_cyt) >= 2 and len(rP_oth) >= 2:
            res = mannwhitneyu(rP_cyt, rP_oth, alternative="greater")
            pvals.append((cohort, float(res.pvalue)))
    # BH
    from statsmodels.stats.multitest import multipletests
    if pvals:
        ps = [p for _, p in pvals]
        rej, p_adj, _, _ = multipletests(ps, method="fdr_bh")
        log["paired_class_comparison_BH"] = [
            {"cohort": c, "p_raw": p, "p_BH": float(pa), "reject_BH": bool(r)}
            for (c, p), pa, r in zip(pvals, p_adj, rej)
        ]

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(OUT_DIR / "multi_drug_axis.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    pd.DataFrame(table_rows).to_csv(
        OUT_DIR / "multi_drug_axis_table.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"A2 done. Delta per cohort: {delta_per_cohort}")
    print(f"Cross-cohort weighted delta = {cross_delta:.3f}")


if __name__ == "__main__":
    main()
