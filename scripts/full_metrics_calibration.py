"""Phase 1 §4.4 — Full metrics + calibration.

Per (signature × cohort × pCR):
  - PR-AUC + average precision
  - F1 at Youden-optimal threshold
  - Brier score + reliability data (10 bins)
  - Permutation null:
      * MARGINAL for GSE16446 (uniformly ER-/TN)
      * WITHIN-PAM50 STRATUM for GSE25066, GSE22226
      Report BOTH within-stratum and marginal empirical p (where applicable).
  - CV protocol clarification: group-by-patient (documented; primary already used this)

Signatures evaluated:
  - DOX framework
  - MKI67+ESR1 two-gene (PRIMARY COMPARATOR, NEW v0.3 headline)
  - random uniform
  - MKI67 single-gene
  - HALLMARK_E2F_TARGETS
  - HALLMARK_G2M_CHECKPOINT

Outputs:
  results/phase_1/full_metrics.json
  results/phase_1/full_metrics_table.tsv
  results/phase_1/reliability_data.json
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, average_precision_score, precision_recall_curve,
    roc_curve, brier_score_loss, f1_score, precision_score, recall_score,
)

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


def youden_threshold(y, s):
    fpr, tpr, thr = roc_curve(y, s)
    j = tpr - fpr
    idx = int(np.argmax(j))
    return thr[idx]


def to_prob(s):
    """Min-max scale s to [0,1] for Brier."""
    s = np.asarray(s, dtype=float)
    lo, hi = np.nanmin(s), np.nanmax(s)
    if hi - lo < 1e-12:
        return np.full_like(s, 0.5)
    return (s - lo) / (hi - lo)


def reliability_bins(y, p, n_bins=10):
    edges = np.linspace(0, 1, n_bins + 1)
    rows = []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        mask = (p >= lo) & (p < hi if i < n_bins - 1 else p <= hi)
        if mask.sum() == 0:
            rows.append({"bin": i, "lo": float(lo), "hi": float(hi),
                         "n": 0, "obs_frac": None, "mean_pred": None})
        else:
            rows.append({"bin": i, "lo": float(lo), "hi": float(hi),
                         "n": int(mask.sum()),
                         "obs_frac": float(y[mask].mean()),
                         "mean_pred": float(p[mask].mean())})
    return rows


def perm_p_marginal(y, s, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    obs = roc_auc_score(y, s)
    aucs = []
    yy = y.copy()
    for _ in range(n):
        rng.shuffle(yy)
        if len(np.unique(yy)) >= 2:
            aucs.append(roc_auc_score(yy, s))
    aucs = np.array(aucs)
    return float((np.sum(aucs >= obs) + 1) / (len(aucs) + 1))


def perm_p_within_stratum(y, s, strata, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    obs = roc_auc_score(y, s)
    aucs = []
    for _ in range(n):
        yp = y.copy()
        for ss in np.unique(strata):
            idx = np.where(strata == ss)[0]
            yp[idx] = rng.permutation(yp[idx])
        if len(np.unique(yp)) >= 2:
            aucs.append(roc_auc_score(yp, s))
    aucs = np.array(aucs)
    return float((np.sum(aucs >= obs) + 1) / (len(aucs) + 1))


def load_hallmark_set(name: str):
    with open(REPO / "h.all.v2025.1.Hs.symbols.gmt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == name:
                return set(g.upper() for g in parts[2:])
    return set()


def singscore_set(expr: pd.DataFrame, gene_set: set) -> pd.Series:
    ranks = expr.rank(axis=1, method="average", pct=True)
    present = sorted(gene_set & set(expr.columns))
    if len(present) < 3:
        return pd.Series(dtype=float, index=expr.index)
    n = expr.shape[1]
    mn = (1 + len(present)) / (2 * n); mx = 1 - mn
    mean_rank = ranks[present].mean(axis=1)
    return 2 * (mean_rank - mn) / (mx - mn) - 1


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {"script": str(me), "script_sha256": sha256_of(me),
           "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
           "cv_protocol": "group-by-patient (each patient appears in only one fold); primary pipeline already uses this",
           "cohorts": {}}
    reliability = {}
    rows = []

    e2f_set = load_hallmark_set("HALLMARK_E2F_TARGETS")
    g2m_set = load_hallmark_set("HALLMARK_G2M_CHECKPOINT")

    for cohort in COHORTS:
        expr = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_expression.tsv",
                           sep="\t", index_col=0)
        expr.columns = [c.upper() for c in expr.columns]
        pcr = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_pcr.tsv", sep="\t")
        dox = pd.read_csv(REPO / "external_data" / cohort / "score_DOXORUBICIN.tsv", sep="\t")
        pam50 = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_pam50.tsv", sep="\t")

        common = sorted(set(expr.index.astype(str)) & set(pcr["patient_id"].astype(str))
                        & set(dox["patient_id"].astype(str)))
        expr_c = expr.loc[common]
        y = pcr.set_index("patient_id").loc[common]["pcr"].astype(int).values

        # Stratifier: PAM50 (NA for GSE16446 fallback to marginal-only)
        pam50_map = pam50.set_index("patient_id")["pam50_subtype"].astype(str).to_dict()
        strata = np.array([pam50_map.get(p, "NA") for p in common])
        single_stratum = len(np.unique(strata)) < 2 or cohort == "GSE16446"

        framework = dox.set_index("patient_id").loc[common]["signature_score"].astype(float).values
        rng = np.random.default_rng(42)
        random_score = rng.uniform(0, 1, size=len(common))
        mki67 = expr_c["MKI67"].astype(float).values if "MKI67" in expr_c.columns else None
        e2f = singscore_set(expr_c, e2f_set).loc[common].astype(float).values
        g2m = singscore_set(expr_c, g2m_set).loc[common].astype(float).values
        if "MKI67" in expr_c.columns and "ESR1" in expr_c.columns:
            zmki = (expr_c["MKI67"] - expr_c["MKI67"].mean()) / (expr_c["MKI67"].std(ddof=0) + 1e-12)
            zesr = (expr_c["ESR1"] - expr_c["ESR1"].mean()) / (expr_c["ESR1"].std(ddof=0) + 1e-12)
            mki_esr = (zmki - zesr).loc[common].astype(float).values
        else:
            mki_esr = None

        sigs = {
            "PRIMARY_DOXORUBICIN_framework": framework,
            "MKI67_minus_ESR1_two_gene_z_HEADLINE": mki_esr,
            "random_uniform_seed42": random_score,
            "MKI67_single_gene": mki67,
            "HALLMARK_E2F_TARGETS": e2f,
            "HALLMARK_G2M_CHECKPOINT": g2m,
        }

        cohort_block = {"n_total": len(common),
                        "n_pcr_positive": int(y.sum()),
                        "single_stratum_mode": bool(single_stratum),
                        "per_signature": {}}
        reliability[cohort] = {}

        for label, s in sigs.items():
            if s is None:
                cohort_block["per_signature"][label] = {"missing": True}
                continue
            mask = np.isfinite(s)
            ys, ss = y[mask], s[mask]
            strata_ss = strata[mask]
            if len(np.unique(ys)) < 2:
                continue
            auroc = float(roc_auc_score(ys, ss))
            pr_auc = float(average_precision_score(ys, ss))
            thr = youden_threshold(ys, ss)
            y_pred = (ss >= thr).astype(int)
            f1 = float(f1_score(ys, y_pred))
            prec = float(precision_score(ys, y_pred, zero_division=0))
            rec = float(recall_score(ys, y_pred, zero_division=0))
            p_prob = to_prob(ss)
            brier = float(brier_score_loss(ys, p_prob))
            rel = reliability_bins(ys, p_prob, 10)
            reliability[cohort][label] = rel

            # Permutation null
            perm_marginal = perm_p_marginal(ys, ss, 1000, 42)
            if single_stratum:
                perm_within = None
            else:
                perm_within = perm_p_within_stratum(ys, ss, strata_ss, 1000, 42)

            block = {
                "AUROC": auroc, "PR_AUC": pr_auc, "average_precision": pr_auc,
                "youden_threshold": float(thr), "F1_at_youden": f1,
                "precision_at_youden": prec, "recall_at_youden": rec,
                "brier_score_on_minmax_prob": brier,
                "permutation_p_marginal": perm_marginal,
                "permutation_p_within_pam50": perm_within,
                "n_used": int(mask.sum()),
            }
            cohort_block["per_signature"][label] = block
            rows.append({"cohort": cohort, "signature": label,
                         "AUROC": auroc, "PR_AUC": pr_auc, "F1_at_youden": f1,
                         "brier": brier, "perm_p_marginal": perm_marginal,
                         "perm_p_within_pam50": perm_within,
                         "n": int(mask.sum())})
        log["cohorts"][cohort] = cohort_block

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(OUT_DIR / "full_metrics.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    pd.DataFrame(rows).to_csv(OUT_DIR / "full_metrics_table.tsv", sep="\t",
                              index=False, float_format="%.6g")
    with open(OUT_DIR / "reliability_data.json", "w") as f:
        json.dump(reliability, f, indent=2, default=str)
    print("A4 done")
    for c in COHORTS:
        cb = log["cohorts"][c]
        fw = cb["per_signature"].get("PRIMARY_DOXORUBICIN_framework", {})
        print(f"  {c}: fw AUROC={fw.get('AUROC', 'NA'):.3f} PR-AUC={fw.get('PR_AUC', 'NA'):.3f} "
              f"Brier={fw.get('brier_score_on_minmax_prob', 'NA'):.3f} "
              f"perm_p_marg={fw.get('permutation_p_marginal', 'NA')} "
              f"perm_p_within={fw.get('permutation_p_within_pam50', 'NA')}")


if __name__ == "__main__":
    main()
