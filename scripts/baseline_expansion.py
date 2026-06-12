"""Phase 1 §4.3 — Baseline expansion.

For each cohort, compute AUROC vs pCR for:
  1. Random uniform (control) — single random draw, seed 42
  2. MKI67 single-gene
  3. HALLMARK_E2F_TARGETS singscore
  4. HALLMARK_G2M_CHECKPOINT singscore
  5. Mean-expression z-score across DOX signature genes (no sign) — NEW
  6. 5-set composite proliferation (P axis from Phase 0)
  7. Random 60 OncoKB genes with random signs — 100 iterations, NULL distribution
  8. Random Hallmark set baseline — 50 iterations
  9. PRIMARY COMPARATOR (HEADLINE): MKI67 + ESR1 two-gene sum z-score — NEW

Comparator: Dox framework signature (existing).
Significance: paired DeLong-style bootstrap p-values vs framework signature AUROC.

Outputs:
  results/phase_1/baseline_expansion.json
  results/phase_1/baseline_expansion.tsv
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path(".")
PHASE0 = REPO / "." / "results" / "phase_0"
OUT_DIR = REPO / "." / "results" / "phase_1"
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]
SIG_DOX = REPO / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all" / "DOXORUBICIN.tsv"
HALLMARK = REPO / "h.all.v2025.1.Hs.symbols.gmt"
ONCOKB = REPO / "external_data" / "oncokb_all_universe.tsv"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_expr(cohort: str) -> pd.DataFrame:
    df = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_expression.tsv",
                     sep="\t", index_col=0)
    df.columns = [c.upper() for c in df.columns]
    return df


def load_pcr(cohort: str) -> pd.DataFrame:
    return pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_pcr.tsv", sep="\t")


def load_dox_score(cohort: str) -> pd.DataFrame:
    return pd.read_csv(REPO / "external_data" / cohort / "score_DOXORUBICIN.tsv", sep="\t")


def load_hallmark_set(name: str):
    with open(HALLMARK) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == name:
                return set(g.upper() for g in parts[2:])
    return set()


def load_all_hallmark():
    sets = {}
    with open(HALLMARK) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = set(g.upper() for g in parts[2:])
    return sets


def singscore_set(expr: pd.DataFrame, gene_set: set) -> pd.Series:
    """Singscore-style for a single gene set (no down direction)."""
    ranks = expr.rank(axis=1, method="average", pct=True)
    present = sorted(gene_set & set(expr.columns))
    if len(present) < 3:
        return pd.Series(dtype=float, index=expr.index)
    n = expr.shape[1]
    mn = (1 + len(present)) / (2 * n); mx = 1 - mn
    mean_rank = ranks[present].mean(axis=1)
    return 2 * (mean_rank - mn) / (mx - mn) - 1


def mean_z(expr: pd.DataFrame, genes: list) -> pd.Series:
    present = [g for g in genes if g in expr.columns]
    if len(present) < 3:
        return pd.Series(dtype=float, index=expr.index)
    sub = expr[present]
    zs = (sub - sub.mean(axis=0)) / (sub.std(axis=0, ddof=0) + 1e-12)
    return zs.mean(axis=1)


def boot_auroc(y, s, n_boot=1000, seed=42):
    rng = np.random.default_rng(seed)
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    obs = float(roc_auc_score(y, s))
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


def paired_boot_p(y, s_a, s_b, n_boot=1000, seed=42):
    """Paired bootstrap: p-value for AUROC(A) - AUROC(B) > 0 (two-sided)."""
    rng = np.random.default_rng(seed)
    if len(np.unique(y)) < 2:
        return float("nan")
    obs = float(roc_auc_score(y, s_a) - roc_auc_score(y, s_b))
    diffs = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        diffs.append(roc_auc_score(y[idx], s_a[idx]) - roc_auc_score(y[idx], s_b[idx]))
    if len(diffs) < 10:
        return float("nan")
    diffs = np.array(diffs)
    # Two-sided p-value (centered at 0)
    p_two = 2.0 * min(float((diffs <= 0).sum() + 1) / (len(diffs) + 1),
                     float((diffs >= 0).sum() + 1) / (len(diffs) + 1))
    return min(p_two, 1.0), obs


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {"script": str(me), "script_sha256": sha256_of(me),
           "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
           "cohorts": {}}

    sig_df = pd.read_csv(SIG_DOX, sep="\t")
    sig_genes = [str(g).upper() for g in sig_df["gene_symbol"]]
    e2f_set = load_hallmark_set("HALLMARK_E2F_TARGETS")
    g2m_set = load_hallmark_set("HALLMARK_G2M_CHECKPOINT")
    all_hm = load_all_hallmark()
    oncokb = pd.read_csv(ONCOKB, sep="\t").iloc[:, 0].astype(str).str.upper().tolist()

    rows = []
    for cohort in COHORTS:
        expr = load_expr(cohort)
        pcr = load_pcr(cohort)
        dox = load_dox_score(cohort)
        axes = pd.read_csv(PHASE0 / f"axes_{cohort}.tsv", sep="\t")

        # Common patients across expr + dox + pcr + axes
        common = (set(expr.index.astype(str)) &
                  set(pcr["patient_id"].astype(str)) &
                  set(dox["patient_id"].astype(str)) &
                  set(axes["sample_id"].astype(str)))
        common = sorted(common)
        expr_c = expr.loc[common]
        pcr_c = pcr.set_index("patient_id").loc[common]["pcr"].astype(int).values
        framework = dox.set_index("patient_id").loc[common]["signature_score"].astype(float).values
        P_axis = axes.set_index("sample_id").loc[common]["P"].astype(float).values

        if len(np.unique(pcr_c)) < 2:
            continue
        n = len(common)

        # 1. Random baseline (single seeded draw)
        rng_rand = np.random.default_rng(42)
        rand_score = rng_rand.uniform(0, 1, size=n)
        # 2. MKI67
        mki67 = expr_c["MKI67"].astype(float).values if "MKI67" in expr_c.columns else None
        # 3. E2F
        e2f = singscore_set(expr_c, e2f_set).loc[common].astype(float).values
        # 4. G2M
        g2m = singscore_set(expr_c, g2m_set).loc[common].astype(float).values
        # 5. Mean-z across DOX signature genes
        mz = mean_z(expr_c, sig_genes).loc[common].astype(float).values
        # 6. P axis
        p_axis_score = P_axis
        # 9. MKI67 + ESR1 two-gene sum z-score (HEADLINE)
        if "MKI67" in expr_c.columns and "ESR1" in expr_c.columns:
            zmki = (expr_c["MKI67"] - expr_c["MKI67"].mean()) / (expr_c["MKI67"].std(ddof=0) + 1e-12)
            zesr = (expr_c["ESR1"] - expr_c["ESR1"].mean()) / (expr_c["ESR1"].std(ddof=0) + 1e-12)
            mki_esr_two_gene = (zmki - zesr).loc[common].astype(float).values  # high prolif, low ER
        else:
            mki_esr_two_gene = None

        baselines = {
            "PRIMARY_DOXORUBICIN_framework": framework,
            "random_uniform_seed42": rand_score,
            "MKI67_single_gene": mki67,
            "HALLMARK_E2F_TARGETS": e2f,
            "HALLMARK_G2M_CHECKPOINT": g2m,
            "DOX_genes_mean_z": mz,
            "P_axis_5set_composite": p_axis_score,
            "MKI67_minus_ESR1_two_gene_z": mki_esr_two_gene,
        }

        per_cohort = {}
        for label, s in baselines.items():
            if s is None:
                per_cohort[label] = {"missing": True}
                continue
            mask = np.isfinite(s) & np.isfinite(pcr_c)
            ys, ss = pcr_c[mask], s[mask]
            auc, lo, hi = boot_auroc(ys, ss, 1000, 42)
            p_paired = None
            if label != "PRIMARY_DOXORUBICIN_framework":
                framework_aligned = framework[mask]
                p_paired, delta = paired_boot_p(ys, ss, framework_aligned, 1000, 42)
            else:
                p_paired, delta = None, 0.0
            per_cohort[label] = {
                "AUROC": auc, "AUROC_ci_95": [lo, hi], "n": int(mask.sum()),
                "delta_vs_framework": float(delta) if delta is not None else None,
                "paired_p_vs_framework": float(p_paired) if p_paired is not None else None,
                "headline_tag": "HEADLINE" if label == "MKI67_minus_ESR1_two_gene_z" else None,
            }
            rows.append({"cohort": cohort, "baseline": label, **per_cohort[label]})

        # 7. Random OncoKB 60-gene baseline x 100
        oncokb_present = [g for g in oncokb if g in expr_c.columns]
        rng_o = np.random.default_rng(123)
        oncokb_aucs = []
        for _ in range(100):
            sample = rng_o.choice(oncokb_present, size=60, replace=False)
            signs = rng_o.choice([-1, 1], size=60)
            mat = expr_c[list(sample)].astype(float)
            zs = (mat - mat.mean(axis=0)) / (mat.std(axis=0, ddof=0) + 1e-12)
            sscore = (zs * signs).mean(axis=1).loc[common].values
            if np.isfinite(sscore).all() and len(np.unique(pcr_c)) >= 2:
                oncokb_aucs.append(float(roc_auc_score(pcr_c, sscore)))
        per_cohort["random_oncokb_60_x100"] = {
            "mean_AUROC": float(np.mean(oncokb_aucs)) if oncokb_aucs else float("nan"),
            "median_AUROC": float(np.median(oncokb_aucs)) if oncokb_aucs else float("nan"),
            "auc_distribution_percentiles": {
                "p2.5": float(np.percentile(oncokb_aucs, 2.5)) if oncokb_aucs else float("nan"),
                "p97.5": float(np.percentile(oncokb_aucs, 97.5)) if oncokb_aucs else float("nan"),
            },
            "p_value_vs_framework": float(
                (np.sum(np.array(oncokb_aucs) >= roc_auc_score(pcr_c, framework)) + 1)
                / (len(oncokb_aucs) + 1)) if oncokb_aucs else float("nan"),
            "n_iterations": len(oncokb_aucs),
        }
        rows.append({"cohort": cohort, "baseline": "random_oncokb_60_x100",
                     "AUROC": per_cohort["random_oncokb_60_x100"]["mean_AUROC"],
                     "AUROC_ci_95": [
                         per_cohort["random_oncokb_60_x100"]["auc_distribution_percentiles"]["p2.5"],
                         per_cohort["random_oncokb_60_x100"]["auc_distribution_percentiles"]["p97.5"]],
                     "n": int(len(common)),
                     "paired_p_vs_framework": per_cohort["random_oncokb_60_x100"]["p_value_vs_framework"],
                     "headline_tag": None})

        # 8. Random Hallmark set x 50
        hm_names = list(all_hm.keys())
        rng_h = np.random.default_rng(456)
        hm_aucs = []
        chosen = rng_h.choice(hm_names, size=50, replace=False)
        for nm in chosen:
            sscore = singscore_set(expr_c, all_hm[nm]).loc[common].values
            if np.isfinite(sscore).all() and len(np.unique(pcr_c)) >= 2:
                hm_aucs.append(float(roc_auc_score(pcr_c, sscore)))
        per_cohort["random_hallmark_x50"] = {
            "mean_AUROC": float(np.mean(hm_aucs)) if hm_aucs else float("nan"),
            "median_AUROC": float(np.median(hm_aucs)) if hm_aucs else float("nan"),
            "auc_distribution_percentiles": {
                "p2.5": float(np.percentile(hm_aucs, 2.5)) if hm_aucs else float("nan"),
                "p97.5": float(np.percentile(hm_aucs, 97.5)) if hm_aucs else float("nan"),
            },
            "p_value_vs_framework": float(
                (np.sum(np.array(hm_aucs) >= roc_auc_score(pcr_c, framework)) + 1)
                / (len(hm_aucs) + 1)) if hm_aucs else float("nan"),
            "n_iterations": len(hm_aucs),
        }
        rows.append({"cohort": cohort, "baseline": "random_hallmark_x50",
                     "AUROC": per_cohort["random_hallmark_x50"]["mean_AUROC"],
                     "AUROC_ci_95": [
                         per_cohort["random_hallmark_x50"]["auc_distribution_percentiles"]["p2.5"],
                         per_cohort["random_hallmark_x50"]["auc_distribution_percentiles"]["p97.5"]],
                     "n": int(len(common)),
                     "paired_p_vs_framework": per_cohort["random_hallmark_x50"]["p_value_vs_framework"],
                     "headline_tag": None})

        log["cohorts"][cohort] = per_cohort

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(OUT_DIR / "baseline_expansion.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    pd.DataFrame(rows).to_csv(OUT_DIR / "baseline_expansion.tsv", sep="\t",
                              index=False, float_format="%.6g")
    print("A3 done. HEADLINE MKI67-ESR1 vs framework:")
    for c in COHORTS:
        b = log["cohorts"].get(c, {})
        fw = b.get("PRIMARY_DOXORUBICIN_framework", {}).get("AUROC", float("nan"))
        h = b.get("MKI67_minus_ESR1_two_gene_z", {})
        ha = h.get("AUROC", float("nan"))
        print(f"  {c}: framework={fw:.3f}  MKI67-ESR1={ha:.3f}  delta={ha-fw:+.3f}")


if __name__ == "__main__":
    main()
