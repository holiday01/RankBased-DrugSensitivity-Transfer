"""EXPLORATORY: Incremental value of drug signature OVER proliferation baseline.

NOT in pre-reg primary family. Reported as supplementary mechanism evidence.

For each (cohort, drug signature) combination, fit:
    glm(pcr ~ drug_signature_z + E2F_TARGETS_z + PAM50_basal_z)

Report:
  - drug coefficient, Wald p, OR per SD
  - joint model AUROC
  - baseline-only model AUROC (pcr ~ E2F + Basal)
  - Delta AUROC (joint - baseline-only)
  - F-test/likelihood-ratio comparing joint vs baseline-only

Interpretation:
  - If drug coefficient p < 0.05 in ≥ 2 of 3 cohorts → framework provides
    incremental information independent of proliferation
  - If joint AUROC > baseline AUROC by ≥ +0.03 → meaningful clinical
    contribution beyond what proliferation alone provides
"""

from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
from scipy.stats import norm, chi2

REVISE_ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv")
EXT = REVISE_ROOT / "external_data"
HALLMARK = Path("/home/holiday01/2026_ISMB_code/enrichment_test/h.all.v2025.1.Hs.symbols.gmt")
FROZEN = REVISE_ROOT / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all"
PANDOX = REVISE_ROOT / "frozen_signatures" / "quantile_30_oncokb_all"  # pan-cancer


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


def zscore(x):
    sd = np.std(x, ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def fit_glm(y, X_cols, design):
    """design = dict of name→array. Add intercept automatically. Return params dict."""
    X = np.column_stack([np.ones(len(y))] + [design[c] for c in X_cols])
    mod = sm.GLM(y, X, family=sm.families.Binomial())
    res = mod.fit(disp=0)
    pred = res.predict(X)
    out = {"loglik": float(res.llf), "df": int(res.df_resid), "auroc": float(roc_auc_score(y, pred))}
    for i, name in enumerate(["Intercept"] + X_cols):
        out[f"{name}_coef"] = float(res.params[i])
        out[f"{name}_se"] = float(res.bse[i])
        out[f"{name}_p"] = float(2 * (1 - norm.cdf(abs(res.params[i] / res.bse[i])))) if res.bse[i] > 0 else float("nan")
        out[f"{name}_or"] = float(np.exp(res.params[i]))
    return out


def lr_test(loglik_full, loglik_reduced, df_diff=1):
    stat = 2 * (loglik_full - loglik_reduced)
    p = 1 - chi2.cdf(stat, df_diff)
    return float(stat), float(p)


COHORTS = {
    "GSE16446": dict(expr=EXT/"GSE16446/GSE16446_expression.tsv",
                     pcr=EXT/"GSE16446/GSE16446_pcr.tsv",
                     pam50=EXT/"GSE16446/GSE16446_pam50.tsv"),
    "GSE25066": dict(expr=EXT/"GSE25066/GSE25066_expression.tsv",
                     pcr=EXT/"GSE25066/GSE25066_pcr.tsv",
                     pam50=EXT/"GSE25066/GSE25066_pam50.tsv"),
    "GSE22226": dict(expr=EXT/"GSE22226/GSE22226_expression.tsv",
                     pcr=EXT/"GSE22226/GSE22226_pcr.tsv",
                     pam50=EXT/"GSE22226/GSE22226_pam50.tsv"),
}

DRUG_SIGS = {
    "DOXORUBICIN":         FROZEN / "DOXORUBICIN.tsv",
    "PACLITAXEL":          FROZEN / "PACLITAXEL.tsv",
    "TOP2_consensus":      FROZEN / "TOP2_POISON_CONSENSUS.tsv",
    "ETOPOSIDE":           FROZEN / "ETOPOSIDE.tsv",
}


def main():
    e2f_genes = read_gmt_set(HALLMARK, "HALLMARK_E2F_TARGETS")
    g2m_genes = read_gmt_set(HALLMARK, "HALLMARK_G2M_CHECKPOINT")
    print(f"Loaded E2F_TARGETS ({len(e2f_genes)} genes), G2M_CHECKPOINT ({len(g2m_genes)} genes)")

    results = []

    for cohort_name, paths in COHORTS.items():
        print(f"\n=== {cohort_name} ===")
        expr = pd.read_csv(paths["expr"], sep="\t", index_col=0)
        expr.columns = [c.upper() for c in expr.columns]
        pcr_df = pd.read_csv(paths["pcr"], sep="\t").dropna(subset=["pcr"])
        pam50 = pd.read_csv(paths["pam50"], sep="\t")
        df = pcr_df[["patient_id", "pcr"]].merge(pam50[["patient_id", "basal_score"]], on="patient_id", how="inner")
        df = df.dropna(subset=["pcr", "basal_score"])
        common = sorted(set(df["patient_id"]) & set(expr.index.astype(str)))
        df = df.set_index("patient_id").loc[common]
        expr_aligned = expr.loc[common]

        # E2F and G2M singscores per patient
        e2f_score = singscore_unidirectional(expr_aligned, e2f_genes).loc[common].to_numpy()
        g2m_score = singscore_unidirectional(expr_aligned, g2m_genes).loc[common].to_numpy()

        y = df["pcr"].astype(int).to_numpy()
        basal_z = zscore(df["basal_score"].to_numpy())
        e2f_z = zscore(e2f_score)
        g2m_z = zscore(g2m_score)
        print(f"  n={len(y)}, pcr_rate={y.mean():.3f}")

        # Baseline-only models (no drug signature)
        baseline_e2f = fit_glm(y, ["E2F", "BASAL"], {"E2F": e2f_z, "BASAL": basal_z})
        baseline_g2m = fit_glm(y, ["G2M", "BASAL"], {"G2M": g2m_z, "BASAL": basal_z})
        # Pick whichever proliferation baseline gives higher AUROC (consistent with score_baselines.py "best baseline" rule)
        if baseline_e2f["auroc"] >= baseline_g2m["auroc"]:
            baseline_label, baseline_results, prolif_z, prolif_name = "E2F_TARGETS", baseline_e2f, e2f_z, "E2F"
        else:
            baseline_label, baseline_results, prolif_z, prolif_name = "G2M_CHECKPOINT", baseline_g2m, g2m_z, "G2M"
        print(f"  Best proliferation baseline = {baseline_label}, AUROC = {baseline_results['auroc']:.3f}")

        for drug_name, sig_path in DRUG_SIGS.items():
            if not sig_path.exists():
                continue
            sens, res = load_signature(sig_path)
            try:
                drug_score = singscore_bidirectional(expr_aligned, sens, res).loc[common].to_numpy()
            except Exception as e:
                print(f"    {drug_name}: ERROR {e}")
                continue
            drug_z = zscore(drug_score)

            # Joint model: pcr ~ drug + prolif + basal
            joint = fit_glm(y, ["DRUG", prolif_name, "BASAL"],
                            {"DRUG": drug_z, prolif_name: prolif_z, "BASAL": basal_z})
            # LR test: drug adds to baseline?
            lr_stat, lr_p = lr_test(joint["loglik"], baseline_results["loglik"], df_diff=1)
            delta_auc = joint["auroc"] - baseline_results["auroc"]

            entry = dict(
                cohort=cohort_name, drug_signature=drug_name,
                n=int(len(y)), pcr_rate=float(y.mean()),
                baseline_label=baseline_label,
                baseline_auroc=baseline_results["auroc"],
                joint_auroc=joint["auroc"],
                delta_auroc=float(delta_auc),
                drug_coef=joint["DRUG_coef"],
                drug_or_per_sd=joint["DRUG_or"],
                drug_wald_p=joint["DRUG_p"],
                prolif_coef=joint[f"{prolif_name}_coef"],
                prolif_p=joint[f"{prolif_name}_p"],
                basal_coef=joint["BASAL_coef"],
                basal_p=joint["BASAL_p"],
                lr_test_stat=lr_stat,
                lr_test_p=lr_p,
                joint_beats_baseline_03=bool(delta_auc >= 0.03 and joint["auroc"] > baseline_results["auroc"]),
                drug_significant_uncorrected=bool(joint["DRUG_p"] < 0.05),
                drug_significant_bonferroni_per_cohort=bool(joint["DRUG_p"] < 0.05 / len(DRUG_SIGS)),
            )
            results.append(entry)
            print(f"    {drug_name:<18}: joint_AUROC={joint['auroc']:.3f}  ΔAUROC={delta_auc:+.3f}  "
                  f"drug_OR={joint['DRUG_or']:.2f} p={joint['DRUG_p']:.3g}  LR_p={lr_p:.3g}")

    out_path = REVISE_ROOT / "primary_test_results" / "incremental_value.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    # Summary TSV
    df_out = pd.DataFrame(results)
    df_out.to_csv(out_path.with_suffix(".tsv"), sep="\t", index=False, float_format="%.4g")

    # Headline decision: does framework add value in ≥ 2 of 3 cohorts?
    print("\n" + "=" * 80)
    print("INCREMENTAL VALUE SUMMARY")
    print("=" * 80)
    for drug in DRUG_SIGS.keys():
        sub = df_out[df_out["drug_signature"] == drug]
        n_sig = sub["drug_significant_uncorrected"].sum()
        n_03 = sub["joint_beats_baseline_03"].sum()
        n_total = len(sub)
        print(f"  {drug:<20}: drug coef p<0.05 in {n_sig}/{n_total} cohorts; ΔAUROC≥+0.03 in {n_03}/{n_total}")

    print(f"\nResults written to {out_path}")


if __name__ == "__main__":
    main()
