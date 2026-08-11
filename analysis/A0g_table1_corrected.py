#!/usr/bin/env python3
"""A0g -- regenerate Table 1 (60-gene framework, nested-model incremental value)
with the CORRECTED bidirectional singscore.

Table 1 columns that depend on the signature score, and so must be recomputed:
  LR p, added OR/SD [95% CI], dBIC, and the +sig (full-model) AUROC.
The four-axis (reduced-model) AUROC does not depend on the signature and is unchanged.
The pooled added OR/SD is a fixed-effect meta-analysis of the per-cohort log-ORs.

Every value is printed as deployed -> corrected so the change is explicit.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

DR = Path(os.environ.get("DATA_ROOT", "/home/holiday01/2026_ISMB_code/revise_bioadv"))
OUT = Path(os.environ.get("OUT_ROOT",
           "/home/holiday01/2026_ISMB_code/revise_bioadv_386/results")) / "A0"
AXES = ["P", "H", "B", "H_PAM50basal"]

SIG = DR / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv"
COH = {
    "GSE16446": (DR / "external_data/GSE16446/GSE16446_expression.tsv",
                 DR / "resubmission_v2/results/phase_0/axes_GSE16446.tsv", "pCR"),
    "GSE25066": (DR / "external_data/GSE25066/GSE25066_expression.tsv",
                 DR / "resubmission_v2/results/phase_0/axes_GSE25066.tsv", "pCR"),
    "GSE22226": (DR / "external_data/GSE22226/GSE22226_expression.tsv",
                 DR / "resubmission_v2/results/phase_0/axes_GSE22226.tsv", "pCR"),
    "GSE32646": (DR / "resubmission_v2/results/phase_2/GSE32646_expression.tsv",
                 DR / "resubmission_v2/results/phase_2/axes_GSE32646.tsv", "pCR"),
    "METABRIC": (DR / "resubmission_v2/results/phase_4/metabric/METABRIC_expression.tsv",
                 DR / "resubmission_v2/results/phase_4/metabric/axes_METABRIC.tsv", "OS5y"),
}
PCR = ["GSE16446", "GSE25066", "GSE22226", "GSE32646"]


def zs(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def sing(e, up, dn, corrected=True):
    n = e.shape[1]
    R = e.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in e.columns]
    D = [g for g in dn if g in e.columns]
    nm = lambda m, k: 2 * (m - (1 + k) / (2 * n)) / ((1 - (1 + k) / (2 * n)) - (1 + k) / (2 * n)) - 1
    a, b = nm(R[U].mean(1), len(U)), nm((1 - R[D]).mean(1), len(D))
    return (a + b) if corrected else (a - b)


def fit(y, X, s):
    s = zs(s).values
    mr = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
    mf = sm.Logit(y, sm.add_constant(np.column_stack([X, s]))).fit(disp=0)
    beta, se = mf.params[-1], mf.bse[-1]
    lrp = float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1))
    k_r, k_f, n = mr.df_model + 1, mf.df_model + 1, len(y)
    bic_r = -2 * mr.llf + k_r * np.log(n)
    bic_f = -2 * mf.llf + k_f * np.log(n)
    return dict(lr_p=lrp, logor=float(beta), se=float(se),
                or_sd=float(np.exp(beta)),
                or_ci=[float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se))],
                dbic=float(bic_f - bic_r),
                auc_red=roc_auc_score(y, mr.predict(sm.add_constant(X))),
                auc_full=roc_auc_score(y, mf.predict(sm.add_constant(np.column_stack([X, s])))))


def main():
    sig = pd.read_csv(SIG, sep="\t")
    up = sig[sig.direction == "sensitivity"].gene_symbol.tolist()
    dn = sig[sig.direction == "resistance"].gene_symbol.tolist()

    rows = {}
    for c, (ep, ap, endpt) in COH.items():
        e = pd.read_csv(ep, sep="\t", index_col=0)
        if e.shape[1] < 5000:
            e = e.T
        ax = pd.read_csv(ap, sep="\t", index_col=0)
        ycol = "pCR" if "pCR" in ax.columns else [x for x in ax.columns if "os" in x.lower() or "event" in x.lower() or "dead" in x.lower()][0]
        ax = ax.dropna(subset=[a for a in AXES if a in ax.columns] + [ycol])
        cm = e.index.intersection(ax.index)
        e, ax = e.loc[cm], ax.loc[cm]
        y = ax[ycol].values.astype(int)
        X = np.column_stack([zs(ax[a]).values for a in AXES])
        out = {}
        for tag, corr in [("dep", False), ("cor", True)]:
            out[tag] = fit(y, X, sing(e, up, dn, corr))
        rows[c] = dict(n=len(y), events=int(y.sum()), endpoint=endpt, **out)

    def pool(tag):
        w = np.array([1 / rows[c][tag]["se"] ** 2 for c in PCR])
        b = np.array([rows[c][tag]["logor"] for c in PCR])
        lb = float((w * b).sum() / w.sum())
        se = float(np.sqrt(1 / w.sum()))
        z = lb / se
        return dict(or_sd=float(np.exp(lb)),
                    or_ci=[float(np.exp(lb - 1.96 * se)), float(np.exp(lb + 1.96 * se))],
                    p=float(2 * stats.norm.sf(abs(z))))

    print("Table 1 -- 60-gene framework, deployed -> CORRECTED scorer\n")
    print(f"{'cohort':11s} {'n':>4s} {'ev':>3s} | {'LR p':>13s} | {'added OR/SD':>26s} | "
          f"{'dBIC':>11s} | {'AUROC 4ax/+sig':>19s}")
    print("-" * 104)
    for c in list(PCR) + ["METABRIC"]:
        r = rows[c]
        d, co = r["dep"], r["cor"]
        print(f"{c:11s} {r['n']:4d} {r['events']:3d} | {d['lr_p']:5.2f} -> {co['lr_p']:5.2f} | "
              f"{d['or_sd']:.2f}[{d['or_ci'][0]:.2f},{d['or_ci'][1]:.2f}] -> "
              f"{co['or_sd']:.2f}[{co['or_ci'][0]:.2f},{co['or_ci'][1]:.2f}] | "
              f"{d['dbic']:+.1f}->{co['dbic']:+.1f} | "
              f"{co['auc_red']:.3f}/{d['auc_full']:.3f}->{co['auc_full']:.3f}")

    pd_, pc = pool("dep"), pool("cor")
    print(f"\nPooled (4 pCR): added OR/SD  {pd_['or_sd']:.2f} [{pd_['or_ci'][0]:.2f},{pd_['or_ci'][1]:.2f}] "
          f"P={pd_['p']:.2f}  ->  {pc['or_sd']:.2f} [{pc['or_ci'][0]:.2f},{pc['or_ci'][1]:.2f}] P={pc['p']:.2f}")

    print("\n== corrected Table 1 rows (for the manuscript) ==")
    for c in list(PCR) + ["METABRIC"]:
        co = rows[c]["cor"]
        print(f"  {c:11s} LRp={co['lr_p']:.2f}  OR/SD={co['or_sd']:.2f} [{co['or_ci'][0]:.2f},{co['or_ci'][1]:.2f}]  "
              f"dBIC={co['dbic']:+.1f}  AUROC {co['auc_red']:.3f}/{co['auc_full']:.3f}")
    print(f"  {'Pooled':11s} OR/SD={pc['or_sd']:.2f} [{pc['or_ci'][0]:.2f},{pc['or_ci'][1]:.2f}] P={pc['p']:.2f}")

    import json
    json.dump({c: rows[c]["cor"] for c in rows} | {"pooled_cor": pc, "pooled_dep": pd_},
              open(OUT / "A0g_table1_corrected.json", "w"), indent=2, default=float)
    print(f"\nAll LR p > 0.05 (corrected): {all(rows[c]['cor']['lr_p'] > 0.05 for c in PCR)}")


if __name__ == "__main__":
    main()
