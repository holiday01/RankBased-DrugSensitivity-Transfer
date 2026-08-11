#!/usr/bin/env python3
"""A0h -- regenerate the 5 bidirectional-scorer rows of Table 2 (the signature grid)
with the CORRECTED singscore, each on the cohort/stratum Table 2 assigns to it.

Rows 6-10 (TIS, cytolytic, TGFBstrom, HRD, DLDA30) use z-score scorers unaffected by the
bug and are unchanged; this script does not touch them.

Per row: marginal AUROC, residual AUROC with bootstrap 95% CI, nested-LRT p.
"""
import ast
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

DR = Path(os.environ.get("DATA_ROOT", "/home/holiday01/2026_ISMB_code/revise_bioadv"))
SCR = DR / "resubmission_v2/scripts"
AXES = ["P", "H", "B", "H_PAM50basal"]


def zs(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def sing(e, up, dn):
    n = e.shape[1]
    R = e.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in e.columns]
    D = [g for g in dn if g in e.columns]
    nm = lambda m, k: 2 * (m - (1 + k) / (2 * n)) / ((1 - (1 + k) / (2 * n)) - (1 + k) / (2 * n)) - 1
    return nm(R[U].mean(1), len(U)) + nm((1 - R[D]).mean(1), len(D))


def lists(path, *names):
    out = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and (not names or t.id in names):
                    try:
                        v = ast.literal_eval(node.value)
                        if isinstance(v, (list, dict)) and v:
                            out[t.id] = v
                    except Exception:
                        pass
    return out


EXPR = {"GSE16446": DR / "external_data/GSE16446/GSE16446_expression.tsv",
        "GSE22226": DR / "external_data/GSE22226/GSE22226_expression.tsv",
        "GSE25066": DR / "external_data/GSE25066/GSE25066_expression.tsv",
        "GSE32646": DR / "resubmission_v2/results/phase_2/GSE32646_expression.tsv",
        "METABRIC": DR / "resubmission_v2/results/phase_4/metabric/METABRIC_expression.tsv"}
AXF = {"GSE16446": DR / "resubmission_v2/results/phase_0/axes_GSE16446.tsv",
       "GSE22226": DR / "resubmission_v2/results/phase_0/axes_GSE22226.tsv",
       "GSE25066": DR / "resubmission_v2/results/phase_0/axes_GSE25066.tsv",
       "GSE32646": DR / "resubmission_v2/results/phase_2/axes_GSE32646.tsv",
       "METABRIC": DR / "resubmission_v2/results/phase_4/metabric/axes_METABRIC.tsv"}


def load(coh):
    e = pd.read_csv(EXPR[coh], sep="\t", index_col=0)
    if e.shape[1] < 5000:
        e = e.T
    ax = pd.read_csv(AXF[coh], sep="\t", index_col=0)
    yc = "pCR" if "pCR" in ax.columns else [c for c in ax.columns
          if any(k in c.lower() for k in ("os", "event", "dead", "pcr"))][0]
    ax = ax.dropna(subset=[a for a in AXES if a in ax.columns] + [yc])
    cm = e.index.intersection(ax.index)
    return e.loc[cm], ax.loc[cm], yc


def main():
    fw = pd.read_csv(DR / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv", sep="\t")
    dv = pd.read_csv(DR / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    ggi = lists(SCR / "phase_6/gge_GGI_gse25066.py")
    gu = [k for k in ggi if "UP" in k.upper() and "DOWN" not in k.upper()][0]
    gd = [k for k in ggi if "DOWN" in k.upper()][0]
    rs = lists(SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")
    mp = lists(SCR / "phase_6/mammaprint_metabric.py", "MAMMAPRINT")["MAMMAPRINT"]

    SIGS = {
        "framework":  (fw[fw.direction == "sensitivity"].gene_symbol.tolist(),
                       fw[fw.direction == "resistance"].gene_symbol.tolist(), "GSE32646", None),
        "GGI":        (ggi[gu], ggi[gd], "GSE25066", None),
        "RSlike":     (rs["RS_UP"], rs["RS_DOWN"], "GSE25066", "HR+"),
        "DoxVariant": (dv[dv.direction == "sensitivity"].gene_symbol.tolist(),
                       dv[dv.direction == "resistance"].gene_symbol.tolist(), "GSE25066", None),
        "MammaPrint": ([g for g, v in mp.items() if v > 0],
                       [g for g, v in mp.items() if v < 0], "METABRIC", None),
    }
    hrpos = set(pd.read_csv(DR / "external_data/GSE25066/GSE25066_pcr_HR+.tsv", sep="\t")["patient_id"])

    print("Table 2 -- 5 bidirectional rows, corrected scorer, each on its assigned cohort\n")
    print(f"{'signature':12s} {'cohort':14s} {'n':>4s} | {'marginal':>8s} | "
          f"{'residual [95% CI]':>22s} | {'LR p':>6s}")
    print("-" * 78)
    res = {}
    for name, (up, dn, coh, strat) in SIGS.items():
        e, ax, yc = load(coh)
        if strat == "HR+":
            keep = ax.index.intersection(hrpos)
            e, ax = e.loc[keep], ax.loc[keep]
        y = ax[yc].values.astype(int)
        s = zs(sing(e, up, dn).loc[ax.index]).values
        X = np.column_stack([zs(ax[a]).values for a in AXES])
        A = sm.add_constant(X)
        marg = roc_auc_score(y, s)
        resid = s - A @ np.linalg.lstsq(A, s, rcond=None)[0]
        pt = roc_auc_score(y, resid)
        rng = np.random.default_rng(42)
        bs = [roc_auc_score(y[i], resid[i]) for i in
              (rng.integers(0, len(y), len(y)) for _ in range(1000)) if len(np.unique(y[i])) > 1]
        lo, hi = np.percentile(bs, 2.5), np.percentile(bs, 97.5)
        mr = sm.Logit(y, A).fit(disp=0)
        mf = sm.Logit(y, sm.add_constant(np.column_stack([X, s]))).fit(disp=0)
        lrp = float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1))
        res[name] = dict(cohort=coh + (" HR+" if strat else ""), n=len(y),
                         marginal=round(marg, 3), residual=round(pt, 3),
                         ci=[round(lo, 3), round(hi, 3)], lr_p=round(lrp, 3))
        cl = coh + (" HR+" if strat else "")
        print(f"{name:12s} {cl:14s} {len(y):4d} | {marg:8.3f} | "
              f"{pt:.3f} [{lo:.3f}, {hi:.3f}]  | {lrp:6.3f}")

    import json
    json.dump(res, open(Path(os.environ.get("OUT_ROOT",
              "/home/holiday01/2026_ISMB_code/revise_bioadv_386/results")) / "A0" /
              "A0h_table2_corrected.json", "w"), indent=2)
    print("\nall LR p > 0.05:", all(r["lr_p"] > 0.05 for r in res.values()))


if __name__ == "__main__":
    main()
