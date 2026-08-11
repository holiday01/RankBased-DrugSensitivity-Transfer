#!/usr/bin/env python3
"""A0b -- re-score the 60-gene framework with the corrected bidirectional singscore
and re-run the headline inference on all four pCR cohorts.

A0_scorer_validation.py established on synthetic data that the deployed formula
(raw = up_score - down_score, with down_score already built on reversed ranks) is an
UNDIRECTED gene-set score. This script quantifies what that cost.

The four axes are NOT affected and are reused unchanged: P is a unidirectional Hallmark
singscore, H and B are transcript z-scores, Basal is a genefu PAM50 centroid call. Only
the bidirectional signature score is recomputed, so any change below is attributable to
the scorer fix alone.

Outputs: results/A0/ -- corrected scores per cohort, and a published-vs-corrected table.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

import os

# Paths are deposit-relative by default so this runs outside the authors' machine.
# DATA_ROOT: frozen_signatures/, external_data/ and the phase_* result tree.
# OUT_ROOT : where A0/ and A1/ outputs are written.
_DR = Path(os.environ.get("DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
_OR = Path(os.environ.get("OUT_ROOT", Path(__file__).resolve().parents[1] / "results"))

ROOT = _DR
OUT = _OR / "A0"
SIG = ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv"
AXES = ["P", "H", "B", "H_PAM50basal"]

COHORTS = {
    "GSE16446": dict(expr=ROOT / "external_data/GSE16446/GSE16446_expression.tsv",
                     axes=ROOT / "resubmission_v2/results/phase_0/axes_GSE16446.tsv"),
    "GSE25066": dict(expr=ROOT / "external_data/GSE25066/GSE25066_expression.tsv",
                     axes=ROOT / "resubmission_v2/results/phase_0/axes_GSE25066.tsv"),
    "GSE22226": dict(expr=ROOT / "external_data/GSE22226/GSE22226_expression.tsv",
                     axes=ROOT / "resubmission_v2/results/phase_0/axes_GSE22226.tsv"),
    "GSE32646": dict(expr=ROOT / "resubmission_v2/results/phase_2/GSE32646_expression.tsv",
                     axes=ROOT / "resubmission_v2/results/phase_2/axes_GSE32646.tsv"),
}


def zscore(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def singscore(expr, up, down, variant):
    """variant='deployed' reproduces the bug; variant='corrected' is Foroutan TotalScore."""
    n_total = expr.shape[1]
    if n_total < 5000:
        raise ValueError(f"only {n_total} genes")
    ranks = expr.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in expr.columns]
    D = [g for g in down if g in expr.columns]

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1

    up_s = norm(ranks[U].mean(axis=1), len(U))
    dn_s = norm((1 - ranks[D]).mean(axis=1), len(D))
    raw = (up_s - dn_s) if variant == "deployed" else (up_s + dn_s)
    return raw - raw.median(), len(U), len(D)


def analyse(y, Xr, sig):
    """Marginal AUROC, residual AUROC (OLS on axes), nested LRT, added OR/SD."""
    s = zscore(sig).values
    marginal = roc_auc_score(y, s)

    # residual-AUROC: regress signature on the axes (outcome never enters), test residual
    A = sm.add_constant(Xr)
    resid = s - A @ np.linalg.lstsq(A, s, rcond=None)[0]
    residual = roc_auc_score(y, resid)
    r2 = 1 - resid.var() / s.var()

    mr = sm.Logit(y, sm.add_constant(Xr)).fit(disp=0)
    mf = sm.Logit(y, sm.add_constant(np.column_stack([Xr, s]))).fit(disp=0)
    lr_p = float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1))
    beta = mf.params[-1]
    se = mf.bse[-1]
    return dict(marginal=marginal, residual=residual, r2_on_axes=r2, lr_p=lr_p,
                added_or_sd=float(np.exp(beta)),
                or_ci=[float(np.exp(beta - 1.96 * se)), float(np.exp(beta + 1.96 * se))],
                auc_reduced=roc_auc_score(y, mr.predict(sm.add_constant(Xr))),
                auc_full=roc_auc_score(y, mf.predict(sm.add_constant(np.column_stack([Xr, s])))))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.read_csv(SIG, sep="\t")
    up = sig.loc[sig.direction == "sensitivity", "gene_symbol"].tolist()
    down = sig.loc[sig.direction == "resistance", "gene_symbol"].tolist()

    rows = []
    for name, cfg in COHORTS.items():
        expr = pd.read_csv(cfg["expr"], sep="\t", index_col=0)
        if expr.shape[1] < 5000:            # stored genes-as-rows
            expr = expr.T
        ax = pd.read_csv(cfg["axes"], sep="\t", index_col=0)
        ax = ax.dropna(subset=AXES + ["pCR"])

        common = expr.index.intersection(ax.index)
        expr = expr.loc[common]
        ax = ax.loc[common]
        y = ax["pCR"].values.astype(int)
        Xr = np.column_stack([zscore(ax[a]).values for a in AXES])

        out = {}
        for variant in ("deployed", "corrected"):
            s, nU, nD = singscore(expr, up, down, variant)
            out[variant] = analyse(y, Xr, s.loc[common])
            out[variant].update(n_up=nU, n_down=nD)
            pd.DataFrame({"sample_id": common, f"score_{variant}": s.loc[common].values}) \
                .to_csv(OUT / f"score_{variant}_{name}.tsv", sep="\t", index=False)

        d, c = out["deployed"], out["corrected"]
        rows.append(dict(cohort=name, n=len(y), events=int(y.sum()),
                         genes=f"{d['n_up']}/{len(up)} up, {d['n_down']}/{len(down)} down",
                         **{f"dep_{k}": v for k, v in d.items() if k not in ("n_up", "n_down")},
                         **{f"cor_{k}": v for k, v in c.items() if k not in ("n_up", "n_down")}))

    json.dump(rows, open(OUT / "A0b_rescore_comparison.json", "w"), indent=2, default=float)

    print("A0b -- 60-gene framework, deployed (buggy) vs corrected bidirectional singscore\n")
    print(f"{'cohort':10s} {'n':>4s} {'ev':>3s} | {'marginal AUROC':>17s} | {'residual AUROC':>17s} | "
          f"{'R2 on axes':>13s} | {'nested-LRT p':>15s}")
    print(f"{'':10s} {'':>4s} {'':>3s} | {'dep':>8s} {'cor':>8s} | {'dep':>8s} {'cor':>8s} | "
          f"{'dep':>6s} {'cor':>6s} | {'dep':>7s} {'cor':>7s}")
    print("-" * 104)
    for r in rows:
        print(f"{r['cohort']:10s} {r['n']:4d} {r['events']:3d} | "
              f"{r['dep_marginal']:8.4f} {r['cor_marginal']:8.4f} | "
              f"{r['dep_residual']:8.4f} {r['cor_residual']:8.4f} | "
              f"{r['dep_r2_on_axes']:6.3f} {r['cor_r2_on_axes']:6.3f} | "
              f"{r['dep_lr_p']:7.4f} {r['cor_lr_p']:7.4f}")

    print(f"\n{'cohort':10s} | {'added OR/SD [95% CI], corrected':>40s} | crosses 0.55 residual?")
    print("-" * 90)
    for r in rows:
        ci = f"{r['cor_added_or_sd']:.3f} [{r['cor_or_ci'][0]:.3f}, {r['cor_or_ci'][1]:.3f}]"
        flag = "YES" if r["cor_residual"] >= 0.55 else "no"
        print(f"{r['cohort']:10s} | {ci:>40s} | {flag}")

    n_cross = sum(1 for r in rows if r["cor_residual"] >= 0.55)
    n_lrt = sum(1 for r in rows if r["cor_lr_p"] < 0.05)
    print(f"\nresidual AUROC >= 0.55 after correction : {n_cross} of {len(rows)} cohorts")
    print(f"nested-LRT p < 0.05 after correction     : {n_lrt} of {len(rows)} cohorts")
    print(f"\nwrote {OUT}/A0b_rescore_comparison.json")


if __name__ == "__main__":
    main()
