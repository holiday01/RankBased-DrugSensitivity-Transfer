#!/usr/bin/env python3
"""A4 -- does the signature predict within subtype, or only separate subtypes?

Reviewer 3 point 4 asks a factual question that wording cannot answer: does the
signature carry pCR information inside a subtype, or does its marginal discrimination
come entirely from separating subtypes whose pCR rates already differ?

Method: the same nested-model LRT used for the headline result, refitted inside each
IHC-defined stratum (HR+/HER2-, HER2+, TNBC). Within a stratum the subtype contrast is
constant, so any remaining discrimination is within-subtype information.

Scores use the CORRECTED bidirectional singscore established in A0. Strata with fewer
than 5 events in either outcome class are reported as not evaluable, following the parent
pre-registration's own rule.

Reports per cell: n, events, pCR rate, marginal AUROC with bootstrap CI, added OR/SD,
nested-LRT p, and the marginal-minus-within-stratum AUROC gap -- the last being the
quantity that answers the reviewer directly.
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

_DR = Path(os.environ.get("DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
_OR = Path(os.environ.get("OUT_ROOT", Path(__file__).resolve().parents[1] / "results"))
OUT = _OR / "A4"
AXES = ["P", "H", "B", "H_PAM50basal"]
NBOOT, SEED = 1000, 42

COHORTS = {
    "GSE16446": (_DR / "external_data/GSE16446/GSE16446_expression.tsv",
                 _DR / "resubmission_v2/results/phase_0/axes_GSE16446.tsv",
                 _DR / "external_data/GSE16446/GSE16446_pcr.tsv"),
    "GSE22226": (_DR / "external_data/GSE22226/GSE22226_expression.tsv",
                 _DR / "resubmission_v2/results/phase_0/axes_GSE22226.tsv",
                 _DR / "external_data/GSE22226/GSE22226_pcr.tsv"),
    "GSE25066": (_DR / "external_data/GSE25066/GSE25066_expression.tsv",
                 _DR / "resubmission_v2/results/phase_0/axes_GSE25066.tsv",
                 _DR / "external_data/GSE25066/GSE25066_pcr.tsv"),
}


def zscore(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def singscore(expr, up, down):
    """Corrected bidirectional singscore (A0): TotalScore = up + down_reversed."""
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in expr.columns]
    D = [g for g in down if g in expr.columns]

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1

    raw = norm(ranks[U].mean(axis=1), len(U)) + norm((1 - ranks[D]).mean(axis=1), len(D))
    return raw - raw.median()


def strata_of(pcr):
    """IHC strata. HER2+ takes precedence, then TNBC, then HR+/HER2-."""
    er = pcr["er_status"].astype(str).str.lower()
    her2 = pcr["her2_status"].astype(str).str.lower()
    pos = lambda x: x.str.startswith(("p", "1", "y", "t"))
    out = pd.Series("other", index=pcr.index)
    out[pos(er) & ~pos(her2)] = "HR+/HER2-"
    out[pos(her2)] = "HER2+"
    out[~pos(er) & ~pos(her2)] = "TNBC"
    return out


def analyse(y, Xr, s, nboot=NBOOT, seed=SEED):
    s = zscore(s).values
    marg = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nboot):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) > 1:
            bs.append(roc_auc_score(y[i], s[i]))
    ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))] if bs else [np.nan] * 2
    mr = sm.Logit(y, sm.add_constant(Xr)).fit(disp=0)
    mf = sm.Logit(y, sm.add_constant(np.column_stack([Xr, s]))).fit(disp=0)
    b, se = mf.params[-1], mf.bse[-1]
    return dict(marginal=round(marg, 4), marginal_ci=[round(c, 4) for c in ci],
                or_sd=round(float(np.exp(b)), 3),
                or_ci=[round(float(np.exp(b - 1.96 * se)), 3),
                       round(float(np.exp(b + 1.96 * se)), 3)],
                lr_p=round(float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1)), 4))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.read_csv(_DR / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv",
                      sep="\t")
    up = sig.loc[sig.direction == "sensitivity", "gene_symbol"].tolist()
    down = sig.loc[sig.direction == "resistance", "gene_symbol"].tolist()

    rows = []
    for coh, (ep, ap, pp) in COHORTS.items():
        expr = pd.read_csv(ep, sep="\t", index_col=0)
        if expr.shape[1] < 5000:
            expr = expr.T
        ax = pd.read_csv(ap, sep="\t", index_col=0).dropna(subset=AXES + ["pCR"])
        pcr = pd.read_csv(pp, sep="\t").set_index("patient_id")
        common = expr.index.intersection(ax.index).intersection(pcr.index)
        expr, ax, pcr = expr.loc[common], ax.loc[common], pcr.loc[common]

        s_all = singscore(expr, up, down)
        strat = strata_of(pcr)

        for name, mask in [("ALL (marginal)", pd.Series(True, index=common))] + \
                          [(k, strat == k) for k in ["HR+/HER2-", "HER2+", "TNBC"]]:
            sub = common[mask.values]
            if len(sub) < 10:
                continue
            y = ax.loc[sub, "pCR"].values.astype(int)
            ev, ne = int(y.sum()), int((1 - y).sum())
            if min(ev, ne) < 5:
                rows.append(dict(cohort=coh, stratum=name, n=len(y), events=ev,
                                 pcr_rate=round(ev / len(y), 3), evaluable=False,
                                 note="<5 events in one class"))
                continue
            Xr = np.column_stack([zscore(ax.loc[sub, a]).values for a in AXES])
            r = analyse(y, Xr, s_all.loc[sub])
            rows.append(dict(cohort=coh, stratum=name, n=len(y), events=ev,
                             pcr_rate=round(ev / len(y), 3), evaluable=True, **r))

    json.dump(rows, open(OUT / "A4_within_subtype.json", "w"), indent=2, default=float)

    print("A4 -- nested-model LRT within IHC subtype (corrected scorer)\n")
    print(f"{'cohort':10s} {'stratum':13s} {'n':>4s} {'ev':>3s} {'pCR':>6s} "
          f"{'marginal AUROC':>22s} {'added OR/SD':>20s} {'LR p':>7s}")
    print("-" * 96)
    for r in rows:
        if not r.get("evaluable"):
            print(f"{r['cohort']:10s} {r['stratum']:13s} {r['n']:4d} {r['events']:3d} "
                  f"{r['pcr_rate']:6.3f}   NOT EVALUABLE ({r['note']})")
            continue
        ci = f"{r['marginal']:.3f} [{r['marginal_ci'][0]:.3f},{r['marginal_ci'][1]:.3f}]"
        orr = f"{r['or_sd']:.2f} [{r['or_ci'][0]:.2f},{r['or_ci'][1]:.2f}]"
        print(f"{r['cohort']:10s} {r['stratum']:13s} {r['n']:4d} {r['events']:3d} "
              f"{r['pcr_rate']:6.3f} {ci:>22s} {orr:>20s} {r['lr_p']:7.4f}")

    print("\nmarginal minus within-stratum AUROC (the quantity R3-4 asks about):")
    for coh in COHORTS:
        cr = [r for r in rows if r["cohort"] == coh and r.get("evaluable")]
        allr = next((r for r in cr if r["stratum"] == "ALL (marginal)"), None)
        if not allr:
            continue
        for r in cr:
            if r["stratum"] == "ALL (marginal)":
                continue
            print(f"  {coh:10s} {r['stratum']:12s} {allr['marginal']:.3f} -> {r['marginal']:.3f}"
                  f"   gap {allr['marginal'] - r['marginal']:+.3f}")
    sig_cells = [r for r in rows if r.get("evaluable") and r["stratum"] != "ALL (marginal)"
                 and r["lr_p"] < 0.05]
    print(f"\nwithin-stratum cells with LR p < 0.05 (nominal): {len(sig_cells)}")
    for r in sig_cells:
        print(f"  {r['cohort']} {r['stratum']} p={r['lr_p']}")
    print(f"\nwrote {OUT}/A4_within_subtype.json")


if __name__ == "__main__":
    main()
