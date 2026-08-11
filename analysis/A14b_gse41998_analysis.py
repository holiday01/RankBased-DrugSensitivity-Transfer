#!/usr/bin/env python3
"""A14 stage b -- the added-value test on GSE41998.

Refuses to run until the section 7 gate verdict exists in decisions_addendum.md. The gate is
outcome-blind and one-way; this file holds the outcome-dependent code that stage a deliberately
does not contain.

GSE41998 is an OPTIONAL cohort under addendum section 7. It is **not** a member of Family A,
whose 58 cells are fixed in section 4.2 and may not be enlarged after the fact. Its cells are
therefore reported with nominal p, outside every BH family, exactly like A8 and the
wrong-mechanism arm.

Why it matters: Reviewer 2's major point M3 is that only GSE32646 is a genuinely independent
held-out pCR cohort, so the pooled estimate is partly circular. A second independent cohort with
253 pCR calls and zero patient overlap answers that premise rather than conceding it.

GSE41998 also carries a randomised within-cohort treatment contrast -- AC followed by
ixabepilone versus paclitaxel -- which no other cohort in the study has except GSE22226's
11-patient trastuzumab group. The signature is tested within each arm as a secondary analysis.

CIs are patient-level bootstrap percentile intervals, B = 2000, seed 42, matching A0i.
Writes results/A14/A14b_analysis.json.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import A0e_dedup as D  # noqa: E402
from A2_grid import build_defs, score, MIN_EVENTS, MIN_RESOLVED  # noqa: E402
from A0i_cluster_bootstrap import _cell_builders  # noqa: E402

B, SEED = 2000, 42
G = D.ROOT / "external_data/GSE41998"
OUT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A14")
LEDGER = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/decisions_addendum.md")
AXES = D.AXES


def require_gate() -> dict:
    if "GSE41998 gate" not in LEDGER.read_text():
        raise SystemExit("REFUSING TO RUN: the section 7 gate verdict is not in "
                         "decisions_addendum.md. Run A14a and record the verdict first.")
    g = json.loads((OUT / "A14a_gate.json").read_text())
    if g["verdict"] != "PASS":
        raise SystemExit(f"gate verdict is {g['verdict']}; the cohort is abandoned and leaves "
                         "no trace in the manuscript (section 7).")
    print(f"gate verdict {g['verdict']} recorded at {g['utc']} -- proceeding")
    return g


def cell(y, Xr, s):
    A = sm.add_constant(Xr)
    R, F = A, np.column_stack([A, s])
    mr, mf = sm.Logit(y, R).fit(disp=0), sm.Logit(y, F).fit(disp=0)
    chi2 = 2 * (mf.llf - mr.llf)
    b, se = mf.params[-1], mf.bse[-1]
    fit = A @ np.linalg.lstsq(A, s, rcond=None)[0]
    resid = roc_auc_score(y, s - fit)
    return dict(n=int(len(y)), events=int(y.sum()),
                marginal=float(roc_auc_score(y, s)), residual=float(resid),
                r2_on_axes=float(1 - ((s - fit) ** 2).sum() / ((s - s.mean()) ** 2).sum()),
                lr_chi2=float(chi2), lr_p=float(stats.chi2.sf(chi2, 1)),
                added_or_sd=float(np.exp(b)),
                or_ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                dbic=float(np.log(len(y)) - chi2),
                auc_reduced=float(roc_auc_score(y, mr.predict(R))),
                auc_full=float(roc_auc_score(y, mf.predict(F))),
                converged=bool(mf.mle_retvals["converged"]))


def boot_ci(y, Xr, s):
    A = sm.add_constant(Xr)
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(B):
        i = rng.integers(0, len(y), len(y))
        yi = y[i]
        if min(yi.sum(), len(yi) - yi.sum()) < MIN_EVENTS:
            continue
        Ai, si = A[i], s[i]
        vals.append(roc_auc_score(yi, si - Ai @ np.linalg.lstsq(Ai, si, rcond=None)[0]))
    v = np.array(vals)
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))], float(v.std(ddof=1))


def verdict_of(c, ci):
    lo, hi = ci
    if hi < 0.5:
        return "INVERSE-DIRECTION"
    if c["lr_p"] < 0.05 and lo > 0.5:
        return "ADDS-VALUE (nominal; outside Family A, no BH applied)"
    if (c["marginal"] >= 0.60 and c["residual"] <= 0.55 and lo <= 0.5 <= hi
            and c["r2_on_axes"] >= 0.50):
        return "REDUNDANT-WITH-AXES"
    if c["marginal"] < 0.60:
        return "NO-MARGINAL-SIGNAL"
    return "INDETERMINATE"


def main() -> None:
    gate = require_gate()
    ax = pd.read_csv(G / "GSE41998_axes.tsv", sep="\t", index_col=0)
    expr = pd.read_csv(G / "GSE41998_expression.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(G / "GSE41998_phenotype.tsv", sep="\t", index_col=0)
    common = ax.index.intersection(expr.index)
    ax, expr = ax.loc[common], expr.loc[common]
    y = ax["pCR"].values.astype(int)
    Xr = np.column_stack([D.zscore(ax[a]).values for a in AXES])
    print(f"analysis set n={len(y)} events={int(y.sum())} "
          f"(base four-axis AUROC computed below)")

    src = _cell_builders()
    defs = build_defs(src)
    res = {"cohort": "GSE41998", "family": "outside Family A (section 7 optional cohort); "
                                            "nominal p, no BH", "gate": gate, "cells": {}}
    for name in list(defs)[:10]:
        up, dn, scorer, cohdep, fid, restrict = defs[name]
        sel = np.ones(len(ax), dtype=bool)
        if restrict == "HR+":
            sel = (ph["er_status_ihc"].reindex(ax.index).astype(str)
                   .str.lower().eq("positive")).values
            if sel.sum() < 30:
                res["cells"][name] = dict(not_evaluable=True,
                                          clause="too few HR-positive patients")
                continue
        s, nu, nd = score(expr.loc[sel], up, dn, scorer)
        s = D.zscore(s).values
        fu = nu / max(len(up), 1)
        fd = (nd / len(dn)) if dn else None
        if fu < MIN_RESOLVED or (fd is not None and fd < MIN_RESOLVED):
            res["cells"][name] = dict(not_evaluable=True, clause="gene resolution below 0.60",
                                      frac_up=fu, frac_down=fd)
            continue
        c = cell(y[sel], Xr[sel], s)
        ci, se = boot_ci(y[sel], Xr[sel], s)
        c.update(residual_ci=ci, residual_se=se, scorer=scorer, deployment_fidelity=fid,
                 frac_up_resolved=float(fu),
                 frac_down_resolved=(float(fd) if fd is not None else None))
        c["verdict"] = verdict_of(c, ci)
        res["cells"][name] = c
        print(f"  {name:12s} n={c['n']:3d} marg={c['marginal']:.3f} "
              f"resid={c['residual']:.3f} [{ci[0]:.3f},{ci[1]:.3f}] "
              f"R2ax={c['r2_on_axes']:.3f} p={c['lr_p']:.4f}  {c['verdict']}")

    # ---- randomised within-cohort treatment contrast ----
    arm = ph["treatment_arm"].reindex(ax.index).astype(str).str.lower()
    fw_s = D.zscore(score(expr, *defs["Framework"][:2], defs["Framework"][2])[0]).values
    res["by_treatment_arm"] = {}
    for a in ["ixabepilone", "paclitaxel"]:
        m = arm.eq(a).values
        if m.sum() >= 40 and min(y[m].sum(), m.sum() - y[m].sum()) >= MIN_EVENTS:
            r = cell(y[m], Xr[m], fw_s[m])
            res["by_treatment_arm"][a] = r
            print(f"  arm {a:12s} n={r['n']:3d} events={r['events']:3d} "
                  f"resid={r['residual']:.3f} p={r['lr_p']:.4f}")
    res["arm_note"] = ("GSE41998 randomises AC followed by ixabepilone versus paclitaxel. This "
                       "is the only randomised within-cohort treatment contrast available to "
                       "the study, so it is the only place where regimen can be varied with "
                       "everything else held constant.")

    fw = res["cells"]["Framework"]
    res["headline"] = (
        f"On a second genuinely independent neoadjuvant pCR cohort (n = {fw['n']}, "
        f"{fw['events']} events, zero patient overlap with any discovery cohort), the "
        f"cell-line-derived signature returns residual AUROC {fw['residual']:.3f} "
        f"[{fw['residual_ci'][0]:.3f}, {fw['residual_ci'][1]:.3f}] and nested-LRT "
        f"p = {fw['lr_p']:.4f} against the four axes, whose own AUROC is "
        f"{fw['auc_reduced']:.3f}.")
    print("\n" + res["headline"])
    (OUT / "A14b_analysis.json").write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {OUT / 'A14b_analysis.json'}")


if __name__ == "__main__":
    main()
