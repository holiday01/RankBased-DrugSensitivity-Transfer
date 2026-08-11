#!/usr/bin/env python3
"""A5 and A4, extended to the cohorts they were missing.

Addendum §13.11 (A4) and §6/A5. Two gaps the 2026-08-01 audit found, both closed here:

A5  Reviewer 3 point 3 asks for subtype counts in **each** cohort. The deposited A5 covers
    GSE16446, GSE22226 and GSE25066 -- three of five. GSE32646 and METABRIC are added.

A4  Reviewer 3 point 4 asks whether the signature predicts within subtype. The deposited A4
    covers three discovery cohorts and **excludes GSE32646**, the only genuinely independent
    cohort in the study, whose IHC is 100% complete. It also reports no per-stratum EPV and no
    MDE, so "no cell reaches nominal significance" cannot be read as evidence of absence, and
    its GSE16446 "TNBC" row is byte-identical to that cohort's marginal row -- the cohort is
    uniformly triple-negative, so that is not a stratified analysis and it is dropped.

Corrections applied per §13.11: GSE32646 added, the pseudo-stratum removed, per-stratum events,
EPV, MDE and residual AUROC reported, and the §4.1 de-duplication applied to GSE25066.

MDE is obtained by the same inversion used in A10 -- the smallest added OR/SD detectable at 80%
power given that stratum's n and event count -- and is a power calculation, not a measurement.

Writes results/A5/A5_cohort_characteristics_extended.json and
results/A4/A4_within_subtype_extended.json.
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
from A0i_cluster_bootstrap import _cell_builders  # noqa: E402
from A3_baselines import ihc, PHEN  # noqa: E402
import A4_within_subtype as A4  # noqa: E402  -- reuse the DEPOSITED stratum definition verbatim

RES = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results")
MB = D.ROOT / "resubmission_v2/results/phase_4/metabric"
SEED = 42
MIN_EVENTS = 5


def strata_from_pcr(cohort, index):
    """The DEPOSITED stratum definition, imported rather than re-implemented.

    A first draft of this script wrote its own: it used PR as well as ER, and dropped patients
    with missing HER2 instead of treating them as negative. That is a different partition -- it
    moved GSE25066's TNBC stratum from 199 to 148 -- and under it two cells reached p < 0.05
    where the deposited analysis reports none. Changing a stratum definition and reporting the
    resulting significance would be exactly the researcher degree of freedom this paper is
    about. §13.11 authorises adding GSE32646, deleting the pseudo-stratum, reporting EPV/MDE
    and applying the de-duplication -- nothing else. A4.strata_of is used verbatim.
    """
    pcr = pd.read_csv(PCR_PATH[cohort], sep="\t", dtype=str).set_index("patient_id")
    return A4.strata_of(pcr.reindex(index))


PCR_PATH = {
    "GSE16446": D.ROOT / "external_data/GSE16446/GSE16446_pcr.tsv",
    "GSE22226": D.ROOT / "external_data/GSE22226/GSE22226_pcr.tsv",
    "GSE25066": D.ROOT / "external_data/GSE25066/GSE25066_pcr.tsv",
    "GSE32646": D.ROOT / "resubmission_v2/results/phase_2/GSE32646_pcr.tsv",
}


def mde_or(n, events, reps=300, seed=SEED):
    """Smallest added OR/SD detectable at 80% power in a 5-predictor logistic model at this
    stratum's n and event count. Power calculation, not a measurement."""
    if events < MIN_EVENTS or n - events < MIN_EVENTS:
        return None
    rng = np.random.default_rng(seed)
    base = np.log(events / (n - events))

    def power(beta):
        hit = 0
        for _ in range(reps):
            X = rng.normal(size=(n, 4))
            s = rng.normal(size=n)
            lin = 0.6 * X[:, 0] - 0.4 * X[:, 1] + beta * s
            lin += base - lin.mean()
            y = (rng.uniform(size=n) < 1 / (1 + np.exp(-lin))).astype(int)
            if min(y.sum(), n - y.sum()) < MIN_EVENTS:
                continue
            R = sm.add_constant(X)
            try:
                mr = sm.Logit(y, R).fit(disp=0)
                mf = sm.Logit(y, np.column_stack([R, s])).fit(disp=0)
                if stats.chi2.sf(2 * (mf.llf - mr.llf), 1) < 0.05:
                    hit += 1
            except Exception:
                continue
        return hit / reps

    lo, hi = 0.05, 3.0
    for _ in range(9):
        mid = (lo + hi) / 2
        if power(mid) >= 0.80:
            hi = mid
        else:
            lo = mid
    return round(float(np.exp(hi)), 3)


def main() -> None:
    _, b_ov, _, _, _ = D.overlap()
    src = _cell_builders()

    # ---------------- A5: subtype counts, all five cohorts ----------------
    a5 = json.loads((RES / "A5/A5_cohort_characteristics.json").read_text())
    have = {r["cohort"] for r in a5}
    print("A5 already covers:", sorted(have))

    for coh in ["GSE32646"]:
        expr, ax = src["load"](coh)
        t = ihc(coh, ax.index)
        er, pr, her2 = t.get("ER"), t.get("PR"), t.get("HER2")
        st = strata_from_pcr(coh, ax.index)
        y = ax["pCR"].values.astype(int)
        a5.append(dict(cohort=coh, n_profiled=int(len(ax)), n_analysed=int(len(ax)),
                       events=int(y.sum()), pcr_rate=round(float(y.mean()), 3),
                       ER_pos=int((er == 1).sum()), ER_neg=int((er == 0).sum()),
                       HER2_pos=int((her2 == 1).sum()), HER2_neg=int((her2 == 0).sum()),
                       TNBC=int((st == "TNBC").sum()),
                       HRpos_HER2neg=int((st == "HR+/HER2-").sum()),
                       HER2pos=int((st == "HER2+").sum())))
        print(f"  + {coh}: n={len(ax)} events={int(y.sum())} "
              f"TNBC={int((st=='TNBC').sum())} HER2+={int((st=='HER2+').sum())} "
              f"HR+/HER2-={int((st=='HR+/HER2-').sum())}")

    mb_ax = pd.read_csv(MB / "axes_METABRIC.tsv", sep="\t", index_col=0)
    mb_cl = pd.read_csv(MB / "METABRIC_clinical.tsv", sep="\t", dtype=str)
    idc = "SAMPLE_ID" if "SAMPLE_ID" in mb_cl.columns else mb_cl.columns[0]
    mb_cl = mb_cl.set_index(idc).reindex(mb_ax.index)
    er = mb_cl.get("ER_STATUS", pd.Series(dtype=str)).astype(str).str.upper()
    her2 = mb_cl.get("HER2_STATUS", pd.Series(dtype=str)).astype(str).str.upper()
    ev = mb_ax["event5y"].dropna()
    a5.append(dict(cohort="METABRIC", n_profiled=int(len(mb_ax)), n_analysed=int(len(ev)),
                   events=int(ev.sum()), pcr_rate=None,
                   endpoint="dichotomised 5-year OS, not pCR",
                   ER_pos=int((er == "POSITIVE").sum()), ER_neg=int((er == "NEGATIVE").sum()),
                   HER2_pos=int((her2 == "POSITIVE").sum()),
                   HER2_neg=int((her2 == "NEGATIVE").sum())))
    print(f"  + METABRIC: n={len(mb_ax)} 5y-events={int(ev.sum())} "
          f"ER+={int((er=='POSITIVE').sum())} HER2+={int((her2=='POSITIVE').sum())}")
    (RES / "A5/A5_cohort_characteristics_extended.json").write_text(json.dumps(a5, indent=2))

    # ---------------- A4: within-subtype, four cohorts, with EPV and MDE ----------------
    rows = []
    for coh in ["GSE16446", "GSE22226", "GSE25066", "GSE32646"]:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        t = ihc(coh, ax.index)
        st = strata_from_pcr(coh, ax.index)
        uniform = st[st != "other"].nunique() <= 1
        y = ax["pCR"].values.astype(int)
        Xr = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
        s = D.zscore(src["sig_of"]("Framework", expr)[0]).values
        if uniform:
            print(f"  {coh}: uniformly one stratum -- no within-subtype analysis is possible; "
                  f"the deposited 'TNBC' row duplicated the marginal row and is dropped")
            rows.append(dict(cohort=coh, stratum=None, evaluable=False,
                             reason="cohort is uniformly one IHC stratum; a within-subtype test "
                                    "is undefined and the previously reported row duplicated "
                                    "the marginal result"))
            continue
        for lab in ["HR+/HER2-", "HER2+", "TNBC"]:
            m = (st == lab).values
            n, ev_ = int(m.sum()), int(y[m].sum())
            rec = dict(cohort=coh, stratum=lab, n=n, events=ev_,
                       epv=round(ev_ / 5, 2) if ev_ else 0.0)
            if n < 20 or min(ev_, n - ev_) < MIN_EVENTS:
                rec.update(evaluable=False,
                           reason=f"fewer than {MIN_EVENTS} events in a class")
                rows.append(rec); continue
            R = sm.add_constant(Xr[m]); F = np.column_stack([R, s[m]])
            mr = sm.Logit(y[m], R).fit(disp=0); mf = sm.Logit(y[m], F).fit(disp=0)
            fit = R @ np.linalg.lstsq(R, s[m], rcond=None)[0]
            b, se = mf.params[-1], mf.bse[-1]
            rec.update(evaluable=True,
                       marginal=round(float(roc_auc_score(y[m], s[m])), 4),
                       residual=round(float(roc_auc_score(y[m], s[m] - fit)), 4),
                       or_sd=round(float(np.exp(b)), 3),
                       or_ci=[round(float(np.exp(b - 1.96 * se)), 3),
                              round(float(np.exp(b + 1.96 * se)), 3)],
                       lr_p=round(float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1)), 4),
                       mde_or_sd_80pct=mde_or(n, ev_))
            rows.append(rec)
            print(f"  {coh:9s} {lab:10s} n={n:3d} ev={ev_:3d} EPV={rec['epv']:.1f} "
                  f"resid={rec['residual']:.3f} p={rec['lr_p']:.4f} "
                  f"MDE OR/SD={rec['mde_or_sd_80pct']}")

    out = dict(spec="A4 as amended by addendum §13.11",
               note=("Per-stratum MDE is the smallest added OR/SD detectable at 80% power at "
                     "that stratum's n and event count. It is a power calculation, not a "
                     "measurement, and it is reported beside every non-significant cell so "
                     "that 'no cell reaches nominal significance' is not read as evidence of "
                     "absence."), cells=rows)
    (RES / "A4/A4_within_subtype_extended.json").write_text(json.dumps(out, indent=2))
    ev_cells = [r for r in rows if r.get("evaluable")]
    print(f"\nevaluable cells: {len(ev_cells)}; "
          f"reaching p<0.05: {[(r['cohort'], r['stratum']) for r in ev_cells if r['lr_p']<0.05] or 'none'}")
    print(f"cells with EPV <= 3: "
          f"{sum(1 for r in ev_cells if r['epv'] <= 3)} of {len(ev_cells)}")


if __name__ == "__main__":
    main()
