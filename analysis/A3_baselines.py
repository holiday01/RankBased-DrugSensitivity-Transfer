#!/usr/bin/env python3
"""A3 -- alternative and strengthened baselines, and the instrument positive controls.

Addendum section 6 (A3) as amended by section 13.2. Section 6's original A3 was a single
robustness check; section 13.2 splits it because a robustness check and a positive control are
different nested tests and section 6 conflated them.

  A3-r    robustness: does the SIGNATURE still add nothing against a stronger baseline?
          (i)   measured IHC ER/PR/HER2 in place of the transcriptomic axes
          (ii)  four axes + grade + cT + cN
          (iii) union of both
          Governed by section 10 condition 2.

  A3-pc1  instrument positive control, NON-transcriptomic:
          LRT[axes + grade + cT + cN] vs [axes].  Certifies the machinery on a real covariate,
          but in an easy regime -- clinical stage is near-orthogonal to a transcriptome-derived
          basis, whereas the signatures under test are 77-89% collinear with it. Reported with
          that caveat in the same sentence as the result.

  A3-pc2  instrument positive control, TRANSCRIPTOMIC and near-collinear:
          leave-one-axis-out. Reduced = 3 of the 4 axes; test whether the held-out axis is
          recovered by the identical nested LRT. This is the regime the real test operates in.

A3-pc1 and A3-pc2 sit OUTSIDE every BH family and are reported with nominal p (section 13.2).

Covariate coding is fixed by section 13.2 so the control is the same control in every cohort:
grade ordinal 1 df; cT ordinal 1 df; cN harmonised to BINARY (N0 vs N+) 1 df, because GSE32646
records nodal status as positive/negative only. A categorical version is a secondary arm.

Cohort availability, verified against the deposited phenotype tables:
  grade + cT + cN  : GSE16446, GSE25066, GSE32646.  GSE22226 deposits NO nodal field of any
                     kind, so A3-pc1 and A3-r(ii) cannot run there -- section 6 lists it, and
                     that is a section 9 deviation logged in decisions_addendum.md.
  measured IHC     : GSE25066, GSE32646 (ER/PR/HER2); GSE22226 (ER and HER2 only, no PR).
                     NOT GSE16446, which deposits no ER or PR IHC -- `esr1bimod` is an ESR1
                     mRNA bimodal call, not immunohistochemistry.
GSE32646 is the PRE-SPECIFIED PRIMARY CONTROL COHORT: 100% complete on all of grade, cT, cN and
IHC, and the only genuinely independent cohort in the study.

Section 13.3's four-branch interpretation rule is applied per block per cohort and is evaluated
here exactly as written, before the numbers were seen.

Analysis sets are the section 4.1 de-duplicated ones. Writes results/A3/A3_baselines.json.
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

B = 2000
SEED = 42
OUT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A3")
R2_REDUNDANT = 0.50   # section 13.3 branch 3 vs branch 4
PHEN = {
    "GSE16446": "external_data/GSE16446/GSE16446_phenotype.tsv",
    "GSE22226": "external_data/GSE22226/GSE22226_phenotype.tsv",
    "GSE25066": "external_data/GSE25066/GSE25066_phenotype.tsv",
    "GSE32646": "resubmission_v2/results/phase_2/GSE32646_phenotype.tsv",
}


def _col(df: pd.DataFrame, *frags: str):
    for f in frags:
        for c in df.columns:
            if f.lower() in c.lower():
                return c
    return None


def clinical(cohort: str, index) -> pd.DataFrame:
    """Grade / cT / cN on the fixed coding of section 13.2. NaN where not deposited."""
    d = pd.read_csv(D.ROOT / PHEN[cohort], sep="\t", dtype=str).set_index("patient_id")
    d = d.reindex(index)
    out = pd.DataFrame(index=index, dtype=float)

    g = d[_col(d, "histologic_grade", "histological_grade", "grade")]
    out["grade"] = pd.to_numeric(g.where(~g.astype(str).str.contains("4", na=False)),
                                 errors="coerce")

    if cohort == "GSE22226":                      # 1 = <=2cm, 2 = >2-5cm, 3 = >5cm
        t = pd.to_numeric(d[_col(d, "clinical_t_stage")], errors="coerce")
    else:
        t = d[_col(d, "clinical_t_stage")] if _col(d, "clinical_t_stage") else d["t"]
        t = pd.to_numeric(t.astype(str).str.extract(r"(\d)")[0], errors="coerce")
    out["cT"] = t

    ncol = _col(d, "clinical_nodal_status", "lymph_node_status")
    if ncol is None and "n" in d.columns:
        ncol = "n"
    if ncol is None:
        out["cN"] = np.nan                        # GSE22226: not deposited
    else:
        s = d[ncol].astype(str).str.upper()
        out["cN"] = np.where(s.str.contains("NEGATIVE") | s.str.match(r"N?0"), 0.0,
                             np.where(s.isin(["NAN", "NONE", ""]), np.nan, 1.0))
    return out


def ihc(cohort: str, index) -> pd.DataFrame:
    d = pd.read_csv(D.ROOT / PHEN[cohort], sep="\t", dtype=str).set_index("patient_id")
    d = d.reindex(index)
    out = pd.DataFrame(index=index, dtype=float)

    def bin_(col, pos=("P", "POSITIVE", "1")):
        if col is None:
            return pd.Series(np.nan, index=index)
        s = d[col].astype(str).str.upper().str.strip()
        v = pd.Series(np.nan, index=index)
        v[s.isin(pos)] = 1.0
        v[s.isin(["N", "NEGATIVE", "0"])] = 0.0
        return v

    if cohort == "GSE22226":
        out["ER"] = bin_(_col(d, "er_(0="))
        out["HER2"] = bin_(_col(d, "her2_(0="))
    elif cohort == "GSE16446":
        return pd.DataFrame(index=index)           # no ER/PR IHC deposited
    else:
        out["ER"] = bin_(_col(d, "er_status_ihc"))
        out["PR"] = bin_(_col(d, "pr_status_ihc"))
        out["HER2"] = bin_(_col(d, "her2_status_fish", "her2_status"))
    return out


def lrt(y, Xred, add) -> dict:
    """Nested LRT, df = number of added columns. Same routine as the paper's primary test."""
    add = np.asarray(add, dtype=float)
    if add.ndim == 1:
        add = add[:, None]
    R = sm.add_constant(Xred) if Xred is not None and Xred.shape[1] else np.ones((len(y), 1))
    F = np.column_stack([R, add])
    mr, mf = sm.Logit(y, R).fit(disp=0), sm.Logit(y, F).fit(disp=0)
    chi2, df = 2 * (mf.llf - mr.llf), add.shape[1]
    r = dict(lr_chi2=float(chi2), df=int(df), lr_p=float(stats.chi2.sf(chi2, df)),
             auc_reduced=float(roc_auc_score(y, mr.predict(R))),
             auc_full=float(roc_auc_score(y, mf.predict(F))),
             n=int(len(y)), events=int(y.sum()),
             epv=float(y.sum() / (F.shape[1] - 1)),
             converged=bool(mf.mle_retvals["converged"]))
    if df == 1:
        b, se = mf.params[-1], mf.bse[-1]
        r.update(added_or_sd=float(np.exp(b)),
                 or_ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))])
    return r


def unadjusted(y, x) -> dict:
    """Section 13.3 step 1: the planted effect size, always reported."""
    X = sm.add_constant(np.asarray(x, dtype=float)[:, None])
    m = sm.Logit(y, X).fit(disp=0)
    b, se = m.params[1], m.bse[1]
    rng = np.random.default_rng(SEED)
    a = [roc_auc_score(y[i], x[i]) for i in (rng.integers(0, len(y), len(y)) for _ in range(500))
         if len(np.unique(y[i])) > 1] if False else None
    boot = []
    for _ in range(500):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) > 1:
            boot.append(roc_auc_score(y[i], np.asarray(x)[i]))
    auc = float(roc_auc_score(y, x))
    return dict(or_per_unit=float(np.exp(b)),
                or_ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                wald_p=float(m.pvalues[1]), auroc=auc,
                auroc_ci=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))])


def r2_on_axes(block, Xax) -> float:
    A = sm.add_constant(Xax)
    b = np.asarray(block, dtype=float)
    fit = A @ np.linalg.lstsq(A, b, rcond=None)[0]
    ss = ((b - b.mean()) ** 2).sum()
    return float(1 - ((b - fit) ** 2).sum() / ss) if ss > 0 else np.nan


def branch(unadj: dict, nested: dict, r2: float) -> str:
    """Section 13.3 step 2, first match wins. Written before the numbers were seen."""
    null_assoc = (unadj["or_ci"][0] <= 1.0 <= unadj["or_ci"][1]) and unadj["auroc_ci"][0] <= 0.5
    if null_assoc:
        return "1 CONTROL UNINFORMATIVE"
    if nested["lr_p"] < 0.05:
        return "2 INSTRUMENT VALIDATED"
    if not np.isnan(r2) and r2 >= R2_REDUNDANT:
        return "3 REDUNDANCY CONFIRMED"
    return "4 INSTRUMENT UNDERPOWERED"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, b_ov, _, _, _ = D.overlap()
    src = _cell_builders()
    res = {"spec": "ANALYSIS_ADDENDUM_2026-08.md section 6 (A3) as amended by 13.2/13.3",
           "primary_control_cohort": "GSE32646", "cohorts": {}}

    for coh in ["GSE32646", "GSE25066", "GSE16446", "GSE22226"]:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        y = ax["pCR"].values.astype(int)
        Xax = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
        sig, _ = src["sig_of"]("Framework", expr)
        sig = D.zscore(sig).values
        cl, ic = clinical(coh, ax.index), ihc(coh, ax.index)
        blocks = {c: cl[c].values for c in cl.columns if cl[c].notna().sum() > 0}
        entry = {"n": int(len(y)), "events": int(y.sum()),
                 "fields": {"clinical": list(blocks), "ihc": list(ic.columns)},
                 "A3_pc1_blocks": {}, "A3_pc2_axes": {}, "A3_r": {}}
        print(f"\n=== {coh}  n={len(y)} events={int(y.sum())} "
              f"clinical={list(blocks)} ihc={list(ic.columns)} ===")

        # ---- A3-pc1: each clinical block, and the joint block, against the four axes ----
        for name in list(blocks) + (["joint"] if len(blocks) > 1 else []):
            cols = list(blocks) if name == "joint" else [name]
            sub = cl[cols].dropna()
            idx = ax.index.get_indexer(sub.index)
            if len(sub) < 30 or y[idx].sum() < 5:
                continue
            yv, Xv, Bv = y[idx], Xax[idx], sub.values
            n_ = lrt(yv, Xv, Bv)
            if name == "joint":
                u = unadjusted(yv, Bv @ np.linalg.lstsq(sm.add_constant(Bv), yv.astype(float),
                                                        rcond=None)[0][1:])
                r2 = float(np.mean([r2_on_axes(Bv[:, j], Xv) for j in range(Bv.shape[1])]))
            else:
                u, r2 = unadjusted(yv, Bv[:, 0]), r2_on_axes(Bv[:, 0], Xv)
            v = branch(u, n_, r2)
            entry["A3_pc1_blocks"][name] = dict(unadjusted=u, nested=n_, r2_on_axes=r2,
                                                branch=v, n_used=int(len(sub)))
            print(f"  pc1 {name:7s} unadj OR={u['or_per_unit']:.2f}"
                  f"[{u['or_ci'][0]:.2f},{u['or_ci'][1]:.2f}] AUROC={u['auroc']:.3f} | "
                  f"LRT p={n_['lr_p']:.4f} {n_['auc_reduced']:.3f}->{n_['auc_full']:.3f} "
                  f"R2ax={r2:.3f} | {v}")

        # ---- A3-pc2: leave-one-axis-out (transcriptomic, near-collinear) ----
        for j, a in enumerate(D.AXES):
            others = np.delete(Xax, j, axis=1)
            n_ = lrt(y, others, Xax[:, j])
            u = unadjusted(y, Xax[:, j])
            r2 = r2_on_axes(Xax[:, j], others)
            v = branch(u, n_, r2)
            entry["A3_pc2_axes"][a] = dict(unadjusted=u, nested=n_, r2_on_others=r2, branch=v)
            print(f"  pc2 {a:12s} unadj AUROC={u['auroc']:.3f} | LRT p={n_['lr_p']:.4f} "
                  f"{n_['auc_reduced']:.3f}->{n_['auc_full']:.3f} R2oth={r2:.3f} | {v}")

        # ---- A3-r: the signature against each strengthened baseline ----
        variants = {"four axes (reference)": Xax}
        if len(ic.columns):
            sub = ic.dropna()
            variants["(i) measured IHC only"] = ("idx", sub)
        if {"grade", "cT", "cN"} <= set(blocks):
            sub = cl[["grade", "cT", "cN"]].dropna()
            variants["(ii) axes + grade + cT + cN"] = ("idx", sub, "axes")
            if len(ic.columns):
                both = pd.concat([cl[["grade", "cT", "cN"]], ic], axis=1).dropna()
                variants["(iii) union"] = ("idx", both, "axes")
        for label, spec in variants.items():
            if isinstance(spec, np.ndarray):
                idx, Xred = np.arange(len(y)), spec
            else:
                sub = spec[1]
                idx = ax.index.get_indexer(sub.index)
                Xred = (np.column_stack([Xax[idx], sub.values]) if len(spec) > 2
                        else sub.values.astype(float))
            if len(idx) < 30 or y[idx].sum() < 5:
                continue
            n_ = lrt(y[idx], Xred, sig[idx])
            entry["A3_r"][label] = n_
            print(f"  r   {label:28s} n={n_['n']:4d} base AUROC={n_['auc_reduced']:.3f} "
                  f"| +sig LRT p={n_['lr_p']:.4f} OR/SD={n_.get('added_or_sd', float('nan')):.3f} "
                  f"{'FIRES' if n_['lr_p'] < 0.05 else 'null'}")
        res["cohorts"][coh] = entry

    (OUT / "A3_baselines.json").write_text(json.dumps(res, indent=2))
    print(f"\nwrote {OUT / 'A3_baselines.json'}")

    print("\n" + "=" * 78)
    print("SECTION 13.3 VERDICT SUMMARY (positive controls)")
    for coh, e in res["cohorts"].items():
        for nm, d_ in e["A3_pc1_blocks"].items():
            print(f"  pc1  {coh:9s} {nm:7s} -> {d_['branch']}")
    for coh, e in res["cohorts"].items():
        fired = [a for a, d_ in e["A3_pc2_axes"].items() if d_["nested"]["lr_p"] < 0.05]
        print(f"  pc2  {coh:9s} axes recovered by the LRT: {fired or 'none'}")


if __name__ == "__main__":
    main()
