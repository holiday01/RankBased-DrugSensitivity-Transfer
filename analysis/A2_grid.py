#!/usr/bin/env python3
"""A2 -- the full signature x cohort added-value grid.

Addendum section 6 (A2) plus section 13.9. This is the direct answer to Reviewer 2's major
point M1, whose either/or (narrow the title, or produce the grid) the submitted revision
answered with neither.

Family A is a FIXED list of 58 cells, enumerated in addendum section 4.2:
  * 10 signatures x 4 pCR cohorts                                     -> 40
  * 6 CCLE drugs x 3 discovery cohorts (CARBOPLATIN, CYCLOPHOSPHAMIDE,
    FLUOROURACIL, PACLITAXEL, ETOPOSIDE, TOP2_POISON_CONSENSUS)       -> 18
m = 58 is a fixed denominator; a NOT-EVALUABLE cell enters BH at p = 1 and is reported by name
with the clause that fired. m is never relaxed to make a cell reachable.

Outside Family A, per section 13.9: DASATINIB, ERLOTINIB and GEFITINIB on the same four
cohorts as a labelled wrong-mechanism NEGATIVE CONTROL arm, reported with nominal p. Excluding
them left the grid with no internal calibration -- a reader could not tell whether "nothing
beats the axes" means the axes dominate or the test is inert.

Verdict classes follow the six-class rule fixed in advance, evaluated in order, first match
wins: NOT-EVALUABLE / ADDS-VALUE / INVERSE-DIRECTION / REDUNDANT-WITH-AXES /
NO-MARGINAL-SIGNAL / INDETERMINATE. All CIs are the patient-level cluster bootstrap percentile
intervals of section 5.2 (B = 2000, seed 42) -- NOT the naive bootstrap used in A0b-A0e --
which is why A0i had to run first.

Per section 13.9 each cell also records scorer type, cohort-dependence, nominal and resolved
set sizes per direction, and deployment fidelity, because five of the nine published signatures
are scored by cohort-wide per-gene z-scoring and are not cohort-free in the sense the Methods
claim. The grid must not mix the two silently.

Analysis sets are the section 4.1 de-duplicated ones. Writes results/A2/A2_grid.json.
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

B, SEED = 2000, 42
OUT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A2")
COHORTS = ["GSE16446", "GSE22226", "GSE25066", "GSE32646"]
DISCOVERY = ["GSE16446", "GSE22226", "GSE25066"]
CCLE6 = ["CARBOPLATIN", "CYCLOPHOSPHAMIDE", "FLUOROURACIL", "PACLITAXEL", "ETOPOSIDE",
         "TOP2_POISON_CONSENSUS"]
WRONG_MOA = ["DASATINIB", "ERLOTINIB", "GEFITINIB"]
FROZEN = "frozen_signatures/solid_only/quantile_30_oncokb_all"
MIN_EVENTS, MIN_RESOLVED = 5, 0.60


def frozen(drug: str):
    p = D.ROOT / FROZEN / f"{drug}.tsv"
    if not p.exists():
        return None
    t = pd.read_csv(p, sep="\t")
    return (list(t.loc[t.direction == "sensitivity", "gene_symbol"]),
            list(t.loc[t.direction == "resistance", "gene_symbol"]))


def mammaprint_lists():
    """MAMMAPRINT is a {gene: +/-1} dict in the deployed script. Read it by AST rather than
    exec, so the deployed module's import-time side effects do not run."""
    import ast
    tree = ast.parse((D.SCR / "phase_6/mammaprint_metabric.py").read_text())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and node.targets
                and getattr(node.targets[0], "id", None) == "MAMMAPRINT"):
            mp = ast.literal_eval(node.value)
            return ([g for g, w in mp.items() if w == +1],
                    [g for g, w in mp.items() if w == -1])
    raise RuntimeError("MAMMAPRINT dict not found")


def build_defs(src) -> dict:
    """name -> (up, down, scorer, cohort_dependent, fidelity, restrict)"""
    hrd = D.script_lists(D.SCR / "phase_8/HRD_gse22226.py",
                         "HRD_GENES_SEVERSON")["HRD_GENES_SEVERSON"]
    tis = D.script_lists(D.SCR / "phase_7/TIS_gse25066.py", "TIS_GENES")["TIS_GENES"]
    tgfb = D.script_lists(D.SCR / "phase_8/TGFBstrom_gse25066.py",
                          "TGFB_STROMAL_GENES")["TGFB_STROMAL_GENES"]
    dl = D.script_lists(D.SCR / "phase_7/DLDA30_gse16446.py", "DLDA30_UP", "DLDA30_DOWN")
    ggi = D.script_lists(D.SCR / "phase_6/gge_GGI_gse25066.py")
    rs = D.script_lists(D.SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")
    gu = [k for k in ggi if "UP" in k.upper() and "DOWN" not in k.upper()][0]
    gd = [k for k in ggi if "DOWN" in k.upper()][0]
    mp_up, mp_dn = mammaprint_lists()
    fw, dv = frozen("DOXORUBICIN"), None
    t = pd.read_csv(D.ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    dv = (list(t.loc[t.direction == "sensitivity", "gene_symbol"]),
          list(t.loc[t.direction == "resistance", "gene_symbol"]))

    S = "rank-singscore-bidirectional"
    Z = "per-gene-z-mean"
    d = {
        "Framework":  (*fw, S, False, "as-published", None),
        "DoxVariant": (*dv, S, False, "as-published", None),
        "GGI":        (ggi[gu], ggi[gd], S, False,
                       "approximation: 41 of 97 genes, rank-scored not scaled sum", None),
        "RSlike":     (rs["RS_UP"], rs["RS_DOWN"], S, False,
                       "approximation: 16-gene stand-in for the 21-gene RS", "HR+"),
        "MammaPrint": (mp_up, mp_dn, S, False,
                       "approximation: 70-gene list rank-scored, NOT the fixed-centroid assay",
                       None),
        "DLDA30":     (dl["DLDA30_UP"], dl["DLDA30_DOWN"], Z, True,
                       "approximation: mean-z difference, not Hess's diagonal LDA", None),
        "TIS":        (tis, [], Z, True, "approximation: unweighted mean z, no housekeeping "
                                         "normalisation, array not NanoString", None),
        "cytolytic":  (["GZMA", "PRF1"], [], Z, True,
                       "approximation: arithmetic not geometric mean", None),
        "HRD":        (hrd, [], Z, True,
                       "pathway-expression score, not a genomic-scar HRD measure", None),
        "TGFBstrom":  (tgfb, [], Z, True, "as-published", None),
    }
    for drug in CCLE6 + WRONG_MOA:
        f = frozen(drug)
        if f:
            d[drug] = (*f, S, False, "as-published", None)
    return d


def score(expr, up, down, scorer):
    if scorer.startswith("rank-singscore"):
        s, _ = D.singscore(expr, up, down)
        nu, nd = len(D.resolve(up, expr.columns)), len(D.resolve(down, expr.columns))
        return s, nu, nd
    if down:
        u, nu = D.mean_z(expr, up)
        dn, nd = D.mean_z(expr, down)
        return u - dn, nu, nd
    s, nu = D.mean_z(expr, up)
    return s, nu, 0


def resid(y, A, s, idx):
    yi, Ai, si = y[idx], A[idx], s[idx]
    if yi.sum() < MIN_EVENTS or (len(yi) - yi.sum()) < MIN_EVENTS:
        return np.nan
    return roc_auc_score(yi, si - Ai @ np.linalg.lstsq(Ai, si, rcond=None)[0])


def cell_stats(y, Xr, s):
    A = sm.add_constant(Xr)
    R, F = A, np.column_stack([A, s])
    mr, mf = sm.Logit(y, R).fit(disp=0), sm.Logit(y, F).fit(disp=0)
    chi2 = 2 * (mf.llf - mr.llf)
    b, se = mf.params[-1], mf.bse[-1]
    fit = A @ np.linalg.lstsq(A, s, rcond=None)[0]
    r2 = float(1 - ((s - fit) ** 2).sum() / ((s - s.mean()) ** 2).sum())
    return dict(marginal=float(roc_auc_score(y, s)),
                residual=float(resid(y, A, s, np.arange(len(y)))),
                r2_on_axes=r2, lr_chi2=float(chi2), lr_p=float(stats.chi2.sf(chi2, 1)),
                added_or_sd=float(np.exp(b)),
                or_ci=[float(np.exp(b - 1.96 * se)), float(np.exp(b + 1.96 * se))],
                dbic=float(np.log(len(y)) - chi2),
                auc_reduced=float(roc_auc_score(y, mr.predict(R))),
                auc_full=float(roc_auc_score(y, mf.predict(F))),
                converged=bool(mf.mle_retvals["converged"]))


def verdict(c) -> str:
    if c["not_evaluable"]:
        return "NOT-EVALUABLE"
    lo, hi = c["residual_ci"]
    if c.get("q") is not None and c["q"] < 0.10 and lo > 0.5:
        return "ADDS-VALUE"
    if hi < 0.5:
        return "INVERSE-DIRECTION"
    if (c["marginal"] >= 0.60 and c["residual"] <= 0.55 and lo <= 0.5 <= hi
            and c["r2_on_axes"] >= 0.50):
        return "REDUNDANT-WITH-AXES"
    if c["marginal"] < 0.60:
        return "NO-MARGINAL-SIGNAL"
    return "INDETERMINATE"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, b_ov, _, _, _ = D.overlap()
    src = _cell_builders()
    defs = build_defs(src)
    hrpos = src["hrpos"]

    data = {}
    for coh in COHORTS:
        expr, ax = src["load"](coh)
        if coh == "GSE25066":
            k = ~ax.index.isin(b_ov)
            expr, ax = expr.loc[k], ax.loc[k]
        data[coh] = dict(expr=expr, ax=ax, y=ax["pCR"].values.astype(int),
                         A=sm.add_constant(np.column_stack(
                             [D.zscore(ax[a]).values for a in D.AXES])),
                         pid=[f"{coh}:{i}" for i in ax.index])
        print(f"{coh}: n={len(ax)} events={int(ax['pCR'].sum())}")

    union = sorted({p for c in data.values() for p in c["pid"]})
    uidx = {p: i for i, p in enumerate(union)}
    for c in data.values():
        c["u"] = np.array([uidx[p] for p in c["pid"]], dtype=np.int64)
    print(f"union: {len(union)} distinct patients")

    plan = ([(s, c, "A") for s in list(defs)[:10] for c in COHORTS]
            + [(d, c, "A") for d in CCLE6 for c in DISCOVERY]
            + [(d, c, "NEG") for d in WRONG_MOA for c in COHORTS])

    cells = []
    for name, coh, fam in plan:
        if name not in defs:
            continue
        up, dn, scorer, cohdep, fid, restrict = defs[name]
        d = data[coh]
        expr, ax = d["expr"], d["ax"]
        sel = np.ones(len(ax), dtype=bool)
        vacancy = None
        if restrict == "HR+":
            if coh == "GSE25066":
                sel = ax.index.isin(hrpos)
            else:
                er = _er_positive(coh, ax.index)
                sel = er.values if er is not None else np.zeros(len(ax), dtype=bool)
            if sel.sum() < 30:
                vacancy = ("stratum-definition vacancy: the score is defined for HR-positive "
                           "disease and this cohort has too few HR-positive patients")
        rec = dict(signature=name, cohort=coh, family=fam, scorer=scorer,
                   cohort_dependent=cohdep, deployment_fidelity=fid,
                   n_up_nominal=len(up), n_down_nominal=len(dn), restrict=restrict)
        if vacancy:
            rec.update(n=int(sel.sum()), events=None, not_evaluable=True,
                       not_evaluable_clause=vacancy, lr_p=1.0,
                       residual_ci=[np.nan, np.nan], marginal=np.nan, residual=np.nan,
                       r2_on_axes=np.nan)
            cells.append(rec)
            continue
        e2, y2 = expr.loc[sel], d["y"][sel]
        A2m, u2 = d["A"][sel], d["u"][sel]
        s, nu, nd = score(e2, up, dn, scorer)
        s = D.zscore(s).values
        fu = nu / max(len(up), 1)
        fd = nd / max(len(dn), 1) if dn else None
        bad = (min(y2.sum(), len(y2) - y2.sum()) < MIN_EVENTS or fu < MIN_RESOLVED
               or (fd is not None and fd < MIN_RESOLVED))
        rec.update(n=int(len(y2)), events=int(y2.sum()), n_up_resolved=nu, n_down_resolved=nd,
                   frac_up_resolved=float(fu),
                   frac_down_resolved=(float(fd) if fd is not None else None),
                   not_evaluable=bool(bad))
        if bad:
            rec.update(not_evaluable_clause=("fewer than 5 events in a class" if
                                             min(y2.sum(), len(y2) - y2.sum()) < MIN_EVENTS
                                             else "gene resolution below 0.60"),
                       lr_p=1.0, residual_ci=[np.nan, np.nan], marginal=np.nan,
                       residual=np.nan, r2_on_axes=np.nan)
        else:
            rec.update(cell_stats(y2, A2m[:, 1:], s))
            rec["_y"], rec["_A"], rec["_s"], rec["_u"] = y2, A2m, s, u2
        cells.append(rec)
        if rec["not_evaluable"]:
            print(f"  {name:22s} {coh:9s} NOT-EVALUABLE ({rec['not_evaluable_clause']})")
        else:
            print(f"  {name:22s} {coh:9s} marg={rec['marginal']:.3f} "
                  f"resid={rec['residual']:.3f} R2ax={rec['r2_on_axes']:.3f} "
                  f"p={rec['lr_p']:.4f}")

    # ---- cluster bootstrap: one draw, every cell restricted to it (section 5.2) ----
    print(f"\ncluster bootstrap B={B} over {len(union)} patients ...")
    live = [c for c in cells if not c["not_evaluable"]]
    rep = np.full((B, len(live)), np.nan)
    rng = np.random.default_rng(SEED)
    for b in range(B):
        mult = np.bincount(rng.integers(0, len(union), len(union)), minlength=len(union))
        for j, c in enumerate(live):
            m = mult[c["_u"]]
            if m.sum() == 0:
                continue
            rep[b, j] = resid(c["_y"], c["_A"], c["_s"],
                              np.repeat(np.arange(len(m)), m))
        if (b + 1) % 500 == 0:
            print(f"  {b + 1}/{B}")
    for j, c in enumerate(live):
        col = rep[:, j][~np.isnan(rep[:, j])]
        c["residual_ci"] = [float(np.percentile(col, 2.5)), float(np.percentile(col, 97.5))]
        c["residual_se_cluster"] = float(col.std(ddof=1))
        c["bootstrap_replicates"] = int(len(col))

    # ---- BH across Family A, m = 58 fixed ----
    famA = [c for c in cells if c["family"] == "A"]
    m = len(famA)
    order = sorted(range(m), key=lambda i: famA[i]["lr_p"])
    q_prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        q_prev = min(q_prev, famA[i]["lr_p"] * m / k)
        famA[i]["q"] = float(min(q_prev, 1.0))
    for c in cells:
        if c["family"] != "A":
            c["q"] = None
        c["verdict"] = verdict(c)

    clean = [{k: v for k, v in c.items() if not k.startswith("_")} for c in cells]
    counts = {}
    for c in clean:
        if c["family"] == "A":
            counts[c["verdict"]] = counts.get(c["verdict"], 0) + 1
    res = dict(spec="ANALYSIS_ADDENDUM_2026-08.md section 6 (A2) + section 13.9",
               family_A_m=m, B=B, seed=SEED, n_union=len(union),
               verdict_counts=counts, cells=clean)
    (OUT / "A2_grid.json").write_text(json.dumps(res, indent=2, default=float))

    print(f"\nFamily A: m = {m}")
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k:24s} {v:3d}  ({v/m*100:.0f}%)")
    adds = [c for c in clean if c["verdict"] == "ADDS-VALUE"]
    inv = [c for c in clean if c["verdict"] == "INVERSE-DIRECTION"]
    print(f"\nADDS-VALUE cells (condition 1): {[(c['signature'], c['cohort']) for c in adds] or 'none'}")
    print(f"INVERSE-DIRECTION cells:        {[(c['signature'], c['cohort']) for c in inv] or 'none'}")
    print(f"REDUNDANT-WITH-AXES majority?   "
          f"{counts.get('REDUNDANT-WITH-AXES', 0)}/{m} -> "
          f"{'yes' if counts.get('REDUNDANT-WITH-AXES', 0) > m/2 else 'NO -- class counts enter the abstract'}")
    neg = [c for c in clean if c["family"] == "NEG" and not c["not_evaluable"]]
    print("\nwrong-mechanism negative controls (outside Family A, nominal p):")
    for c in sorted(neg, key=lambda x: x["lr_p"])[:6]:
        print(f"  {c['signature']:12s} {c['cohort']:9s} marg={c['marginal']:.3f} "
              f"resid={c['residual']:.3f} p={c['lr_p']:.4f}"
              f"{'  <-- fires' if c['lr_p'] < 0.05 else ''}")
    print(f"\nwrote {OUT / 'A2_grid.json'}")


def _er_positive(cohort: str, index):
    from A3_baselines import ihc
    t = ihc(cohort, index)
    return t["ER"] == 1.0 if "ER" in t.columns else None


if __name__ == "__main__":
    main()
