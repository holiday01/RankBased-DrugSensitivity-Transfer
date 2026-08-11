#!/usr/bin/env python3
"""A0i -- patient-level cluster bootstrap (addendum section 5.2).

This is the PRIMARY estimator for any equivalence claim, and it supplies the per-cell
confidence intervals that section 6's verdict classes 2, 3 and 4 turn on. It must run before
A2 and A13, because every interval those analyses report is defined as a percentile interval
of these replicates.

Section 5.2 fixes the parts that would otherwise be analyst choices:
  * resampling unit  = the distinct patient of section 4.1 (de-duplicated)
  * B                = 2000, numpy default_rng seed 42
  * draw size        = |union| patients, the de-duplicated union itself
  * multiplicity     = a patient drawn m times contributes m rows to each of their cells
  * weights          = each cell's observed, un-resampled SE, held fixed across replicates
  * cross-cohort     = a patient in two cohorts is drawn once and enters every cell holding them
  * <5 events either class in a replicate -> that cell drops for that replicate, the replicate
    is kept, and the realised k is recorded
  * pooled SE        = SD of the 2000 replicate pools;  TOST p = Phi((pooled - margin)/SE_boot)

One reading is not pinned by section 5.2 and is recorded here rather than left implicit:
**scores and axes are computed once on the de-duplicated analysis set and then restricted to
the drawn patients.** Section 5.2 says each cell is "restricted to the drawn patients", not
rescored on them. Recomputing the five per-gene z-scored signatures inside each replicate
would fold scorer instability into an interval meant to carry sampling variability, and would
make the interval depend on who else was drawn. Documented as an interpretation, not a
deviation.

Also corrects the rho = 1 column, which section 5.2 flags as wrong in A0e_dedup.py:
  Var = sum_c (sum_{i in c} 1/se_i)^2 / (sum_i w_i)^2,  w_i = 1/se_i^2
not SE_naive * sqrt(k/C), which holds only under equal cell weights.

Writes results/A0/A0i_cluster_bootstrap.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import A0e_dedup as D  # noqa: E402  -- reuse the deposited cell definitions verbatim

B = 2000
SEED = 42
MARGINS = [0.53, 0.55, 0.575, 0.60]
MIN_EVENTS = 5
OUT = D.ROOT / "revise_bioadv_386/results/A0" if False else Path(
    "/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A0")


def residual_point(y: np.ndarray, A: np.ndarray, s: np.ndarray, idx: np.ndarray) -> float:
    """Residual AUROC on a row subset. Point estimate only -- the outer cluster bootstrap
    supplies the uncertainty, so the inner bootstrap of A0e.residual_auroc is not run."""
    yi, Ai, si = y[idx], A[idx], s[idx]
    if yi.sum() < MIN_EVENTS or (len(yi) - yi.sum()) < MIN_EVENTS:
        return np.nan
    beta = np.linalg.lstsq(Ai, si, rcond=None)[0]
    return roc_auc_score(yi, si - Ai @ beta)


def fe_pool_rho1(points: np.ndarray, ses: np.ndarray, cohorts: list[str]) -> dict:
    """Naive inverse-variance pool plus the exact rho = 1 block-worst-case SE."""
    w = 1.0 / ses ** 2
    pooled = float((w * points).sum() / w.sum())
    se_naive = float(np.sqrt(1.0 / w.sum()))
    var = 0.0
    for c in sorted(set(cohorts)):
        m = np.array([x == c for x in cohorts])
        var += (1.0 / ses[m]).sum() ** 2
    se_rho1 = float(np.sqrt(var) / w.sum())
    return dict(pooled=pooled, se_naive=se_naive, se_rho1=se_rho1,
                rho1_over_naive=se_rho1 / se_naive)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a_ov, b_ov, shared, _, _ = D.overlap()
    print(f"overlap: {len(shared)} shared trial identifiers "
          f"({len(a_ov)} GSE22226 ids, {len(b_ov)} GSE25066 ids)")

    # ---- rebuild the PRIMARY arm's cells, keeping patient identity this time ----
    # A0e.main() builds these inline; re-derive here so each cell carries its patient index.
    src = _cell_builders()
    cells = []
    for name, coh, stratum in D_CELLS:
        expr, ax = src["load"](coh)
        if stratum == "HR+":
            k = ax.index.intersection(src["hrpos"])
            expr, ax = expr.loc[k], ax.loc[k]
        # PRIMARY arm: the 65 shared patients are removed from GSE25066, never from GSE22226
        keep = ~ax.index.isin(b_ov)
        expr, ax = expr.loc[keep], ax.loc[keep]
        s, ng = src["sig_of"](name, expr)
        y = ax["pCR"].values.astype(int)
        Xr = np.column_stack([D.zscore(ax[a]).values for a in D.AXES])
        A = sm.add_constant(Xr)
        sv = D.zscore(s).values
        pt, se = D.residual_auroc(y, Xr, s)
        cells.append(dict(name=name, cohort=coh, stratum=stratum or "all",
                          n=int(len(y)), events=int(y.sum()), genes=int(ng),
                          point=float(pt), se=float(se),
                          _y=y, _A=A, _s=sv, _pid=[f"{coh}:{i}" for i in ax.index]))
        print(f"  {name:11s} {coh:9s} n={len(y):4d} ev={int(y.sum()):3d} "
              f"resid={pt:.4f} se={se:.4f}")

    # ---- the resampling universe: distinct patients across all cells ----
    union = sorted({p for c in cells for p in c["_pid"]})
    uidx = {p: i for i, p in enumerate(union)}
    n_union = len(union)
    for c in cells:
        c["_u"] = np.array([uidx[p] for p in c["_pid"]], dtype=np.int64)
    print(f"union of de-duplicated analysis sets: {n_union} distinct patients")

    points = np.array([c["point"] for c in cells])
    ses = np.array([c["se"] for c in cells])
    cohorts = [c["cohort"] for c in cells]
    obs = fe_pool_rho1(points, ses, cohorts)
    w = 1.0 / ses ** 2

    # ---- B replicates ----
    rng = np.random.default_rng(SEED)
    rep_pool = np.full(B, np.nan)
    rep_cell = np.full((B, len(cells)), np.nan)
    rep_k = np.zeros(B, dtype=int)
    for b in range(B):
        draw = rng.integers(0, n_union, n_union)
        mult = np.bincount(draw, minlength=n_union)
        vals, wts = [], []
        for j, c in enumerate(cells):
            m = mult[c["_u"]]
            if m.sum() == 0:
                continue
            idx = np.repeat(np.arange(len(m)), m)
            v = residual_point(c["_y"], c["_A"], c["_s"], idx)
            if np.isnan(v):
                continue  # <5 events either class: cell drops, replicate kept
            rep_cell[b, j] = v
            vals.append(v)
            wts.append(w[j])
        if vals:
            vals, wts = np.array(vals), np.array(wts)
            rep_pool[b] = (wts * vals).sum() / wts.sum()
            rep_k[b] = len(vals)
        if (b + 1) % 250 == 0:
            print(f"  replicate {b + 1}/{B}")

    ok = ~np.isnan(rep_pool)
    se_boot = float(np.nanstd(rep_pool[ok], ddof=1))
    pooled = obs["pooled"]

    result = dict(
        spec="ANALYSIS_ADDENDUM_2026-08.md section 5.2",
        arm="PRIMARY: remove from GSE25066",
        B=B, seed=SEED, n_union=n_union,
        replicates_usable=int(ok.sum()),
        realised_k=dict(mean=float(rep_k[ok].mean()), min=int(rep_k[ok].min()),
                        max=int(rep_k[ok].max()),
                        full_k_fraction=float((rep_k[ok] == len(cells)).mean())),
        observed=obs,
        se_boot=se_boot,
        boot_over_naive=se_boot / obs["se_naive"],
        boot_over_rho1=se_boot / obs["se_rho1"],
        ci_boot_percentile=[float(np.nanpercentile(rep_pool[ok], 2.5)),
                            float(np.nanpercentile(rep_pool[ok], 97.5))],
        tost_boot={f"{m}": float(stats.norm.cdf((pooled - m) / se_boot)) for m in MARGINS},
        tost_naive={f"{m}": float(stats.norm.cdf((pooled - m) / obs["se_naive"]))
                    for m in MARGINS},
        tost_rho1={f"{m}": float(stats.norm.cdf((pooled - m) / obs["se_rho1"]))
                   for m in MARGINS},
        cells=[{k: v for k, v in c.items() if not k.startswith("_")} for c in cells],
    )
    # per-cell percentile CIs -- section 5.2: "the same 2000 replicates supply the per-cell CIs"
    for j, c in enumerate(result["cells"]):
        col = rep_cell[:, j]
        col = col[~np.isnan(col)]
        c["ci_cluster_percentile"] = [float(np.percentile(col, 2.5)),
                                      float(np.percentile(col, 97.5))]
        c["se_cluster"] = float(col.std(ddof=1))
        c["replicates_evaluable"] = int(len(col))
        c["tost_cell_0.55"] = float(stats.norm.cdf((c["point"] - 0.55) / c["se_cluster"]))

    (OUT / "A0i_cluster_bootstrap.json").write_text(json.dumps(result, indent=2))

    # ---- report ----
    print(f"\npooled (PRIMARY arm)            {pooled:.4f}")
    print(f"SE naive / rho=1 / bootstrap    {obs['se_naive']:.5f} / "
          f"{obs['se_rho1']:.5f} / {se_boot:.5f}")
    print(f"bootstrap / naive = {result['boot_over_naive']:.2f}x ; "
          f"bootstrap / rho=1 = {result['boot_over_rho1']:.2f}x")
    print(f"realised k: mean {result['realised_k']['mean']:.2f}, "
          f"all-{len(cells)} in {result['realised_k']['full_k_fraction']*100:.1f}% of replicates")
    print("\nTOST one-sided, by margin:")
    print(f"  {'margin':>7s} {'naive':>10s} {'rho=1':>10s} {'BOOTSTRAP':>12s}   verdict(a=0.05)")
    for m in MARGINS:
        pb = result["tost_boot"][f"{m}"]
        print(f"  {m:>7.3f} {result['tost_naive'][f'{m}']:10.2e} "
              f"{result['tost_rho1'][f'{m}']:10.2e} {pb:12.2e}   "
              f"{'PASS' if pb < 0.05 else 'FAIL'}")
    print("\nper-cell residual AUROC with cluster-bootstrap percentile CI:")
    for c in result["cells"]:
        lo, hi = c["ci_cluster_percentile"]
        print(f"  {c['name']:11s} {c['cohort']:9s} {c['point']:.4f} [{lo:.4f}, {hi:.4f}] "
              f"n={c['n']:4d} ev={c['events']:3d}  TOST@0.55 p={c['tost_cell_0.55']:.3f} "
              f"{'PASS' if c['tost_cell_0.55'] < 0.05 else 'FAIL'}")
    print(f"\nwrote {OUT / 'A0i_cluster_bootstrap.json'}")


D_CELLS = [("Framework", "GSE32646", None), ("DoxVariant", "GSE25066", None),
           ("GGI", "GSE25066", None), ("RSlike", "GSE25066", "HR+"),
           ("DLDA30", "GSE16446", None), ("TIS", "GSE25066", None),
           ("cytolytic", "GSE22226", None), ("HRD", "GSE22226", None),
           ("TGFBstrom", "GSE25066", None)]


def _cell_builders() -> dict:
    """Rebuild A0e's gene lists and scorer dispatch, so cells are byte-identical to the
    deposited definitions and only the patient bookkeeping differs."""
    hrd = D.script_lists(D.SCR / "phase_8/HRD_gse22226.py",
                         "HRD_GENES_SEVERSON")["HRD_GENES_SEVERSON"]
    tis = D.script_lists(D.SCR / "phase_7/TIS_gse25066.py", "TIS_GENES")["TIS_GENES"]
    tgfb = D.script_lists(D.SCR / "phase_8/TGFBstrom_gse25066.py",
                          "TGFB_STROMAL_GENES")["TGFB_STROMAL_GENES"]
    dl = D.script_lists(D.SCR / "phase_7/DLDA30_gse16446.py", "DLDA30_UP", "DLDA30_DOWN")
    ggi = D.script_lists(D.SCR / "phase_6/gge_GGI_gse25066.py")
    rs = D.script_lists(D.SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")
    fw = pd.read_csv(D.ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv",
                     sep="\t")
    dv = pd.read_csv(D.ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    ggi_up = [k for k in ggi if "UP" in k.upper() and "DOWN" not in k.upper()]
    ggi_dn = [k for k in ggi if "DOWN" in k.upper()]
    hrpos = set(pd.read_csv(D.ROOT / "external_data/GSE25066/GSE25066_pcr_HR+.tsv",
                            sep="\t")["patient_id"])

    def sig_of(name, expr):
        if name == "Framework":
            return D.singscore(expr, fw.loc[fw.direction == "sensitivity", "gene_symbol"],
                               fw.loc[fw.direction == "resistance", "gene_symbol"])
        if name == "DoxVariant":
            return D.singscore(expr, dv.loc[dv.direction == "sensitivity", "gene_symbol"],
                               dv.loc[dv.direction == "resistance", "gene_symbol"])
        if name == "GGI":
            return D.singscore(expr, ggi[ggi_up[0]], ggi[ggi_dn[0]])
        if name == "RSlike":
            return D.singscore(expr, rs["RS_UP"], rs["RS_DOWN"])
        if name == "DLDA30":
            u, _ = D.mean_z(expr, dl["DLDA30_UP"])
            d, _ = D.mean_z(expr, dl["DLDA30_DOWN"])
            return u - d, len(dl["DLDA30_UP"]) + len(dl["DLDA30_DOWN"])
        return D.mean_z(expr, {"TIS": tis, "cytolytic": ["GZMA", "PRF1"],
                               "HRD": hrd, "TGFBstrom": tgfb}[name])

    return dict(load=D.load, sig_of=sig_of, hrpos=hrpos)


if __name__ == "__main__":
    main()
