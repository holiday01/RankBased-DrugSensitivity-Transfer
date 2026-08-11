#!/usr/bin/env python3
"""A0e (v2) -- the GSE22226 / GSE25066 patient overlap, and which side to remove.

GSE22226's `i-spy_id` and GSE25066's `sample_id` are the same I-SPY 1 trial identifier.
65 patients are shared; every matched GSE25066 row carries source == 'ISPY'. The manuscript
calls these "three independent cohorts" (main.tex L111, L661) and pools them as independent.

REMOVAL DIRECTION IS NOT ARBITRARY. GSE25066 (Hatzis 2011) enrolled HER2-negative disease
only, so the overlap IS the HER2-negative half of GSE22226. Removing it from GSE22226 leaves
62 patients that are 62.5% HER2-positive against 31.7% in the full analysis set, retains all
11 trastuzumab-treated patients and removes none, and no longer matches the cohort the
manuscript describes at L156-159. Removing from GSE25066 instead takes 61 patients out of a
cohort that is HER2-negative throughout (485 negative / 6 positive), which leaves HER2 and
the pCR rate essentially unmoved (18.0% vs 20.6%, p = 0.74). It is NOT case-mix neutral in
every respect -- clinical T stage does shift (T3 50.8% vs 25.5%, p = 6.4e-04), a consequence
of I-SPY 1 requiring tumours >= 3 cm -- and it removes 61 of the 79 ISPY-source patients in
the analysis set rather than the whole source layer, because removal is by patient id.
Direction B is nonetheless primary, and direction A is reported as a sensitivity arm.

Both pools are additionally reported with the shared-patient structure honoured: after
de-duplication five cells still sit on the same GSE25066 patients (RSlike on a nested HR+
subset of them) and two on the same GSE22226 patients, so the naive fixed-effect SE is
understated. That is the circularity raised in review, and it is not fixed by
de-duplication.

Every cell is recomputed here from source, including the undeduplicated baseline, so the
before/after difference is attributable to de-duplication and nothing else. Gene lists are
read from the deployed scripts rather than retyped. Bidirectional signatures use the
CORRECTED singscore established in A0.
"""
import ast
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
SCR = ROOT / "resubmission_v2/scripts"
RES = ROOT / "resubmission_v2/results"
OUT = _OR / "A0"
AXES = ["P", "H", "B", "H_PAM50basal"]
MARGINS = [0.53, 0.55, 0.575, 0.60]
NBOOT, SEED = 1000, 42


def zscore(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def script_lists(path, *names):
    """Read module-level list constants without executing the script."""
    out = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and (not names or t.id in names):
                    try:
                        v = ast.literal_eval(node.value)
                        if isinstance(v, list) and v and all(isinstance(x, str) for x in v):
                            out[t.id] = v
                        elif isinstance(v, dict):
                            out[t.id] = v
                    except Exception:
                        pass
    return out


_ALIAS = None


def alias_map():
    """HGNC alias -> canonical, the same table the deployed scorers use."""
    global _ALIAS
    if _ALIAS is None:
        p = ROOT / "external_data/hgnc_alias_table.tsv"
        _ALIAS = {}
        if p.exists():
            df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
            cols = {c.lower(): c for c in df.columns}
            if "alias_symbol" in cols and "canonical_symbol" in cols:
                for _, r in df.iterrows():
                    a, c = r[cols["alias_symbol"]].upper(), r[cols["canonical_symbol"]].upper()
                    if a and c:
                        _ALIAS[a] = c
    return _ALIAS


def resolve(genes, columns):
    """Match by symbol, then by HGNC alias -- without the alias step, CTGF->CCN2 and three
    DoxVariant genes are silently dropped.

    Addendum section 13.1: an alias row X->Y is applied ONLY when X is not itself a canonical
    symbol somewhere in the table. 605 rows fail that test, and four of them collide with
    members of the frozen signature -- HGF->IL6, SMO->SMOX, CDH1->FZR1, MET->SLTM -- all four
    approved HGNC symbols mapped onto unrelated genes. Applying the table unguarded moves the
    headline result; the filter reproduces the exact-match resolution that Table 1 was built on.
    """
    up = {c.upper(): c for c in columns}
    am = alias_map()
    canonical = set(am.values())
    out = []
    for g in genes:
        gU = str(g).upper()
        if gU in up:
            out.append(up[gU])
        elif gU in am and am[gU] in up and gU not in canonical:
            out.append(up[am[gU]])
    return list(dict.fromkeys(out))


def mean_z(expr, genes):
    g = resolve(genes, expr.columns)
    sub = expr[g]
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1), len(g)


def singscore(expr, up, down):
    """Corrected bidirectional singscore: TotalScore = up + down_reversed (see A0)."""
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    U = resolve(up, expr.columns)
    D = resolve(down, expr.columns)

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1

    raw = norm(ranks[U].mean(axis=1), len(U)) + norm((1 - ranks[D]).mean(axis=1), len(D))
    return raw - raw.median(), len(U) + len(D)


def residual_auroc(y, Xr, s, nboot=NBOOT, seed=SEED):
    s = zscore(s).values
    A = sm.add_constant(Xr)

    def one(idx):
        Ai, si, yi = A[idx], s[idx], y[idx]
        if len(np.unique(yi)) < 2:
            return np.nan
        return roc_auc_score(yi, si - Ai @ np.linalg.lstsq(Ai, si, rcond=None)[0])

    pt = one(np.arange(len(y)))
    rng = np.random.default_rng(seed)
    bs = np.array([one(rng.integers(0, len(y), len(y))) for _ in range(nboot)])
    bs = bs[~np.isnan(bs)]
    return float(pt), float(bs.std(ddof=1))


def fe_pool(pts, ses, cohorts=None, label=""):
    """Inverse-variance fixed effect. Unrounded inputs only -- rounding before pooling
    moved the 4th decimal.

    With `cohorts`, a second SE is reported for the rho = 1 worst case: cells within a cohort
    perfectly correlated, cohorts independent. Five cells sit on the same GSE25066 patients
    and two on the same GSE22226 patients, so the naive SE is understated.

    An earlier version used SE_naive * sqrt(k/C) for this. That is the rho = 1 design effect
    only when the cells carry EQUAL weight; here 1/se ranges 13.7 to 33.0, and the shortcut
    lands 26% below the exact value. It also happened to coincide with the "clustered" column
    that version printed, so the two were the same number reported twice. Both defects are
    fixed here: the exact block form is used, and no figure is labelled "clustered" -- the
    patient-level cluster bootstrap of addendum section 5.2 has not been run, and nothing in
    this script estimates it.

    Exact form, with w_i = 1/se_i^2 and Cov(theta_i, theta_j) = se_i*se_j within a cohort:

        Var(pooled) = sum_c ( sum_{i in c} 1/se_i )^2 / (sum_i w_i)^2
    """
    p, se = np.asarray(pts, float), np.asarray(ses, float)
    w = 1 / se ** 2
    W = w.sum()
    pooled = float((w * p).sum() / W)
    sep = float(np.sqrt(1 / W))
    out = dict(label=label, k=len(p), pooled=pooled, se_naive=sep,
               ci=[pooled - 1.96 * sep, pooled + 1.96 * sep],
               **{f"tost_{m}": float(stats.norm.cdf((pooled - m) / sep)) for m in MARGINS})
    if cohorts is not None:
        co = np.asarray(cohorts)
        var = sum((np.sum(1 / se[co == c])) ** 2 for c in np.unique(co)) / W ** 2
        se_r1 = float(np.sqrt(var))
        out.update(se_rho1=se_r1, rho1_over_naive=se_r1 / sep,
                   ci_rho1=[pooled - 1.96 * se_r1, pooled + 1.96 * se_r1],
                   cluster_bootstrap="NOT RUN -- see addendum section 5.2 item 2",
                   **{f"tost_rho1_{m}": float(stats.norm.cdf((pooled - m) / se_r1))
                      for m in MARGINS})
    return out


def overlap():
    a = pd.read_csv(ROOT / "external_data/GSE22226/GSE22226_phenotype.tsv", sep="\t", dtype=str)
    b = pd.read_csv(ROOT / "external_data/GSE25066/GSE25066_phenotype.tsv", sep="\t", dtype=str)
    a["k"] = pd.to_numeric(a["i-spy_id"], errors="coerce")
    b["k"] = pd.to_numeric(b["sample_id"], errors="coerce")
    sh = set(a["k"].dropna()) & set(b["k"].dropna())
    return (set(a.loc[a["k"].isin(sh), "patient_id"]),
            set(b.loc[b["k"].isin(sh), "patient_id"]), sh, a, b)


def load(cohort):
    P = {"GSE16446": (ROOT / "external_data/GSE16446/GSE16446_expression.tsv",
                      RES / "phase_0/axes_GSE16446.tsv"),
         "GSE22226": (ROOT / "external_data/GSE22226/GSE22226_expression.tsv",
                      RES / "phase_0/axes_GSE22226.tsv"),
         "GSE25066": (ROOT / "external_data/GSE25066/GSE25066_expression.tsv",
                      RES / "phase_0/axes_GSE25066.tsv"),
         "GSE32646": (RES / "phase_2/GSE32646_expression.tsv", RES / "phase_2/axes_GSE32646.tsv")}
    ep, ap = P[cohort]
    e = pd.read_csv(ep, sep="\t", index_col=0)
    if e.shape[1] < 5000:
        e = e.T
    ax = pd.read_csv(ap, sep="\t", index_col=0).dropna(subset=AXES + ["pCR"])
    c = e.index.intersection(ax.index)
    return e.loc[c], ax.loc[c]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a_ov, b_ov, shared, phen22, phen25 = overlap()

    # ---------- gene lists, read from the deployed scripts ----------
    hrd = script_lists(SCR / "phase_8/HRD_gse22226.py", "HRD_GENES_SEVERSON")["HRD_GENES_SEVERSON"]
    tis = script_lists(SCR / "phase_7/TIS_gse25066.py", "TIS_GENES")["TIS_GENES"]
    tgfb = script_lists(SCR / "phase_8/TGFBstrom_gse25066.py", "TGFB_STROMAL_GENES")["TGFB_STROMAL_GENES"]
    dl = script_lists(SCR / "phase_7/DLDA30_gse16446.py", "DLDA30_UP", "DLDA30_DOWN")
    ggi = script_lists(SCR / "phase_6/gge_GGI_gse25066.py")
    rs = script_lists(SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")
    fw = pd.read_csv(ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv", sep="\t")
    dv = pd.read_csv(ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    ggi_up = [k for k in ggi if "UP" in k.upper() and "DOWN" not in k.upper()]
    ggi_dn = [k for k in ggi if "DOWN" in k.upper()]
    hrpos = set(pd.read_csv(ROOT / "external_data/GSE25066/GSE25066_pcr_HR+.tsv",
                            sep="\t")["patient_id"])
    print(f"gene lists: HRD={len(hrd)} (Severson), TIS={len(tis)}, TGFB={len(tgfb)}, "
          f"DLDA30={len(dl['DLDA30_UP'])}/{len(dl['DLDA30_DOWN'])}")

    def sig_of(name, expr):
        if name == "Framework":
            return singscore(expr, fw.loc[fw.direction == "sensitivity", "gene_symbol"],
                             fw.loc[fw.direction == "resistance", "gene_symbol"])
        if name == "DoxVariant":
            return singscore(expr, dv.loc[dv.direction == "sensitivity", "gene_symbol"],
                             dv.loc[dv.direction == "resistance", "gene_symbol"])
        if name == "GGI":
            return singscore(expr, ggi[ggi_up[0]], ggi[ggi_dn[0]])
        if name == "RSlike":
            return singscore(expr, rs["RS_UP"], rs["RS_DOWN"])
        if name == "DLDA30":
            u, _ = mean_z(expr, dl["DLDA30_UP"])
            d, _ = mean_z(expr, dl["DLDA30_DOWN"])
            return u - d, len(dl["DLDA30_UP"]) + len(dl["DLDA30_DOWN"])
        return mean_z(expr, {"TIS": tis, "cytolytic": ["GZMA", "PRF1"],
                             "HRD": hrd, "TGFBstrom": tgfb}[name])

    CELLS = [("Framework", "GSE32646", None), ("DoxVariant", "GSE25066", None),
             ("GGI", "GSE25066", None), ("RSlike", "GSE25066", "HR+"),
             ("DLDA30", "GSE16446", None), ("TIS", "GSE25066", None),
             ("cytolytic", "GSE22226", None), ("HRD", "GSE22226", None),
             ("TGFBstrom", "GSE25066", None)]

    cache = {}
    def cell(name, coh, stratum, drop):
        if coh not in cache:
            cache[coh] = load(coh)
        expr, ax = cache[coh]
        if stratum == "HR+":
            k = ax.index.intersection(hrpos)
            expr, ax = expr.loc[k], ax.loc[k]
        if drop:
            k = ~ax.index.isin(drop)
            expr, ax = expr.loc[k], ax.loc[k]
        s, ng = sig_of(name, expr)
        y = ax["pCR"].values.astype(int)
        Xr = np.column_stack([zscore(ax[a]).values for a in AXES])
        pt, se = residual_auroc(y, Xr, s)
        return dict(name=name, cohort=coh, stratum=stratum or "all", n=len(y),
                    events=int(y.sum()), genes=ng, point=pt, se=se)

    arms = {"baseline (overlap retained)": set(),
            "PRIMARY: remove from GSE25066": b_ov,
            "sensitivity: remove from GSE22226": a_ov}

    out = {}
    for label, drop in arms.items():
        rows = [cell(n, c, s, drop) for n, c, s in CELLS]
        out[label] = dict(cells=rows,
                          pool=fe_pool([r["point"] for r in rows], [r["se"] for r in rows],
                                       [r["cohort"] for r in rows], label))

    # ---- how much of each analysis set the overlap actually occupies ----
    base_cells = {r["name"]: r for r in out["baseline (overlap retained)"]["cells"]}
    occupancy = {}
    for coh, ids in [("GSE22226", a_ov), ("GSE25066", b_ov)]:
        _, ax = cache[coh]
        occupancy[coh] = dict(analysis_n=len(ax), overlap_in_set=int(ax.index.isin(ids).sum()))
    _, ax25 = cache["GSE25066"]
    hr = ax25.index.intersection(hrpos)
    occupancy["GSE25066 HR+ (RSlike)"] = dict(analysis_n=len(hr),
                                              overlap_in_set=int(hr.isin(b_ov).sum()))

    print(f"\noverlap: {len(shared)} shared I-SPY ids; "
          f"{len(a_ov)} GSE22226 GSMs, {len(b_ov)} GSE25066 GSMs")
    print("occupancy inside each analysis set (these differ -- do not quote one figure):")
    for k, v in occupancy.items():
        print(f"  {k:24s} {v['overlap_in_set']:3d} of {v['analysis_n']:4d}")
    print()
    base = {r["name"]: r for r in out["baseline (overlap retained)"]["cells"]}
    print(f"{'cell':12s} {'cohort':10s} | {'baseline':>16s} | {'PRIMARY (-25066)':>18s} | "
          f"{'sens (-22226)':>18s}")
    print(f"{'':12s} {'':10s} | {'n':>4s} {'resid':>10s} | {'n':>4s} {'resid':>12s} | "
          f"{'n':>4s} {'resid':>12s}")
    print("-" * 88)
    for n, c, s in CELLS:
        b = base[n]
        p = [r for r in out["PRIMARY: remove from GSE25066"]["cells"] if r["name"] == n][0]
        a = [r for r in out["sensitivity: remove from GSE22226"]["cells"] if r["name"] == n][0]
        print(f"{n:12s} {c:10s} | {b['n']:4d} {b['point']:10.4f} | "
              f"{p['n']:4d} {p['point']:12.4f} | {a['n']:4d} {a['point']:12.4f}")

    print(f"\n{'arm':38s} {'k':>2s} {'pooled':>8s} {'95% CI':>18s} | "
          + "  ".join(f"TOST<{m}" for m in MARGINS))
    print("-" * 112)
    for label, d in out.items():
        p = d["pool"]
        ci = f"[{p['ci'][0]:.4f},{p['ci'][1]:.4f}]"
        ts = "  ".join(f"{p[f'tost_{m}']:8.2e}" for m in MARGINS)
        print(f"{label:38s} {p['k']:2d} {p['pooled']:8.4f} {ci:>18s} | {ts}")

    print("\nSame pools at the rho = 1 worst case: cells within a cohort perfectly correlated,")
    print("cohorts independent. Five cells sit on the same GSE25066 patients and two on the")
    print("same GSE22226 patients, so the naive SE above is understated -- this is the")
    print("circularity the reviewer already named, and de-duplication does not fix it.")
    print(f"\n{'arm':38s} {'SE x':>5s} {'rho=1 CI':>20s} | {'<0.53':>9s} {'<0.55':>9s} "
          f"{'<0.575':>9s} {'<0.60':>9s}")
    print("-" * 106)
    for label, d in out.items():
        p = d["pool"]
        ci = f"[{p['ci_rho1'][0]:.4f},{p['ci_rho1'][1]:.4f}]"
        ts = " ".join(f"{p[f'tost_rho1_{m}']:9.2e}" for m in MARGINS)
        print(f"{label:38s} {p['rho1_over_naive']:5.2f} {ci:>20s} | {ts}")
    print("\nNOTE: the patient-level cluster bootstrap that addendum section 5.2 makes the")
    print("primary figure for any equivalence claim HAS NOT BEEN RUN. Nothing above estimates")
    print("it. The rho = 1 column is an upper bound on it, not a substitute for it.")

    # ---- strictest subset: one cell per cohort, every combination ----
    print("\nStrictest independent subset -- one cell per cohort, all combinations:")
    from itertools import product
    for label in ["PRIMARY: remove from GSE25066"]:
        rows = out[label]["cells"]
        by = {}
        for r in rows:
            by.setdefault(r["cohort"], []).append(r)
        combos = list(product(*by.values()))
        ps = []
        for cb in combos:
            pl = fe_pool([r["point"] for r in cb], [r["se"] for r in cb])
            ps.append((pl["pooled"], pl["tost_0.55"], "+".join(r["name"] for r in cb)))
        fails = [x for x in ps if x[1] >= 0.05]
        ps.sort(key=lambda x: -x[1])
        print(f"  {label}: k=4, {len(combos)} combinations, "
              f"{len(fails)} FAIL equivalence at 0.55")
        for pooled, p, nm in ps[:3]:
            print(f"    worst: {nm:52s} pooled {pooled:.4f}  p = {p:.3f}")
        out[label]["strict_subset"] = dict(n_combinations=len(combos), n_fail=len(fails),
                                           worst_p=float(ps[0][1]), worst_cells=ps[0][2])

    # ---------- retained vs removed balance on GSE22226 ----------
    _, ax22 = cache["GSE22226"]
    ph = phen22.set_index("patient_id")
    ids = [i for i in ax22.index if i in ph.index]
    rem = [i for i in ids if i in a_ov]
    ret = [i for i in ids if i not in a_ov]
    print(f"\nGSE22226 balance -- removed n={len(rem)} vs retained n={len(ret)}")
    fields = {"er_(0=negative;_1=positive;)": "ER+",
              "intrinsic_subtype_by_pam50": "PAM50",
              "pathological_complete_response_(pcr)": "pCR",
              "neoadjuvant_chemotherapy*_(1_=_ac_only;_2_=_ac_+_t;_3_=_ac_+_t_+_herceptin;_4_=_ac_+_t_+_other;)": "regimen"}
    bal = {}
    for col, lab in fields.items():
        if col not in ph.columns:
            continue
        r1 = ph.loc[rem, col].value_counts().to_dict()
        r2 = ph.loc[ret, col].value_counts().to_dict()
        bal[lab] = dict(removed=r1, retained=r2)
        print(f"  {lab:8s} removed={r1}  retained={r2}")

    json.dump(dict(shared=len(shared), gse22226_gsm=len(a_ov), gse25066_gsm=len(b_ov),
                   arms=out, gse22226_balance=bal),
              open(OUT / "A0e_dedup.json", "w"), indent=2, default=float)
    print(f"\nwrote {OUT}/A0e_dedup.json")


if __name__ == "__main__":
    main()
