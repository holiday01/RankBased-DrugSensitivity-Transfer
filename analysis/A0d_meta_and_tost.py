#!/usr/bin/env python3
"""A0d -- rebuild the residual-AUROC meta-pool and the TOST equivalence claim.

Three defects are corrected together, because they all land on the same pooled number:

 1. SCORER. Five of the eleven cells use the bidirectional singscore that A0 showed to be
    sign-inverted (framework, Dox_variant, MammaPrint, GGI, RSlike). Those are re-scored.
    The other six (DLDA30, TIS, cytolytic, HRD, TGFBstrom) use z-score based scorers that
    are unaffected, and are carried over from the deposit unchanged.

 2. METABRIC DOUBLE COUNT. The deposited meta_cells.tsv carries MammaPrint on METABRIC
    TWICE -- "Mammaprint*" (n=399) and "Mammaprint_metabric" (n=0). The k=11 pool counts
    the same 399 patients in two cells.

 3. PATIENT OVERLAP. GSE22226 and GSE25066 share 65 patients (58 inside the analysis
    sets), so cells drawn from those two cohorts are not independent.

The pool is reported four ways so the reader can see what each assumption buys:
    naive fixed-effect        -- as published, every cell independent
    METABRIC de-duplicated    -- the double count removed
    cohort-clustered          -- variance inflated for cells sharing a cohort
    one-cell-per-cohort       -- the strictest independent subset

TOST is run across the pre-specified margin sensitivity curve {0.53, 0.55, 0.575, 0.60}.
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
SCR = ROOT / "resubmission_v2/scripts"
RES = ROOT / "resubmission_v2/results"
OUT = _OR / "A0"
AXES = ["P", "H", "B", "H_PAM50basal"]
MARGINS = [0.53, 0.55, 0.575, 0.60]
NBOOT = 1000
SEED = 42


def zscore(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def sing(expr, up, down, corrected=True):
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in expr.columns]
    D = [g for g in down if g in expr.columns]

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1

    up_s = norm(ranks[U].mean(axis=1), len(U))
    dn_s = norm((1 - ranks[D]).mean(axis=1), len(D))
    raw = (up_s + dn_s) if corrected else (up_s - dn_s)
    return raw - raw.median()


def residual_auroc(y, Xr, s, nboot=NBOOT, seed=SEED):
    """Residual AUROC with a bootstrap SE. The outcome never enters the OLS step."""
    s = zscore(s).values
    A = sm.add_constant(Xr)

    def one(idx):
        Ai, si, yi = A[idx], s[idx], y[idx]
        if len(np.unique(yi)) < 2:
            return np.nan
        r = si - Ai @ np.linalg.lstsq(Ai, si, rcond=None)[0]
        return roc_auc_score(yi, r)

    point = one(np.arange(len(y)))
    rng = np.random.default_rng(seed)
    boots = np.array([one(rng.integers(0, len(y), len(y))) for _ in range(nboot)])
    boots = boots[~np.isnan(boots)]
    return float(point), float(boots.std(ddof=1)), [float(np.percentile(boots, 2.5)),
                                                    float(np.percentile(boots, 97.5))]


def load_sigs():
    import ast
    sigs = {}
    fw = pd.read_csv(ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv", sep="\t")
    sigs["Framework"] = (fw.loc[fw.direction == "sensitivity", "gene_symbol"].tolist(),
                         fw.loc[fw.direction == "resistance", "gene_symbol"].tolist())
    dv = pd.read_csv(ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    sigs["DoxVariant"] = (dv.loc[dv.direction == "sensitivity", "gene_symbol"].tolist(),
                          dv.loc[dv.direction == "resistance", "gene_symbol"].tolist())
    for key, path, upn, dnn in [("GGI", SCR / "phase_6/gge_GGI_gse25066.py", None, None),
                                ("RSlike", SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")]:
        cands = {}
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        try:
                            v = ast.literal_eval(node.value)
                            if isinstance(v, list) and v and all(isinstance(x, str) for x in v):
                                cands[t.id] = v
                        except Exception:
                            pass
        if upn and upn in cands:
            sigs[key] = (cands[upn], cands[dnn])
        else:
            u = [k for k in cands if "UP" in k.upper() and "DOWN" not in k.upper()]
            d = [k for k in cands if "DOWN" in k.upper()]
            sigs[key] = (cands[u[0]], cands[d[0]])
    for node in ast.parse((SCR / "phase_6/mammaprint_metabric.py").read_text()).body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "MAMMAPRINT" for t in node.targets):
            w = ast.literal_eval(node.value)
            sigs["MammaPrint"] = ([g for g, v in w.items() if v > 0],
                                  [g for g, v in w.items() if v < 0])
    return sigs


def load_cohort(name):
    paths = {
        "GSE25066": (ROOT / "external_data/GSE25066/GSE25066_expression.tsv",
                     RES / "phase_0/axes_GSE25066.tsv"),
        "GSE32646": (RES / "phase_2/GSE32646_expression.tsv", RES / "phase_2/axes_GSE32646.tsv"),
        "METABRIC": (RES / "phase_4/metabric/METABRIC_expression.tsv",
                     RES / "phase_4/metabric/axes_METABRIC.tsv"),
    }
    ep, ap = paths[name]
    expr = pd.read_csv(ep, sep="\t", index_col=0)
    if expr.shape[1] < 5000:
        expr = expr.T
    ax = pd.read_csv(ap, sep="\t", index_col=0)
    ycol = "pCR" if "pCR" in ax.columns else [c for c in ax.columns if "pcr" in c.lower() or "event" in c.lower()][0]
    ax = ax.dropna(subset=[c for c in AXES if c in ax.columns] + [ycol])
    common = expr.index.intersection(ax.index)
    return expr.loc[common], ax.loc[common], ycol


def fe_pool(points, ses, label, cluster=None):
    """Inverse-variance fixed effect. If cluster given, inflate SEs for shared cohorts."""
    p = np.asarray(points, float)
    se = np.asarray(ses, float)
    if cluster is not None:
        cl = np.asarray(cluster)
        for c in np.unique(cl):
            m = cl == c
            if m.sum() > 1:
                se[m] = se[m] * np.sqrt(m.sum())     # cells share the same patients
    w = 1 / se ** 2
    pooled = float((w * p).sum() / w.sum())
    se_pool = float(np.sqrt(1 / w.sum()))
    ci = [pooled - 1.96 * se_pool, pooled + 1.96 * se_pool]
    Q = float((w * (p - pooled) ** 2).sum())
    df = len(p) - 1
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    # One-sided equivalence: H0 theta >= margin, H1 theta < margin.
    # p = P(observe this or lower | theta = margin) = Phi((pooled - margin)/se).
    tost = {f"tost_{m}": float(stats.norm.cdf((pooled - m) / se_pool)) for m in MARGINS}
    return dict(label=label, k=len(p), pooled=round(pooled, 4),
                ci=[round(ci[0], 4), round(ci[1], 4)], se=round(se_pool, 5),
                I2=round(I2, 1), **{k: v for k, v in tost.items()})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sigs = load_sigs()
    dep = pd.read_csv(RES / "phase_9_revalidation/meta_cells.tsv", sep="\t")

    # ---- cells that must be re-scored (signature, cohort, stratum) ----
    todo = [("Framework", "GSE32646", None), ("DoxVariant", "GSE25066", None),
            ("GGI", "GSE25066", None), ("RSlike", "GSE25066", "HR+"),
            ("MammaPrint", "METABRIC", None)]

    hrpos = set(pd.read_csv(ROOT / "external_data/GSE25066/GSE25066_pcr_HR+.tsv",
                            sep="\t")["patient_id"])

    cells = []
    cache = {}
    for signame, coh, stratum in todo:
        if coh not in cache:
            cache[coh] = load_cohort(coh)
        expr, ax, ycol = cache[coh]
        if stratum == "HR+":
            keep = ax.index.intersection(hrpos)
            expr, ax = expr.loc[keep], ax.loc[keep]
        y = ax[ycol].values.astype(int)
        Xr = np.column_stack([zscore(ax[a]).values for a in AXES])
        up, down = sigs[signame]
        pt_c, se_c, ci_c = residual_auroc(y, Xr, sing(expr, up, down, True))
        pt_d, se_d, _ = residual_auroc(y, Xr, sing(expr, up, down, False))
        cells.append(dict(name=signame, cohort=coh, stratum=stratum or "all",
                          n=len(y), events=int(y.sum()), affected=True,
                          point=round(pt_c, 4), se=round(se_c, 5),
                          ci=[round(ci_c[0], 4), round(ci_c[1], 4)],
                          deployed_point=round(pt_d, 4)))
        print(f"  re-scored {signame:11s} x {coh:9s} {stratum or '':4s} "
              f"n={len(y):4d}  residual {pt_d:.4f} -> {pt_c:.4f}")

    # ---- unaffected cells, carried over from the deposit ----
    UNAFFECTED = {"DLDA30_gse16446": "GSE16446", "TIS_gse25066": "GSE25066",
                  "cytolytic_gse22226": "GSE22226", "HRD_gse22226": "GSE22226",
                  "TGFBstrom_gse25066": "GSE25066"}
    for _, r in dep.iterrows():
        if r["name"] in UNAFFECTED:
            cells.append(dict(name=r["name"], cohort=UNAFFECTED[r["name"]], stratum="all",
                              n=int(r["n"]) if pd.notna(r["n"]) else 0, events=None,
                              affected=False, point=round(float(r["point"]), 4),
                              se=round(float(r["se"]), 5), ci=None, deployed_point=None))

    df = pd.DataFrame(cells)
    df.to_csv(OUT / "A0d_meta_cells_corrected.tsv", sep="\t", index=False)

    print(f"\n{len(df)} cells "
          f"({df.affected.sum()} re-scored, {(~df.affected).sum()} carried over unchanged)")
    print("METABRIC MammaPrint double count in the deposit: removed (was 2 cells, now 1)\n")

    pools = []
    pools.append(fe_pool(df.point, df.se, f"ALL cells, naive (k={len(df)})"))
    pcr = df[df.cohort != "METABRIC"]
    pools.append(fe_pool(pcr.point, pcr.se, f"pCR-endpoint only, naive (k={len(pcr)})"))
    pools.append(fe_pool(pcr.point, pcr.se, f"pCR-endpoint, cohort-clustered (k={len(pcr)})",
                         cluster=pcr.cohort))
    one = pcr.sort_values("n", ascending=False).groupby("cohort", as_index=False).first()
    pools.append(fe_pool(one.point, one.se, f"one cell per cohort (k={len(one)})"))

    print(f"{'pool':44s} {'k':>3s} {'pooled':>8s} {'95% CI':>18s} {'I2':>5s} | "
          + " ".join(f"TOST<{m}" for m in MARGINS))
    print("-" * 118)
    for p in pools:
        ci = f"[{p['ci'][0]:.4f},{p['ci'][1]:.4f}]"
        ts = " ".join(f"{p[f'tost_{m}']:9.2e}" for m in MARGINS)
        print(f"{p['label']:44s} {p['k']:3d} {p['pooled']:8.4f} {ci:>18s} {p['I2']:5.1f} | {ts}")

    json.dump(dict(cells=cells, pools=pools), open(OUT / "A0d_pools.json", "w"),
              indent=2, default=float)

    print("\nEquivalence verdict (TOST p < 0.05 => equivalence to <margin> is established):")
    for p in pools:
        ok = [str(m) for m in MARGINS if p[f"tost_{m}"] < 0.05]
        fail = [str(m) for m in MARGINS if p[f"tost_{m}"] >= 0.05]
        print(f"  {p['label']:44s} holds at: {', '.join(ok) or 'NONE'}"
              + (f"   fails at: {', '.join(fail)}" if fail else ""))

    print(f"\nPublished claim was: pooled 0.5015 [0.4742, 0.5288], k=9, TOST<0.55 p=2.5e-4")
    print(f"wrote {OUT}/A0d_pools.json and A0d_meta_cells_corrected.tsv")


if __name__ == "__main__":
    main()
