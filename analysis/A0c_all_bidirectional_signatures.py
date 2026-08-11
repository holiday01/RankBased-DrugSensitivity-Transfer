#!/usr/bin/env python3
"""A0c -- is the R^2 collapse universal, or specific to the 60-gene framework?

A0b showed that correcting the scorer makes the framework's R^2 on the four axes fall
from ~0.30 to 0.03-0.09, which would remove the paper's "redundancy / aliasing" mechanism.
Before rewriting that argument, establish whether the same thing happens to every affected
signature.

Exactly five of the ten grid signatures use the buggy bidirectional compute_singscore:
    framework (60-gene Doxorubicin), Dox_variant, GGI, MammaPrint, RSlike
The other five are unaffected and are NOT re-run here:
    cytolytic, TIS, TGFBstrom  -- plain mean z-score, all genes one direction
    DLDA30                     -- mean(UP z) - mean(DOWN z); down is NOT rank-reversed,
                                  so the subtraction is correct
    HRD/Fanconi                -- z-score based (verified separately)

All five are scored on GSE25066 (n=488, the largest cohort, and the one carrying ~4/5 of
the pooled weight) so the comparison is like-for-like. Gene lists are extracted from the
deployed scripts by AST parsing, so no script is executed and no constant is retyped.
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
OUT = _OR / "A0"
EXPR = ROOT / "external_data/GSE25066/GSE25066_expression.tsv"
AXES_F = ROOT / "resubmission_v2/results/phase_0/axes_GSE25066.tsv"
AXES = ["P", "H", "B", "H_PAM50basal"]


def const_lists(path, names):
    """Pull module-level list constants out of a script without importing it."""
    tree = ast.parse(Path(path).read_text())
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id in names:
                    try:
                        found[t.id] = [str(x) for x in ast.literal_eval(node.value)]
                    except Exception:
                        pass
    return found


def load_signatures():
    sigs = {}

    # 1. 60-gene framework and the Dox variant come from the frozen deposit.
    fw = pd.read_csv(ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv", sep="\t")
    sigs["framework"] = (fw.loc[fw.direction == "sensitivity", "gene_symbol"].tolist(),
                         fw.loc[fw.direction == "resistance", "gene_symbol"].tolist())

    # Dox_variant is the quantile_20_protein build, per phase_6/Dox_variant_gse25066.py:39
    pan = ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv"
    if not pan.exists():
        hits = list((ROOT / "frozen_signatures").rglob("quantile_20_protein/DOXORUBICIN.tsv"))
        pan = hits[0] if hits else pan
    if pan.exists():
        dv = pd.read_csv(pan, sep="\t")
        sigs["Dox_variant"] = (dv.loc[dv.direction == "sensitivity", "gene_symbol"].tolist(),
                               dv.loc[dv.direction == "resistance", "gene_symbol"].tolist())
    else:
        print("  ! Dox_variant: quantile_20_protein/DOXORUBICIN.tsv not found -- skipped")

    # MammaPrint is stored as a {gene: +1/-1} weight dict, not two lists.
    mp = SCR / "phase_6/mammaprint_metabric.py"
    if mp.exists():
        for node in ast.parse(mp.read_text()).body:
            if isinstance(node, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == "MAMMAPRINT" for t in node.targets):
                try:
                    w = ast.literal_eval(node.value)
                    sigs["MammaPrint"] = ([g for g, v in w.items() if v > 0],
                                          [g for g, v in w.items() if v < 0])
                except Exception as e:
                    print(f"  ! MammaPrint: {e}")

    # 2. GGI and RSlike are plain gene lists inside the deployed scripts.
    for key, path, upn, dnn in [
        ("GGI", SCR / "phase_6/gge_GGI_gse25066.py", None, None),
        ("RSlike", SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN"),
    ]:
        if not path.exists():
            continue
        tree = ast.parse(path.read_text())
        cands = {}
        for node in tree.body:
            if isinstance(node, ast.Assign):
                for t in node.targets:
                    if isinstance(t, ast.Name):
                        try:
                            v = ast.literal_eval(node.value)
                            if isinstance(v, list) and v and all(isinstance(x, str) for x in v):
                                cands[t.id] = v
                        except Exception:
                            pass
        if upn and dnn and upn in cands and dnn in cands:
            sigs[key] = (cands[upn], cands[dnn])
            continue
        up = [k for k in cands if "UP" in k.upper() and "DOWN" not in k.upper()]
        dn = [k for k in cands if "DOWN" in k.upper()]
        if up and dn:
            sigs[key] = (cands[up[0]], cands[dn[0]])
        else:
            print(f"  ! {key}: could not identify up/down lists in {path.name} "
                  f"(candidates: {sorted(cands)})")
    return sigs


def zscore(s):
    s = pd.Series(s).astype(float)
    return (s - s.mean()) / s.std(ddof=0)


def score(expr, up, down, variant):
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)
    U = [g for g in up if g in expr.columns]
    D = [g for g in down if g in expr.columns]
    if len(U) < 3 or len(D) < 3:
        return None, len(U), len(D)

    def norm(m, n):
        mn = (1 + n) / (2 * n_total)
        return 2 * (m - mn) / ((1 - mn) - mn) - 1

    up_s = norm(ranks[U].mean(axis=1), len(U))
    dn_s = norm((1 - ranks[D]).mean(axis=1), len(D))
    raw = (up_s - dn_s) if variant == "deployed" else (up_s + dn_s)
    return raw - raw.median(), len(U), len(D)


def analyse(y, Xr, s):
    s = zscore(s).values
    A = sm.add_constant(Xr)
    resid = s - A @ np.linalg.lstsq(A, s, rcond=None)[0]
    mr = sm.Logit(y, sm.add_constant(Xr)).fit(disp=0)
    mf = sm.Logit(y, sm.add_constant(np.column_stack([Xr, s]))).fit(disp=0)
    return dict(marginal=roc_auc_score(y, s), residual=roc_auc_score(y, resid),
                r2=1 - resid.var() / s.var(),
                lr_p=float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("A0c -- all bidirectional-singscore signatures on GSE25066 (n=488)\n")
    sigs = load_signatures()
    print(f"signatures recovered: {', '.join(sigs)}\n")

    expr = pd.read_csv(EXPR, sep="\t", index_col=0)
    if expr.shape[1] < 5000:
        expr = expr.T
    ax = pd.read_csv(AXES_F, sep="\t", index_col=0).dropna(subset=AXES + ["pCR"])
    common = expr.index.intersection(ax.index)
    expr, ax = expr.loc[common], ax.loc[common]
    y = ax["pCR"].values.astype(int)
    Xr = np.column_stack([zscore(ax[a]).values for a in AXES])

    rows = []
    for name, (up, down) in sigs.items():
        sd, nU, nD = score(expr, up, down, "deployed")
        sc, _, _ = score(expr, up, down, "corrected")
        if sd is None:
            print(f"  ! {name}: too few genes resolved ({nU} up, {nD} down) -- skipped")
            continue
        d, c = analyse(y, Xr, sd), analyse(y, Xr, sc)
        rows.append(dict(signature=name, n_up=nU, n_down=nD,
                         **{f"dep_{k}": v for k, v in d.items()},
                         **{f"cor_{k}": v for k, v in c.items()},
                         r_dep_cor=float(np.corrcoef(sd, sc)[0, 1])))

    json.dump(rows, open(OUT / "A0c_bidirectional_signatures.json", "w"), indent=2, default=float)

    print(f"{'signature':13s} {'genes':>11s} | {'R2 on 4 axes':>15s} | {'residual AUROC':>16s} | "
          f"{'marginal':>15s} | {'LRT p':>15s}")
    print(f"{'':13s} {'up/down':>11s} | {'dep':>7s} {'cor':>7s} | {'dep':>7s} {'cor':>8s} | "
          f"{'dep':>7s} {'cor':>7s} | {'dep':>7s} {'cor':>7s}")
    print("-" * 106)
    for r in rows:
        print(f"{r['signature']:13s} {r['n_up']:5d}/{r['n_down']:<5d} | "
              f"{r['dep_r2']:7.3f} {r['cor_r2']:7.3f} | "
              f"{r['dep_residual']:7.4f} {r['cor_residual']:8.4f} | "
              f"{r['dep_marginal']:7.4f} {r['cor_marginal']:7.4f} | "
              f"{r['dep_lr_p']:7.4f} {r['cor_lr_p']:7.4f}")

    print("\ncorrelation between the deployed and corrected score, per signature:")
    for r in rows:
        print(f"  {r['signature']:13s} r = {r['r_dep_cor']:+.4f}")

    up_r2 = [r for r in rows if r["cor_r2"] > r["dep_r2"]]
    dn_r2 = [r for r in rows if r["cor_r2"] < r["dep_r2"]]
    print(f"\nR2 on the axes FELL after correction for : "
          f"{', '.join(r['signature'] for r in dn_r2) or 'none'}")
    print(f"R2 on the axes ROSE after correction for : "
          f"{', '.join(r['signature'] for r in up_r2) or 'none'}")
    print(f"\nresidual AUROC >= 0.55 after correction  : "
          f"{', '.join(r['signature'] for r in rows if r['cor_residual'] >= 0.55) or 'none'}")
    print(f"nested-LRT p < 0.05 after correction     : "
          f"{', '.join(r['signature'] for r in rows if r['cor_lr_p'] < 0.05) or 'none'}")
    print(f"\nwrote {OUT}/A0c_bidirectional_signatures.json")


if __name__ == "__main__":
    main()
