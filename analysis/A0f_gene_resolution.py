#!/usr/bin/env python3
"""A0f -- gene-resolution fraction for all 58 Family A cells.

Required by ANALYSIS_ADDENDUM_2026-08.md §12. §4.2 trigger 2 marks a cell not evaluable when
the resolved fraction of its up-set OR its down-set, taken separately, falls below 0.60.
Recording every fraction at lock time is what stops the not-evaluable set being revisited
after the results are seen.

Uses no outcome data. Gene symbols are matched exactly, then through the same HGNC alias
table the deployed scorers use, so the fractions here are the ones the analysis will see.

Output: results/A0/gene_resolution_58.tsv
"""
import ast
from pathlib import Path

import pandas as pd

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
SIG_DIR = ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all"

COHORTS = {
    "GSE16446": ROOT / "external_data/GSE16446/GSE16446_expression.tsv",
    "GSE22226": ROOT / "external_data/GSE22226/GSE22226_expression.tsv",
    "GSE25066": ROOT / "external_data/GSE25066/GSE25066_expression.tsv",
    "GSE32646": RES / "phase_2/GSE32646_expression.tsv",
}
DISCOVERY = ["GSE16446", "GSE22226", "GSE25066"]
DRUGS = ["CARBOPLATIN", "CYCLOPHOSPHAMIDE", "FLUOROURACIL", "PACLITAXEL",
         "ETOPOSIDE", "TOP2_POISON_CONSENSUS"]

_alias = None


def alias_map():
    global _alias
    if _alias is None:
        _alias = {}
        p = ROOT / "external_data/hgnc_alias_table.tsv"
        if p.exists():
            df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
            c = {x.lower(): x for x in df.columns}
            if "alias_symbol" in c and "canonical_symbol" in c:
                for _, r in df.iterrows():
                    a, k = r[c["alias_symbol"]].upper(), r[c["canonical_symbol"]].upper()
                    if a and k:
                        _alias[a] = k
    return _alias


def resolve(genes, columns):
    up = {x.upper(): x for x in columns}
    am = alias_map()
    hit = []
    for g in genes:
        gU = str(g).upper()
        if gU in up or (gU in am and am[gU] in up):
            hit.append(gU)
    return len(set(hit))


def script_lists(path, *names):
    out = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and (not names or t.id in names):
                    try:
                        v = ast.literal_eval(node.value)
                        if isinstance(v, (list, dict)) and v:
                            out[t.id] = v
                    except Exception:
                        pass
    return out


def frozen(name):
    d = pd.read_csv(SIG_DIR / f"{name}.tsv", sep="\t")
    return (d.loc[d.direction == "sensitivity", "gene_symbol"].tolist(),
            d.loc[d.direction == "resistance", "gene_symbol"].tolist())


def signatures():
    s = {}
    s["Framework"] = frozen("DOXORUBICIN")
    dv = pd.read_csv(ROOT / "frozen_signatures/quantile_20_protein/DOXORUBICIN.tsv", sep="\t")
    s["DoxVariant"] = (dv.loc[dv.direction == "sensitivity", "gene_symbol"].tolist(),
                       dv.loc[dv.direction == "resistance", "gene_symbol"].tolist())
    mp = script_lists(SCR / "phase_6/mammaprint_metabric.py", "MAMMAPRINT")["MAMMAPRINT"]
    s["MammaPrint"] = ([g for g, v in mp.items() if v > 0], [g for g, v in mp.items() if v < 0])
    g = script_lists(SCR / "phase_6/gge_GGI_gse25066.py")
    gu = [k for k in g if "UP" in k.upper() and "DOWN" not in k.upper()]
    gd = [k for k in g if "DOWN" in k.upper()]
    s["GGI"] = (g[gu[0]], g[gd[0]])
    r = script_lists(SCR / "phase_6/RSlike_gse25066_HRpos.py", "RS_UP", "RS_DOWN")
    s["RSlike"] = (r["RS_UP"], r["RS_DOWN"])
    d = script_lists(SCR / "phase_7/DLDA30_gse16446.py", "DLDA30_UP", "DLDA30_DOWN")
    s["DLDA30"] = (d["DLDA30_UP"], d["DLDA30_DOWN"])
    s["TIS"] = (script_lists(SCR / "phase_7/TIS_gse25066.py", "TIS_GENES")["TIS_GENES"], [])
    s["cytolytic"] = (["GZMA", "PRF1"], [])
    s["HRD"] = (script_lists(SCR / "phase_8/HRD_gse22226.py",
                             "HRD_GENES_SEVERSON")["HRD_GENES_SEVERSON"], [])
    s["TGFBstrom"] = (script_lists(SCR / "phase_8/TGFBstrom_gse25066.py",
                                   "TGFB_STROMAL_GENES")["TGFB_STROMAL_GENES"], [])
    for d_ in DRUGS:
        s[d_] = frozen(d_)
    return s


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sigs = signatures()
    cols = {c: set(pd.read_csv(p, sep="\t", index_col=0, nrows=0).columns)
            if pd.read_csv(p, sep="\t", index_col=0, nrows=1).shape[1] > 5000
            else set(pd.read_csv(p, sep="\t", index_col=0, usecols=[0]).index)
            for c, p in COHORTS.items()}

    rows = []
    plan = ([(s, c) for s in ["Framework", "DoxVariant", "MammaPrint", "GGI", "RSlike",
                              "DLDA30", "TIS", "cytolytic", "HRD", "TGFBstrom"]
             for c in COHORTS]
            + [(d, c) for d in DRUGS for c in DISCOVERY])

    for sig, coh in plan:
        up, dn = sigs[sig]
        nu, nd = len(set(g.upper() for g in up)), len(set(g.upper() for g in dn))
        ru, rd = resolve(up, cols[coh]), resolve(dn, cols[coh])
        fu = ru / nu if nu else float("nan")
        fd = rd / nd if nd else float("nan")
        worst = min([x for x in (fu, fd) if x == x], default=float("nan"))
        rows.append(dict(signature=sig, cohort=coh, n_up=nu, n_down=nd,
                         resolved_up=ru, resolved_down=rd,
                         frac_up=round(fu, 4) if fu == fu else "",
                         frac_down=round(fd, 4) if fd == fd else "",
                         worst_frac=round(worst, 4),
                         not_evaluable_trigger2=("YES" if worst < 0.60 else "no")))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gene_resolution_58.tsv", sep="\t", index=False)

    print(f"Family A: {len(df)} cells (expected 58)\n")
    print(df.to_string(index=False))
    flagged = df[df.not_evaluable_trigger2 == "YES"]
    print(f"\ncells failing trigger 2 (worst set resolution < 0.60): {len(flagged)}")
    if len(flagged):
        print(flagged[["signature", "cohort", "worst_frac"]].to_string(index=False))
    print(f"\nwrote {OUT}/gene_resolution_58.tsv")


if __name__ == "__main__":
    main()
