#!/usr/bin/env python3
"""A6b -- the endpoint-substitution check (R3-6) on GSE41998.

**Specified and run at revision stage. Not pre-registered.** The addendum's A6 harmonises
the pCR endpoint to RCB-0/I on the cohorts that deposit an RCB class, but GSE41998 was not
in the study when A6 was written. The manuscript stated that GSE41998 deposits an RCB-0/I
call on the same 253 patients and then declined to run the check, calling it the natural
next step; Reviewer 3's point 6 asked for exactly this, and the data are in the deposit,
so declining costs more than running it.

Why GSE41998 is the cleanest place to ask. The other cohorts carrying an RCB field are
GSE22226 and GSE25066, which share 65 I-SPY 1 patients, so an endpoint comparison there is
confounded with the overlap. GSE41998 has zero overlap with any cohort in the study, and
the two endpoints are measured on the *same* patients, so the only thing that changes
between the two fits is the definition of a response.

Everything else is held to A14b: the same four axes, the same frozen framework signature,
the same scorer, and patient-level bootstrap percentile intervals at B = 2000, seed 42.

Writes results/A6/A6b_gse41998_endpoint.json.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import A0e_dedup as D                                                     # noqa: E402
from A2_grid import build_defs, score                                     # noqa: E402
from A0i_cluster_bootstrap import _cell_builders                          # noqa: E402
from A14b_gse41998_analysis import cell, boot_ci, G, AXES                  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "results" / "A6"


def main() -> None:
    ax = pd.read_csv(G / "GSE41998_axes.tsv", sep="\t", index_col=0)
    expr = pd.read_csv(G / "GSE41998_expression.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(G / "GSE41998_phenotype.tsv", sep="\t", index_col=0)
    common = ax.index.intersection(expr.index)
    ax, expr = ax.loc[common], expr.loc[common]

    defs = build_defs(_cell_builders())
    up, dn, scorer = defs["Framework"][:3]
    s_all = D.zscore(score(expr, up, dn, scorer)[0]).values
    Xr_all = np.column_stack([D.zscore(ax[a]).values for a in AXES])

    rcb = ph["pcr_rcb01_raw"].reindex(ax.index).astype(str).str.strip().str.lower()
    endpoints = {
        "pCR (deposited call, as used throughout)": ax["pCR"].values.astype(int),
        "RCB-0/I (harmonised, R3-6)": rcb.map({"yes": 1, "no": 0}).values,
    }

    res = {
        "cohort": "GSE41998",
        "provenance": "SPECIFIED AT REVISION STAGE -- not pre-registered; see decisions_addendum entry 13",
        "question": "R3-6: does the conclusion depend on how response is defined?",
        "why_this_cohort": ("GSE22226 and GSE25066 share 65 I-SPY 1 patients, so an endpoint "
                            "comparison on them is confounded with the overlap. GSE41998 has "
                            "zero overlap and carries both calls on the same patients."),
        "bootstrap": {"B": 2000, "seed": 42, "kind": "patient-level percentile"},
        "endpoints": {},
    }

    for label, y in endpoints.items():
        keep = ~pd.isna(y)
        yi = np.asarray(y[keep], dtype=int)
        c = cell(yi, Xr_all[keep], s_all[keep])
        ci, se = boot_ci(yi, Xr_all[keep], s_all[keep])
        c.update(residual_ci=ci, residual_se=se, n_dropped_missing=int((~keep).sum()))
        res["endpoints"][label] = c
        print(f"{label:44s} n={c['n']:3d} events={c['events']:3d} "
              f"marg={c['marginal']:.3f} resid={c['residual']:.3f} "
              f"[{ci[0]:.3f},{ci[1]:.3f}] LR P={c['lr_p']:.4f}")

    a, b = list(res["endpoints"].values())
    res["comparison"] = {
        "residual_delta_rcb_minus_pcr": round(b["residual"] - a["residual"], 4),
        "both_null_at_0.05": bool(a["lr_p"] >= 0.05 and b["lr_p"] >= 0.05),
        "reading": ("The verdict does not depend on the endpoint definition."
                    if a["lr_p"] >= 0.05 and b["lr_p"] >= 0.05 else
                    "The two endpoints disagree; the endpoint definition is load-bearing "
                    "and must be reported as such."),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "A6b_gse41998_endpoint.json").write_text(json.dumps(res, indent=1))
    print("\n" + res["comparison"]["reading"])
    print(f"wrote {OUT / 'A6b_gse41998_endpoint.json'}")


if __name__ == "__main__":
    main()
