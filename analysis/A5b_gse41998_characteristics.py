#!/usr/bin/env python3
"""A5b -- the cohort-characteristics row for GSE41998 (R3-3).

A5 predates GSE41998, so its deposited table covers five cohorts and not this one. R3-3
asked for subtype counts per cohort; supplying five of six would answer the comment only
partly. This computes the missing row from the phenotype file rather than transcribing the
prose already in section 2.2, so the published table and the text have one source.

Counts follow A5's own conventions: HER2 missing is treated as negative, matching the
deposited `A4.strata_of`, and the analysis set is the one the four-axis model can be fitted
on -- the same set A14b uses.

Writes results/A5/A5b_gse41998_row.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from A14b_gse41998_analysis import G                                       # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "results" / "A5"


def main() -> None:
    ax = pd.read_csv(G / "GSE41998_axes.tsv", sep="\t", index_col=0)
    expr = pd.read_csv(G / "GSE41998_expression.tsv", sep="\t", index_col=0)
    ph = pd.read_csv(G / "GSE41998_phenotype.tsv", sep="\t", index_col=0)
    keep = ax.index.intersection(expr.index)
    p = ph.reindex(keep)

    low = lambda col: p[col].astype(str).str.strip().str.lower()
    er = low("er_status_ihc").eq("positive")
    pr = low("pr_status_ihc").eq("positive")
    her2 = low("her2_status").eq("positive")          # missing -> not positive, as in A4
    hr = er | pr
    y = ax.loc[keep, "pCR"].astype(int)

    row = {
        "cohort": "GSE41998",
        "n_profiled": int(len(ph)),
        "n_analysed": int(len(keep)),
        "events": int(y.sum()),
        "pcr_rate": round(float(y.mean()), 3),
        "ER_pos": int(er.sum()),
        "ER_neg": int((~er).sum()),
        "HER2_pos": int(her2.sum()),
        "HER2_neg": int((~her2).sum()),
        "TNBC": int((~hr & ~her2).sum()),
        "HRpos_HER2neg": int((hr & ~her2).sum()),
        "conventions": ("HER2 missing counted negative, as in the deposited A4.strata_of; "
                        "TNBC = ER-negative and PR-negative and HER2-negative; "
                        "analysis set = samples with both an expression profile and all four axes"),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "A5b_gse41998_row.json").write_text(json.dumps(row, indent=1))
    for k, v in row.items():
        print(f"  {k:16s} {v}")


if __name__ == "__main__":
    main()
