"""A14 BH family, applied as section 7 of the addendum specifies it.

Why this file exists rather than an edit to A14b. `analysis/A14b_gse41998_analysis.py`
is listed in the LOCK's deposit_manifest by sha256, so editing it would break the
verification a third party is told to run. Its output
`results/A14/A14b_analysis.json` carries `family = "outside Family A ...; nominal p,
no BH"`, which contradicts the locked addendum. The correction is therefore written
as a new artefact and recorded in the deviation log; the pinned script and its output
are left byte-identical.

What the locked addendum actually says, section 7, verbatim:

    If it passes, **A14 forms its own BH family**: its cell list and integer m are
    written to `decisions_addendum.md` before its first p value is computed, and
    section 10 conditions 1 and 6 apply to its cells exactly as to Family A cells.

So a family was pre-specified. Two things follow, and they point in opposite
directions, so both are reported:

  - m is not a free parameter. Section 4.2 fixes the signature list at ten, and
    confines the six-drug CCLE arm to the three discovery cohorts, which GSE41998 is
    not. Ten signatures on one cohort is the family's maximum extent under section
    4.2's own enumeration rule, so m = 10 with no researcher latitude.

  - The recording requirement was not met. No entry in decisions_addendum.md names
    this cell list or this m, and entry 5 asserted the opposite of section 7. The
    deviation is real and is logged; what it does not do is make m adjustable, since
    the lists that determine it were locked.

Folding these cells into Family A at m = 68 is the one option the locked document
positively forbids: section 4.2 says "m = 58 is a fixed denominator" and "No
substitution into or out of these lists after this document is locked".
"""

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "results" / "A14" / "A14b_analysis.json"
OUT = ROOT / "results" / "A14" / "A14c_bh_family.json"

Q_THRESHOLD = 0.10  # section 4.2

# Section 4.2's signature list, verbatim and in its order. m is derived from this
# list rather than from whatever happens to be in the results file.
SIGNATURES_4_2 = [
    "Framework", "DoxVariant", "MammaPrint", "GGI", "RSlike",
    "DLDA30", "TIS", "cytolytic", "HRD", "TGFBstrom",
]


def benjamini_hochberg(pvalues):
    """Step-up. Returns q in the input order, monotone non-decreasing in rank."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    q = [None] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        q[i] = running
    return q


def verdict_class(cell, q):
    """Section 6's six classes, evaluated in the listed order, first match wins."""
    up, down = cell["frac_up_resolved"], cell["frac_down_resolved"]
    resolved = [f for f in (up, down) if f is not None]
    if min(cell["events"], cell["n"] - cell["events"]) < 5 or any(f < 0.60 for f in resolved):
        return "NOT-EVALUABLE"
    lo, hi = cell["residual_ci"]
    if q < Q_THRESHOLD and lo > 0.5:
        return "ADDS-VALUE"
    if hi < 0.5:
        return "INVERSE-DIRECTION"
    if cell["marginal"] >= 0.60 and cell["residual"] <= 0.55 and lo <= 0.5 <= hi \
            and cell["r2_on_axes"] >= 0.50:
        return "REDUNDANT-WITH-AXES"
    if cell["marginal"] < 0.60:
        return "NO-MARGINAL-SIGNAL"
    return "INDETERMINATE"


def main():
    src = json.loads(SRC.read_text())
    cells = src["cells"]

    missing = [s for s in SIGNATURES_4_2 if s not in cells]
    extra = [s for s in cells if s not in SIGNATURES_4_2]
    if missing or extra:
        raise SystemExit(
            f"cell list does not match section 4.2: missing={missing} extra={extra}. "
            "m is derived from the locked list, so this must be reconciled, not coerced."
        )

    m = len(SIGNATURES_4_2)
    # A cell that is not evaluable enters BH at p = 1 with the denominator unreduced
    # (section 4.2). None is, here, but the rule is applied rather than assumed.
    pvals = []
    for s in SIGNATURES_4_2:
        c = cells[s]
        resolved = [f for f in (c["frac_up_resolved"], c["frac_down_resolved"]) if f is not None]
        not_evaluable = (min(c["events"], c["n"] - c["events"]) < 5
                         or any(f < 0.60 for f in resolved)
                         or not c["converged"])
        pvals.append(1.0 if not_evaluable else c["lr_p"])

    qs = benjamini_hochberg(pvals)

    out_cells = {}
    for s, p, q in zip(SIGNATURES_4_2, pvals, qs):
        c = cells[s]
        out_cells[s] = {
            "lr_p": c["lr_p"],
            "p_entering_bh": p,
            "bh_q": q,
            "residual": c["residual"],
            "residual_ci": c["residual_ci"],
            "marginal": c["marginal"],
            "r2_on_axes": c["r2_on_axes"],
            "verdict": verdict_class(c, q),
            "verdict_in_A14b": c["verdict"],
        }

    adds = [s for s, v in out_cells.items() if v["verdict"] == "ADDS-VALUE"]

    result = {
        "cohort": src["cohort"],
        "family": (
            "A14 — its own BH family, per addendum section 7. m = 10: the ten "
            "signatures of section 4.2 on one cohort. The six-drug CCLE arm is "
            "confined by section 4.2 to the three discovery cohorts and cannot "
            "extend here, so m is fixed by the locked lists, not chosen."
        ),
        "m": m,
        "cell_list": SIGNATURES_4_2,
        "q_threshold": Q_THRESHOLD,
        "supersedes": {
            "file": "results/A14/A14b_analysis.json",
            "field": "family",
            "stale_value": src["family"],
            "reason": (
                "Contradicts addendum section 7, which gives A14 its own BH family "
                "and applies section 10 conditions 1 and 6 to its cells. A14b and its "
                "producing script are hash-pinned in the LOCK deposit_manifest and are "
                "left unmodified; this file supersedes the field."
            ),
        },
        "deviation": (
            "Section 7 required the cell list and integer m to be written to "
            "decisions_addendum.md before the first p value was computed. No entry "
            "does so. The deviation is logged as entry 11. It does not confer "
            "latitude over m, which section 4.2's locked lists determine."
        ),
        "cells": out_cells,
        "adds_value_cells": adds,
        "condition_1": {
            "fires": bool(adds),
            "trigger": "any cell classed ADDS-VALUE: BH q < 0.10 and residual-AUROC 95% CI lower bound > 0.5",
            "consequence": (
                "The abstract names that signature and cohort explicitly and drops the "
                "unqualified \"no signature added detectable value\"."
            ),
        },
        "ci_provenance": (
            "Cluster-bootstrap percentile intervals carried through from "
            "A14b_gse41998_analysis.py (B = 2000, seed 42, section 5.2), unmodified."
        ),
    }

    OUT.write_text(json.dumps(result, indent=2) + "\n")

    width = max(len(s) for s in SIGNATURES_4_2)
    print(f"A14 family: m = {m}, q threshold {Q_THRESHOLD}")
    print(f"{'signature':{width}}  {'LRT p':>8}  {'BH q':>8}  residual [95% CI]        verdict")
    for s in sorted(SIGNATURES_4_2, key=lambda x: out_cells[x]["bh_q"]):
        v = out_cells[s]
        lo, hi = v["residual_ci"]
        print(f"{s:{width}}  {v['lr_p']:8.4f}  {v['bh_q']:8.4f}  "
              f"{v['residual']:.3f} [{lo:.3f}, {hi:.3f}]  {v['verdict']}")
    print()
    print(f"ADDS-VALUE: {', '.join(adds) if adds else 'none'}")
    print(f"section 10 condition 1 fires: {bool(adds)}")
    print(f"written: {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
