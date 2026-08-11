#!/usr/bin/env python3
"""A0 -- bidirectional scorer validation.

Accompanies the analysis addendum deposited at 10.5281/zenodo.20726403.

Does the deployed bidirectional score do what the Methods say it does?

main.tex L127-129/L138 describe a bidirectional 60-gene signature scored by rank-based
singscore. A bidirectional singscore must RISE when the up genes are highly ranked AND the
down genes are lowly ranked, and must be near zero when both move together.

The deployed implementation (resubmission_v2/scripts/phase_6/_common.py:73-96 and 13 other
copies, plus 14 deposited copies under upload/zenodo_v3/) computes

    down_score = normalize((1 - rank[down]).mean(), n_down)     # already reversed
    raw        = up_score - down_score                          # reversed AGAIN

Subtracting a score that is already built on reversed ranks is a double negation. This
script decides the question on synthetic data where the right answer is known by
construction, so the verdict does not depend on any real-data interpretation.

Run:  python3 A0_scorer_validation.py
Exit: 0 if the deployed formula is correct, 1 if it is not.
"""
import sys

import numpy as np
import pandas as pd

N_GENES = 6000          # above the >=5000 guard the real scorers enforce
N_SET = 30              # matches the 30-up/30-down frozen signature
RNG = np.random.default_rng(0)


def normalize(mean_rank, n_set, n_total):
    """The exact normalisation used by every deployed scorer."""
    mn = (1 + n_set) / (2 * n_total)
    mx = 1 - mn
    return 2 * (mean_rank - mn) / (mx - mn) - 1


def deployed(ranks, up, down, n_total):
    """What the code does: down built on reversed ranks, then SUBTRACTED."""
    up_s = normalize(ranks[up].mean(axis=1), len(up), n_total)
    dn_s = normalize((1 - ranks[down]).mean(axis=1), len(down), n_total)
    return up_s - dn_s


def corrected(ranks, up, down, n_total):
    """Foroutan singscore TotalScore: down built on reversed ranks, then ADDED."""
    up_s = normalize(ranks[up].mean(axis=1), len(up), n_total)
    dn_s = normalize((1 - ranks[down]).mean(axis=1), len(down), n_total)
    return up_s + dn_s


def build():
    """Four samples whose correct ordering is fixed by construction.

    A: up HIGH, down LOW   -> signature-positive, correct score should be MAXIMAL
    B: up LOW,  down HIGH  -> signature-negative, correct score should be MINIMAL
    C: up HIGH, down HIGH  -> components cancel, correct score should be ~0
    D: up LOW,  down LOW   -> components cancel, correct score should be ~0
    """
    genes = [f"G{i:05d}" for i in range(N_GENES)]
    up = genes[:N_SET]
    down = genes[N_SET:2 * N_SET]
    rows, names = [], []

    def make(up_high, down_high):
        x = RNG.normal(0, 1, N_GENES)
        x[:N_SET] = 8.0 if up_high else -8.0
        x[N_SET:2 * N_SET] = 8.0 if down_high else -8.0
        return x

    for nm, uh, dh in [("A_up-hi_down-lo", True, False),
                       ("B_up-lo_down-hi", False, True),
                       ("C_both-hi", True, True),
                       ("D_both-lo", False, False)]:
        rows.append(make(uh, dh))
        names.append(nm)

    expr = pd.DataFrame(rows, index=names, columns=genes)
    return expr, up, down


def main():
    expr, up, down = build()
    n_total = expr.shape[1]
    ranks = expr.rank(axis=1, method="average", pct=True)

    dep = deployed(ranks, up, down, n_total)
    cor = corrected(ranks, up, down, n_total)

    print("A0 -- bidirectional scorer validation on synthetic data")
    print(f"{N_GENES} genes, {N_SET} up / {N_SET} down, ordering known by construction\n")
    print(f"{'sample':20s} {'deployed (up - down)':>22s} {'corrected (up + down)':>23s}")
    print("-" * 68)
    for s in expr.index:
        print(f"{s:20s} {dep[s]:22.4f} {cor[s]:23.4f}")

    A, B, C, D = "A_up-hi_down-lo", "B_up-lo_down-hi", "C_both-hi", "D_both-lo"
    failures = []

    print("\nRequirements a bidirectional singscore must satisfy:\n")

    # 1. The signature-positive sample must score highest.
    ok1 = dep[A] == max(dep)
    print(f"  [{'PASS' if ok1 else 'FAIL'}] deployed: A (up high, down low) is the maximum")
    if not ok1:
        failures.append(f"deployed ranks {dep.idxmax()} above A")

    # 2. The two cancelling samples must sit between A and B, near zero.
    ok2 = abs(dep[C]) < abs(dep[A]) and abs(dep[D]) < abs(dep[A])
    print(f"  [{'PASS' if ok2 else 'FAIL'}] deployed: C and D (components cancel) are near zero")
    if not ok2:
        failures.append("deployed gives cancelling samples an extreme score")

    # 3. Direction: A must outrank B.
    ok3 = dep[A] > dep[B]
    print(f"  [{'PASS' if ok3 else 'FAIL'}] deployed: A outranks B")
    if not ok3:
        failures.append("deployed inverts the signature direction")

    print("\nSame three requirements against the corrected formula:\n")
    print(f"  [{'PASS' if cor[A] == max(cor) else 'FAIL'}] corrected: A is the maximum")
    print(f"  [{'PASS' if abs(cor[C]) < abs(cor[A]) and abs(cor[D]) < abs(cor[A]) else 'FAIL'}] "
          f"corrected: C and D are near zero")
    print(f"  [{'PASS' if cor[A] > cor[B] else 'FAIL'}] corrected: A outranks B")

    # What the deployed formula is actually measuring.
    r_up = ranks[up].mean(axis=1)
    r_dn = ranks[down].mean(axis=1)
    undirected = r_up + r_dn
    contrast = r_up - r_dn
    print("\nWhat each formula actually measures (Pearson r on these 4 samples):")
    print(f"  deployed  vs undirected sum  (r_up + r_dn) : {np.corrcoef(dep, undirected)[0, 1]:+.8f}")
    print(f"  deployed  vs contrast        (r_up - r_dn) : {np.corrcoef(dep, contrast)[0, 1]:+.8f}")
    print(f"  corrected vs undirected sum                : {np.corrcoef(cor, undirected)[0, 1]:+.8f}")
    print(f"  corrected vs contrast                      : {np.corrcoef(cor, contrast)[0, 1]:+.8f}")

    print("\n" + "=" * 68)
    if failures:
        print("VERDICT: the deployed formula is WRONG.")
        for f in failures:
            print(f"  - {f}")
        print("\nIt scores 'both gene sets highly ranked' as the strongest signature-positive")
        print("sample, which is the defining behaviour of an UNDIRECTED gene-set score, not")
        print("of the bidirectional singscore the Methods describe.")
        print("\nFix: raw = up_score + down_score  (down_score already uses reversed ranks),")
        print("or equivalently drop the (1 - rank) reversal and keep the subtraction.")
        return 1
    print("VERDICT: the deployed formula behaves as a bidirectional singscore.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
