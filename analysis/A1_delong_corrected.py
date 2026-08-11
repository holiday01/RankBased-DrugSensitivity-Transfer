#!/usr/bin/env python3
"""A1 -- corrected paired-DeLong / LRT-power analysis.

Accompanies the analysis addendum deposited at 10.5281/zenodo.20726403.

Corrects two defects in resubmission_v2/scripts/lrt_power_and_delong.py:

  (1) Call-site argument swap (original line 67).
      delong_var(y, p1, p2) returns (AUC(p1), AUC(p2), ...), but was called as
      delong_var(y, pf, pr) -- full predictions in the p1 slot -- and unpacked
      as (a_red, a_full). The deposited JSON therefore carries the full-model
      AUROC under the key "auc_reduced" and vice versa, and every "delta" is
      sign-flipped.

  (2) Inconsistent degenerate-path return order (original line 48).
      The var<=0 branch returned (auc2, auc1, ...) while the normal branch
      returned (auc1, auc2, ...). Never triggered on these four cohorts, but it
      is the same class of defect and would silently swap labels if it fired.

Guard against recurrence: each AUROC is independently recomputed with
sklearn.roc_auc_score from its own model's predictions and asserted equal to
the DeLong structural estimate. A future argument swap now fails loudly.

Everything else -- cohorts, axes, seed, grid, NSIM -- is unchanged from the
original, so the LRT power curve must reproduce the deposited values exactly.
That reproduction is asserted, and is the evidence that this rerun changes
only what the fix was meant to change.
"""
import contextlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from scipy import stats
from sklearn.metrics import roc_auc_score

import os

# Paths are deposit-relative by default so this runs outside the authors' machine.
# DATA_ROOT: frozen_signatures/, external_data/ and the phase_* result tree.
# OUT_ROOT : where A0/ and A1/ outputs are written.
_DR = Path(os.environ.get("DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
_OR = Path(os.environ.get("OUT_ROOT", Path(__file__).resolve().parents[1] / "results"))

P9 = _DR / "resubmission_v2/results/phase_9_revalidation"
_FROZEN = os.environ.get("FROZEN_ROOT", "")
FROZEN_P9 = (Path(_FROZEN) / "03_zenodo_v3/results_summaries/phase_9_revalidation"
             if _FROZEN else None)
OUT = _OR / "A1"

# revalidation.py has no __main__ guard: importing it re-executes the whole pipeline
# and REWRITES these deposited artefacts in the working copy. They are deterministic
# under the current seed and inputs, but a future change would clobber them silently.
# Hash them before and after so an overwrite fails loudly instead.
GUARDED = ["revalidation_summary.json", "meta_cells.tsv",
           "forest_residual_auroc.png", "power_curve.png"]


def _hash(p):
    import hashlib
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


_before = {f: _hash(P9 / f) for f in GUARDED}

sys.path.insert(0, str(P9))
with contextlib.redirect_stdout(io.StringIO()):
    from revalidation import load_cohort, AXES, zscore  # noqa: E402

_changed = [f for f in GUARDED if _hash(P9 / f) != _before[f]]
if _changed:
    raise SystemExit(
        "ABORT: importing revalidation.py altered deposited artefact(s): "
        + ", ".join(_changed)
        + "\nThese are listed in SHA256SUMS.txt. Restore them from "
        + str(FROZEN_P9) + " and investigate before proceeding.")

# Independently confirm the working copies still agree with the frozen deposit.
_drift = [f for f in GUARDED
          if (FROZEN_P9 / f).exists() and _hash(P9 / f) != _hash(FROZEN_P9 / f)]

PCR = ["GSE16446", "GSE25066", "GSE22226", "GSE32646"]
GRID = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75]
NSIM = 1000
rng = np.random.default_rng(42)


def delong_var(y, p_ref, p_new):
    """Paired DeLong for AUC(p_new) - AUC(p_ref).

    Returns (auc_ref, auc_new, z, p, se). Argument and return order agree: the
    reference model is always first in both.

    CAVEAT -- read before using p. Both AUROCs here come from models fitted on
    the SAME data, and the models are nested. DeLong's test is not valid in that
    setting (Demler, Pencina & D'Agostino, Stat Med 2012): the statistic is
    degenerate under the null. Measured on these cohorts by permutation, its
    empirical type-I error at alpha=0.05 is ~0.004 (GSE16446) and ~0.000
    (GSE25066), and its power against a planted axis-orthogonal effect at
    AUROC 0.65 is 0.142 versus 0.502 for the nested LRT. Non-significance is
    therefore close to guaranteed by construction, and these p values must be
    reported as descriptive only. The nested-model LRT is the valid test.
    """
    y = np.asarray(y)
    m = int((y == 1).sum())
    n = int((y == 0).sum())

    def structural(px):
        pos, neg = px[y == 1], px[y == 0]
        V10 = np.array([(np.sum(neg < xp) + 0.5 * np.sum(neg == xp)) / n for xp in pos])
        V01 = np.array([(np.sum(pos > xn) + 0.5 * np.sum(pos == xn)) / m for xn in neg])
        return V10, V01

    V10_ref, V01_ref = structural(p_ref)
    V10_new, V01_new = structural(p_new)
    auc_ref, auc_new = V10_ref.mean(), V10_new.mean()
    S10 = np.cov(np.vstack([V10_ref, V10_new]))
    S01 = np.cov(np.vstack([V01_ref, V01_new]))
    S = S10 / m + S01 / n
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    if var <= 0:
        # The original fabricated z=0, p=1.0 here -- a manufactured non-significant
        # result for a test that is undefined. Degeneracy is reported, not papered over.
        return auc_ref, auc_new, float("nan"), float("nan"), float("nan")
    se = float(np.sqrt(var))
    z = (auc_new - auc_ref) / se
    return auc_ref, auc_new, float(z), float(2 * stats.norm.sf(abs(z))), se


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    deposited = json.load(open(P9 / "lrt_power_and_delong.json"))

    out = {"paired_delong_reduced_vs_full": {}, "lrt_power_curve": {}}
    corrections = []

    for c in PCR:
        df = load_cohort(c)
        y = df["y"].values.astype(int)
        Xr = np.column_stack([zscore(df[a]) for a in AXES])
        Xf = np.column_stack([Xr, zscore(df["sig"])])

        mr = sm.Logit(y, sm.add_constant(Xr)).fit(disp=0)
        mf = sm.Logit(y, sm.add_constant(Xf)).fit(disp=0)
        pr = mr.predict(sm.add_constant(Xr))
        pf = mf.predict(sm.add_constant(Xf))

        a_red, a_full, z, p, se = delong_var(y, pr, pf)

        # Independent check: a swapped argument can no longer pass silently.
        assert abs(a_red - roc_auc_score(y, pr)) < 1e-9, f"{c}: reduced AUROC mislabelled"
        assert abs(a_full - roc_auc_score(y, pf)) < 1e-9, f"{c}: full AUROC mislabelled"

        lrt_p = float(stats.chi2.sf(2 * (mf.llf - mr.llf), 1))
        delta = float(a_full - a_red)
        lo, hi = delta - 1.96 * se, delta + 1.96 * se

        out["paired_delong_reduced_vs_full"][c] = dict(
            n=int(len(y)), events=int(y.sum()),
            auc_reduced=round(float(a_red), 4), auc_full=round(float(a_full), 4),
            delta=round(delta, 4), delta_se=round(se, 4),
            delta_ci95=[round(lo, 4), round(hi, 4)],
            z=round(z, 3), p=round(p, 4), lrt_p=round(lrt_p, 4),
            delong_valid=False,
            delong_caveat=("Nested in-sample models: DeLong is degenerate here "
                           "(Demler 2012). Empirical type-I error ~0.00-0.004; "
                           "power 0.14 vs 0.50 for the nested LRT at planted "
                           "AUROC 0.65. Descriptive only -- the LRT is the valid test."))

        dep = deposited["paired_delong_reduced_vs_full"][c]
        corrections.append(dict(
            cohort=c,
            published_reduced=dep["auc_reduced"], published_full=dep["auc_full"],
            published_delta=dep["delta"],
            corrected_reduced=round(float(a_red), 4), corrected_full=round(float(a_full), 4),
            corrected_delta=round(float(a_full - a_red), 4),
            p_unchanged=abs(dep["p"] - round(p, 4)) < 1e-4))

        # ---- LRT power curve: unchanged from the original ----
        n = len(y)
        H = Xr @ np.linalg.pinv(Xr.T @ Xr) @ Xr.T
        det = []
        for A in GRID:
            mu = np.sqrt(2) * stats.norm.ppf(min(max(A, 1e-4), 1 - 1e-4))
            hits = 0
            for _ in range(NSIM):
                noise = rng.standard_normal(n)
                planted = mu * y + noise
                resid = planted - H @ planted
                Xf2 = np.column_stack([Xr, zscore(resid)])
                try:
                    m0 = sm.Logit(y, sm.add_constant(Xr)).fit(disp=0)
                    m1 = sm.Logit(y, sm.add_constant(Xf2)).fit(disp=0)
                    if stats.chi2.sf(2 * (m1.llf - m0.llf), 1) < 0.05:
                        hits += 1
                except Exception:
                    pass
            det.append(round(hits / NSIM, 3))
        mde = next((GRID[i] for i, d in enumerate(det) if d >= 0.8), ">0.75")
        out["lrt_power_curve"][c] = dict(n=n, grid=GRID, detection=det, mde80_lrt=mde)

        # The fix must not perturb the power curve.
        assert det == deposited["lrt_power_curve"][c]["detection"], \
            f"{c}: power curve changed -- the rerun altered more than the labelling"

    (OUT / "lrt_power_and_delong_CORRECTED.json").write_text(json.dumps(out, indent=2))
    (OUT / "A1_correction_table.json").write_text(json.dumps(corrections, indent=2))

    print("A1 -- paired DeLong, reduced (4 axes) vs full (+60-gene signature)\n")
    if _drift:
        print(f"WARNING: working copy differs from frozen deposit for: {', '.join(_drift)}\n")
    print(f"{'cohort':10s} {'n':>4s} {'ev':>3s} | {'published':>19s} | {'corrected':>19s} | "
          f"{'dAUROC 95% CI':>18s} | {'LRT p':>7s}")
    print(f"{'':10s} {'':>4s} {'':>3s} | {'red':>6s} {'full':>6s} {'d':>5s} | "
          f"{'red':>6s} {'full':>6s} {'d':>5s} | {'(descriptive only)':>18s} |")
    print("-" * 104)
    for r, cor in zip(out["paired_delong_reduced_vs_full"].items(), corrections):
        c, v = r
        ci = f"[{v['delta_ci95'][0]:+.4f},{v['delta_ci95'][1]:+.4f}]"
        print(f"{c:10s} {v['n']:4d} {v['events']:3d} | {cor['published_reduced']:6.4f} "
              f"{cor['published_full']:6.4f} {cor['published_delta']:+.3f} | "
              f"{cor['corrected_reduced']:6.4f} {cor['corrected_full']:6.4f} "
              f"{cor['corrected_delta']:+.3f} | {ci:>18s} | {v['lrt_p']:7.4f}")

    d = [v["delta"] for v in out["paired_delong_reduced_vs_full"].values()]
    ps = [v["p"] for v in out["paired_delong_reduced_vs_full"].values()]
    lps = [v["lrt_p"] for v in out["paired_delong_reduced_vs_full"].values()]
    up = [c for c, v in out["paired_delong_reduced_vs_full"].items() if v["delta"] > 0]
    print(f"\ncorrected dAUROC range : {min(d):+.4f} to {max(d):+.4f}"
          f"   (manuscript states -0.015 to 0.000)")
    print(f"AUROC rises in         : {', '.join(up) if up else 'no cohort'}"
          f"   (manuscript states 'does not increase in any cohort')")
    print(f"DeLong p               : all >= {min(ps):.3f}  -- unchanged by the fix, but see below")
    print(f"nested-LRT p           : all >= {min(lps):.3f}  -- this is the valid test")
    print(f"power curve            : reproduces deposit exactly (asserted)")
    print(f"deposited artefacts    : unaltered by import (hash-guarded)")
    print("\nNOTE: the DeLong p values above are DESCRIPTIVE ONLY. Both models are fitted\n"
          "on the same data and are nested, where DeLong is degenerate (Demler 2012):\n"
          "measured type-I error ~0.00-0.004 at alpha=0.05, power 0.14 vs 0.50 for the\n"
          "nested LRT at planted AUROC 0.65. Non-significance is near-guaranteed by\n"
          "construction, so these p values cannot corroborate the null. Report the\n"
          "dAUROC point estimate with its CI, and rest the inference on the LRT.")
    print(f"\nwrote {OUT}/lrt_power_and_delong_CORRECTED.json")
    print(f"wrote {OUT}/A1_correction_table.json")


if __name__ == "__main__":
    main()
