"""A10 re-run: required n to exclude a small residual AUROC, with both defects fixed.

Why a new file. `analysis/A_remaining.py` is listed in the LOCK's deposit_manifest, so it
is not edited here; this script is the corrected execution and writes to a new output.
`results/A10/A10_required_n.json` is left as deposited, still carrying the withdrawn
figures, so the record of what was originally computed survives.

The two defects, both named in `main.tex:1217-1225` and both reproduced here before being
fixed:

1. **A logit-to-probit constant applied inside a logistic model.**
   `A_remaining.py:246` returned `sqrt(2) * Phi^-1(target) * 1.7`. The relation
   AUROC = Phi(delta / sqrt(2)) already gives the effect on the standard-normal latent
   scale, so `sqrt(2) * Phi^-1(target)` is the whole answer; the extra 1.7 inflates the
   planted signal by about 70 %. A larger planted effect is easier to detect, so the
   required n came out too small — the "understated by roughly a factor of two" in the
   manuscript.

2. **The intercept set the mean of the linear predictor, not the event rate.**
   `A_remaining.py:210` shifted `lin` so that its mean equalled logit(0.23). Because the
   logistic is convex on one side of its inflection and concave on the other, the mean of
   sigmoid(lin) is not sigmoid(mean(lin)) once lin has variance, and the realised rate came
   out near 0.28. Here the intercept is solved numerically so the realised rate is 0.23.

Rather than trusting the analytic relation after correcting it, the realised residual AUROC
is measured back out of the simulation and reported alongside. If the plant does not land on
its target, the required n is answering a different question and the reader can see that.
"""

import json
import pathlib

import numpy as np
import statsmodels.api as sm
from scipy import optimize, stats

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "A10" / "A10b_required_n_corrected.json"

SEED = 42
EVENT_RATE = 0.23
REPS = 400          # matches the deposited run, so only the defects differ
ALPHA = 0.05
POWER = 0.80
TARGETS = [0.55, 0.60]
N_LO, N_HI = 100, 12000   # upper bound raised: the corrected requirement is larger


def beta_analytic(target):
    """Effect on the standard-normal latent scale. No logit-to-probit rescaling.

    AUROC = Phi(delta / sqrt(2)) for a single standard-normal predictor, so
    delta = sqrt(2) * Phi^-1(AUROC). Defect 1 was multiplying this by 1.7.
    """
    return float(np.sqrt(2) * stats.norm.ppf(target))


def realised_auroc_at(beta, target_seed, n=4000, reps=60):
    """Mean residual AUROC actually produced by planting `beta`."""
    rng = np.random.default_rng([SEED, target_seed, 7777])
    vals = []
    for _ in range(reps):
        X = rng.normal(size=(n, 4))
        lin = 0.9 * X[:, 0] - 0.6 * X[:, 1] + 0.5 * X[:, 3]
        s = rng.normal(size=n)
        lin = lin + beta * s
        lin = lin - lin.mean()
        lin = lin + intercept_for(lin, EVENT_RATE)
        y = (rng.uniform(size=n) < 1.0 / (1.0 + np.exp(-lin))).astype(int)
        a = residual_auroc(s, X, y)
        if not np.isnan(a):
            vals.append(a)
    return float(np.mean(vals))


def calibrate_beta(target, target_seed):
    """Solve for the beta whose REALISED residual AUROC is the target.

    The analytic value is only a starting point. The quantity the paper reports is the
    residual AUROC after axis adjustment on a binary outcome at a 23 % event rate, and
    that is not exactly the latent-scale relation — the uncorrected run's own diagnostic
    showed the plant landing at 0.540 when 0.55 was intended. Calibrating removes the
    residual gap instead of carrying it into the answer.
    """
    lo, hi = 0.05, 4.0
    for _ in range(24):
        mid = (lo + hi) / 2
        if realised_auroc_at(mid, target_seed) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-3:
            break
    return (lo + hi) / 2


def intercept_for(lin_centred, rate):
    """Solve for a with mean(sigmoid(a + lin)) == rate. Defect 2 was skipping this."""
    def gap(a):
        return float(np.mean(1.0 / (1.0 + np.exp(-(a + lin_centred)))) - rate)
    return float(optimize.brentq(gap, -20.0, 20.0))


def residual_auroc(s, X, y):
    """The paper's own instrument: AUROC of the signature after OLS-adjusting for axes."""
    Xc = sm.add_constant(X)
    resid = s - Xc @ np.linalg.lstsq(Xc, s, rcond=None)[0]
    pos, neg = resid[y == 1], resid[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float((stats.rankdata(np.r_[pos, neg])[:len(pos)].sum()
                  - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def sim_power(n, target_seed, beta, reps=REPS):
    """Deterministic in (n, target): the binary search below must not chase Monte Carlo
    noise, which it does if every evaluation draws from one advancing stream."""
    rng = np.random.default_rng([SEED, target_seed, n])
    hits = used = 0
    realised = []
    for _ in range(reps):
        X = rng.normal(size=(n, 4))
        lin = 0.9 * X[:, 0] - 0.6 * X[:, 1] + 0.5 * X[:, 3]
        s = rng.normal(size=n)
        lin = lin + beta * s
        lin = lin - lin.mean()
        lin = lin + intercept_for(lin, EVENT_RATE)
        y = (rng.uniform(size=n) < 1.0 / (1.0 + np.exp(-lin))).astype(int)
        if min(y.sum(), n - y.sum()) < 5:
            continue
        R = sm.add_constant(X)
        F = np.column_stack([R, s])
        try:
            mr = sm.Logit(y, R).fit(disp=0)
            mf = sm.Logit(y, F).fit(disp=0)
        except Exception:
            continue
        used += 1
        realised.append((float(y.mean()), residual_auroc(s, X, y)))
        if stats.chi2.sf(2 * (mf.llf - mr.llf), 1) < ALPHA:
            hits += 1
    if used == 0:
        return 0.0, {}
    rates = [r for r, _ in realised]
    aurocs = [a for _, a in realised if not np.isnan(a)]
    return hits / used, {
        "replicates_used": used,
        "realised_event_rate": round(float(np.mean(rates)), 4),
        "realised_residual_auroc": round(float(np.mean(aurocs)), 4),
    }


def main():
    rng = np.random.default_rng(SEED)
    res = {
        "supersedes": "results/A10/A10_required_n.json",
        "withdrawn_figures": {"0.55": 560, "0.60": 146},
        "defects_corrected": [
            "logit-to-probit constant 1.7 removed from the planted effect "
            "(A_remaining.py:246); it inflated the signal and understated required n",
            "intercept solved so the realised event rate is 0.23 rather than setting the "
            "mean of the linear predictor (A_remaining.py:210), which realised ~0.28",
        ],
        "settings": {"seed": SEED, "event_rate": EVENT_RATE, "replicates": REPS,
                     "alpha": ALPHA, "power": POWER, "test": "nested-model LRT, 1 df"},
        "targets": {},
    }

    for target in TARGETS:
        target_seed = int(round(target * 1000))
        beta = calibrate_beta(target, target_seed)
        lo, hi = N_LO, N_HI
        while hi - lo > 50:
            mid = (lo + hi) // 2
            p, _ = sim_power(mid, target_seed, beta)
            if p >= POWER:
                hi = mid
            else:
                lo = mid
        power_at, diag = sim_power(hi, target_seed, beta)
        # A single crossing point from a 400-replicate search is worth about +/- 0.02 in
        # power, so report the curve around it. A reader can then see how sharp the
        # threshold is rather than taking one integer on trust.
        curve = {}
        for n in sorted({max(N_LO, hi + d) for d in (-150, -100, -50, 0, 50, 100)}):
            p, _ = sim_power(n, target_seed, beta)
            curve[str(n)] = round(float(p), 3)
        res["targets"][str(target)] = {
            "required_n_80pct": int(hi),
            "power_at_required": round(float(power_at), 4),
            "search_granularity": 50,
            "power_curve": curve,
            "planted_beta_calibrated": round(beta, 4),
            "planted_beta_analytic": round(beta_analytic(target), 4),
            **diag,
        }
        print(f"  residual AUROC {target}: required n ~ {hi} "
              f"(power {power_at:.3f}, realised event rate "
              f"{diag.get('realised_event_rate')}, realised residual AUROC "
              f"{diag.get('realised_residual_auroc')})")

    res["largest_public_cohort"] = {
        "name": "GSE25066", "n": 488,
        "note": "No public neoadjuvant pCR cohort reaches the n required to exclude a "
                "residual AUROC of 0.55 at 80% power.",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, indent=2) + "\n")
    print(f"written: {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
