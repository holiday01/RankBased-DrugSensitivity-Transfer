"""Phase 0 SHARED MODULE — GATE (Generalization Artifact Test) diagnostic.

Imported by Phase 1 and Phase 3 (NOT run as a standalone script).

The GATE diagnostic decomposes a marginal AUROC into:
    marginal_AUROC : AUROC of (signature_score, outcome) ignoring strata
    within_AUROC   : per-stratum AUROC, weighted by stratum sample size
    GAP            = marginal_AUROC - weighted_mean(within_AUROC)
    p_value        : paired bootstrap p-value for GAP > 0

GAP threshold convention (locked in this study pre-reg §7.1):
    GAP >= 0.10  -> BASE-RATE ARTIFACT. Marginal signal is explained by the
                    stratum distribution differing between class 0 / class 1
                    (e.g., subtype-prevalence effect). NEGATIVE evidence for
                    a within-biology mechanism.
    GAP <  0.05  -> NO ARTIFACT. Marginal AUROC is recapitulated within
                    each stratum -> genuine within-stratum signal.
    0.05 <= GAP < 0.10 -> AMBIGUOUS. Requires sensitivity analyses
                    (e.g., alternate stratifiers, larger cohorts).

The original this study GATE diagnostic on doxorubicin GSE25066 was
0.752 with GAP ~0.18 -> NEGATIVE per this convention.

For cohorts with no valid PAM50 stratification (GSE16446, all ER−/basal),
the diagnostic auto-falls back to "marginal-only" mode and returns
within_AUROC = {single_stratum: marginal_AUROC} with weighted_mean equal to
marginal_AUROC, GAP = 0, and emits a warning that this is the GSE16446-style
fallback (see pre-reg §6.2).
"""
from __future__ import annotations

import warnings
from typing import Iterable, Mapping

import numpy as np
from sklearn.metrics import roc_auc_score


# Locked thresholds from pre-reg §7.1 — do not modify in module code
GAP_THRESHOLD_ARTIFACT = 0.10
GAP_THRESHOLD_NO_ARTIFACT = 0.05


def _bootstrap_auroc(y: np.ndarray, s: np.ndarray, n_boot: int = 1000, seed: int = 42):
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    obs = float(roc_auc_score(y, s))
    rng = np.random.default_rng(seed)
    n = len(y)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        boot.append(roc_auc_score(y[idx], s[idx]))
    if len(boot) < 10:
        return obs, float("nan"), float("nan")
    return obs, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def gate(
    signature_score: Iterable[float],
    outcome: Iterable[int],
    stratum: Iterable | None,
    n_boot: int = 1000,
    seed: int = 42,
) -> dict:
    """Compute the GATE diagnostic.

    Parameters
    ----------
    signature_score : array-like of float
    outcome         : array-like of int {0, 1}
    stratum         : array-like of categorical or None
                      If None or single-unique-value, falls back to
                      "marginal-only" mode (e.g., GSE16446 all-basal cohort).
    n_boot          : paired-bootstrap iterations (default 1000)
    seed            : RNG seed for bootstrap (default 42)

    Returns
    -------
    dict with keys:
        marginal_AUROC, marginal_CI (tuple (low, high))
        within_stratum_AUROC : dict {stratum_label: auroc, ..., '_weighted_mean': float}
        GAP                  : float = marginal_AUROC - weighted_within_mean
        p_value              : float = paired-bootstrap P(GAP_boot <= 0)
        mode                 : 'stratified' | 'marginal_only'
        n_total, stratum_sizes
        verdict              : 'BASE_RATE_ARTIFACT' | 'NO_ARTIFACT' | 'AMBIGUOUS' | 'INSUFFICIENT'
    """
    y = np.asarray(list(outcome), dtype=float)
    s = np.asarray(list(signature_score), dtype=float)
    if stratum is None:
        strata = None
    else:
        strata = np.asarray(list(stratum), dtype=object)

    # Filter NaNs in score / outcome
    finite = np.isfinite(y) & np.isfinite(s)
    if strata is not None:
        # treat 'nan'-string or None as missing stratum -> drop
        finite_strat = np.array([
            (x is not None) and (not (isinstance(x, float) and np.isnan(x))) and (str(x).lower() != "nan")
            for x in strata
        ])
        finite = finite & finite_strat
    y = y[finite]
    s = s[finite]
    if strata is not None:
        strata = strata[finite]

    n_total = int(len(y))
    if n_total < 10 or len(np.unique(y)) < 2:
        return {
            "marginal_AUROC": float("nan"),
            "marginal_CI": (float("nan"), float("nan")),
            "within_stratum_AUROC": {},
            "GAP": float("nan"),
            "p_value": float("nan"),
            "mode": "insufficient",
            "n_total": n_total,
            "stratum_sizes": {},
            "verdict": "INSUFFICIENT",
        }

    marg, ci_lo, ci_hi = _bootstrap_auroc(y, s, n_boot=n_boot, seed=seed)

    # Decide mode
    if strata is None or len(np.unique(strata)) < 2:
        warnings.warn(
            "GATE: single-stratum cohort (e.g., GSE16446 all-basal); "
            "falling back to marginal-only mode per pre-reg §6.2. "
            "GAP is forced to 0; verdict is 'INSUFFICIENT' for artifact detection."
        )
        return {
            "marginal_AUROC": marg,
            "marginal_CI": (ci_lo, ci_hi),
            "within_stratum_AUROC": {"_single_stratum": marg, "_weighted_mean": marg},
            "GAP": 0.0,
            "p_value": float("nan"),
            "mode": "marginal_only",
            "n_total": n_total,
            "stratum_sizes": {"_single_stratum": n_total},
            "verdict": "INSUFFICIENT",
        }

    # Stratified mode
    within = {}
    sizes = {}
    weighted_sum = 0.0
    total_weight = 0
    for lab in np.unique(strata):
        idx = strata == lab
        n_k = int(idx.sum())
        sizes[str(lab)] = n_k
        y_k = y[idx]
        s_k = s[idx]
        if n_k < 5 or len(np.unique(y_k)) < 2:
            within[str(lab)] = float("nan")
            continue
        a_k = float(roc_auc_score(y_k, s_k))
        within[str(lab)] = a_k
        weighted_sum += a_k * n_k
        total_weight += n_k

    weighted_mean = weighted_sum / total_weight if total_weight > 0 else float("nan")
    within["_weighted_mean"] = weighted_mean
    gap = marg - weighted_mean if np.isfinite(weighted_mean) else float("nan")

    # Paired bootstrap for GAP > 0
    rng = np.random.default_rng(seed)
    boot_gaps = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        y_b = y[idx]
        s_b = s[idx]
        st_b = strata[idx]
        if len(np.unique(y_b)) < 2:
            continue
        marg_b = roc_auc_score(y_b, s_b)
        ws, tw = 0.0, 0
        for lab in np.unique(st_b):
            ii = st_b == lab
            if ii.sum() < 5 or len(np.unique(y_b[ii])) < 2:
                continue
            a = roc_auc_score(y_b[ii], s_b[ii])
            ws += a * ii.sum()
            tw += ii.sum()
        if tw == 0:
            continue
        boot_gaps.append(marg_b - ws / tw)

    if boot_gaps:
        boot_gaps_arr = np.array(boot_gaps)
        # p-value: H0 GAP <= 0; one-sided P(GAP_boot <= 0 | obs)
        p_value = float((np.sum(boot_gaps_arr <= 0) + 1) / (len(boot_gaps_arr) + 1))
    else:
        p_value = float("nan")

    if not np.isfinite(gap):
        verdict = "INSUFFICIENT"
    elif gap >= GAP_THRESHOLD_ARTIFACT:
        verdict = "BASE_RATE_ARTIFACT"
    elif gap < GAP_THRESHOLD_NO_ARTIFACT:
        verdict = "NO_ARTIFACT"
    else:
        verdict = "AMBIGUOUS"

    return {
        "marginal_AUROC": marg,
        "marginal_CI": (ci_lo, ci_hi),
        "within_stratum_AUROC": within,
        "GAP": gap,
        "p_value": p_value,
        "mode": "stratified",
        "n_total": n_total,
        "stratum_sizes": sizes,
        "verdict": verdict,
    }


def gate_from_dataframe(df, score_col: str, outcome_col: str, stratum_col: str | None,
                       n_boot: int = 1000, seed: int = 42) -> dict:
    """Convenience wrapper: same as gate() but takes a DataFrame + col names."""
    s = df[score_col]
    y = df[outcome_col]
    st = df[stratum_col] if stratum_col is not None and stratum_col in df.columns else None
    return gate(s, y, st, n_boot=n_boot, seed=seed)


# Module sentinel — declare no main entry point
def _no_main_block():
    """gate_diagnostic.py is a SHARED MODULE; it has no __main__ block by design.
    Phase 1 and Phase 3 scripts import `gate` and `gate_from_dataframe`."""
    raise RuntimeError("gate_diagnostic.py is a shared module, not a script.")
