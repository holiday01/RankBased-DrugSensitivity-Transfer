"""Phase 1 §5.2.1 — GATE diagnostic worked example #1 (Doxorubicin).

For each cohort (GSE16446, GSE25066, GSE22226):
  - Signature: frozen Doxorubicin signature score
  - Outcome: pCR
  - Stratum: PAM50 5-class (GSE25066, GSE22226); single-stratum fallback (GSE16446)

Reports marginal AUROC, within-stratum AUROC per stratum,
GAP = marginal - weighted within-mean, GAP CI, GATE p-value, verdict.

Output: results/phase_1/gate_worked_example_1_dox.json
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(".")
OUT_DIR = REPO / "." / "results" / "phase_1"
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]

sys.path.insert(0, str(REPO / "." / "scripts"))
from gate_diagnostic import gate, gate_from_dataframe  # noqa: E402


def sha256_of(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gap_ci(score, outcome, stratum, n_boot=1000, seed=42):
    """Paired bootstrap CI for GAP."""
    rng = np.random.default_rng(seed)
    s = np.asarray(score, dtype=float)
    y = np.asarray(outcome, dtype=int)
    st = np.asarray(stratum)
    n = len(y)
    gaps = []
    from sklearn.metrics import roc_auc_score
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb, sb, stb = y[idx], s[idx], st[idx]
        if len(np.unique(yb)) < 2:
            continue
        try:
            marg = roc_auc_score(yb, sb)
            ws, tw = 0.0, 0
            for lab in np.unique(stb):
                ii = stb == lab
                if ii.sum() < 5 or len(np.unique(yb[ii])) < 2:
                    continue
                ws += roc_auc_score(yb[ii], sb[ii]) * ii.sum()
                tw += ii.sum()
            if tw == 0:
                continue
            gaps.append(marg - ws / tw)
        except Exception:
            continue
    if len(gaps) < 10:
        return float("nan"), float("nan")
    return float(np.percentile(gaps, 2.5)), float(np.percentile(gaps, 97.5))


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {"script": str(me), "script_sha256": sha256_of(me),
           "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
           "cohorts": {}}

    for cohort in COHORTS:
        dox = pd.read_csv(REPO / "external_data" / cohort / "score_DOXORUBICIN.tsv", sep="\t")
        pcr = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_pcr.tsv", sep="\t")
        pam50 = pd.read_csv(REPO / "external_data" / cohort / f"{cohort}_pam50.tsv", sep="\t")
        common = sorted(set(dox["patient_id"].astype(str)) & set(pcr["patient_id"].astype(str))
                        & set(pam50["patient_id"].astype(str)))
        df = pd.DataFrame({"patient_id": common})
        df["signature_score"] = dox.set_index("patient_id").loc[common]["signature_score"].astype(float).values
        df["pCR"] = pcr.set_index("patient_id").loc[common]["pcr"].astype(int).values
        df["PAM50"] = pam50.set_index("patient_id").loc[common]["pam50_subtype"].astype(str).values

        # GSE16446 single-stratum fallback per pre-reg §6.2
        if cohort == "GSE16446" or df["PAM50"].nunique() < 2:
            res = gate(df["signature_score"].values, df["pCR"].values, None,
                       n_boot=1000, seed=42)
            ci_lo, ci_hi = float("nan"), float("nan")
        else:
            res = gate(df["signature_score"].values, df["pCR"].values,
                       df["PAM50"].values, n_boot=1000, seed=42)
            ci_lo, ci_hi = gap_ci(df["signature_score"].values,
                                  df["pCR"].values, df["PAM50"].values, 1000, 42)
        log["cohorts"][cohort] = {
            "n_total": int(len(df)),
            "n_pcr_positive": int(df["pCR"].sum()),
            "marginal_AUROC": res["marginal_AUROC"],
            "marginal_CI": list(res["marginal_CI"]),
            "within_stratum_AUROC": res["within_stratum_AUROC"],
            "stratum_sizes": res["stratum_sizes"],
            "GAP": res["GAP"],
            "GAP_CI_95": [ci_lo, ci_hi],
            "gate_p_value": res["p_value"],
            "mode": res["mode"],
            "verdict": res["verdict"],
        }

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(OUT_DIR / "gate_worked_example_1_dox.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    print("A7 done")
    for c in COHORTS:
        b = log["cohorts"][c]
        print(f"  {c}: marg={b['marginal_AUROC']:.3f} GAP={b['GAP']:.3f} "
              f"verdict={b['verdict']} mode={b['mode']}")


if __name__ == "__main__":
    main()
