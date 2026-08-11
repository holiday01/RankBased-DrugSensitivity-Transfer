"""Track D robustness — GATE within-stratum, permutation null, and paired
bootstrap vs MKI67 proliferation baseline for the heme-matched bortezomib
signature on GSE9782. Decides whether the AUROC 0.66 is REAL within-stratum
signal (unlike the original 0.752 base-rate artifact).
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path("/home/holiday01/2026_ISMB_code")
sys.path.insert(0, str(REPO / "revise_bioadv" / "scripts"))
from parse_gse_to_matrix import read_series_matrix, read_gpl_annot, collapse_probes_max_iqr
from score_external_cohort import compute_singscore, load_signature

COHORT = REPO / "revise_bioadv" / "external_data" / "GSE9782"
PHASE_D = REPO / "revise_bioadv" / "resubmission_v2" / "results" / "phase_D"
OUT = PHASE_D / "gse9782_eval"


def parse_ragged(meta):
    ids = meta["sample_ids"]; lines = meta["sample_chars"]["characteristics"]
    rows = []
    for i, sid in enumerate(ids):
        row = {"patient_id": sid}
        for line in lines:
            v = line[i] if i < len(line) else ""
            if "=" in v:
                k, val = v.split("=", 1); row[k.strip().lower().replace(" ", "_")] = val.strip()
        rows.append(row)
    return pd.DataFrame(rows).set_index("patient_id")


def pooled_within_stratum_auroc(score, y, strata):
    """Average within-stratum AUROC weighted by stratum n (strata with both classes)."""
    aucs, ws = [], []
    for s in pd.unique(strata):
        idx = np.where(strata == s)[0]
        if len(idx) >= 8 and len(np.unique(y[idx])) == 2:
            aucs.append(roc_auc_score(y[idx], score[idx])); ws.append(len(idx))
    if not aucs:
        return float("nan"), {}
    aucs = np.array(aucs); ws = np.array(ws)
    return float(np.average(aucs, weights=ws)), {"per_stratum": list(map(float, aucs)), "n": list(map(int, ws))}


def main():
    res = {"started_utc": datetime.now(timezone.utc).isoformat()}
    expr_raw, meta = read_series_matrix(COHORT / "series_matrix.txt.gz")
    annot = read_gpl_annot(COHORT / "platform_annot.annot.gz")
    p2s = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    collapsed, _ = collapse_probes_max_iqr(expr_raw, p2s)
    if "sample_ids" not in meta:
        meta["sample_ids"] = list(expr_raw.columns)
    phen = parse_ragged(meta)

    treat = phen[[c for c in phen.columns if c.startswith("treatment")][0]].astype(str).str.upper()
    resp = phen[[c for c in phen.columns if "pgx_responder" in c][0]].astype(str).str.upper()
    tc_col = next((c for c in phen.columns if "tc_class" in c), None)
    keep = [s for s in phen.index if "PS341" in treat.get(s, "") and resp.get(s) in ("R", "NR")
            and s in collapsed.index]
    expr = collapsed.loc[keep].rename(columns={c: c.upper() for c in collapsed.columns})
    y = (resp.loc[keep] == "R").astype(int).to_numpy()
    rank_pct = expr.rank(axis=1, method="average", pct=True)
    res["n"] = len(keep); res["response_rate"] = float(y.mean())

    # heme-matched signature score
    sens, resg, _, info = load_signature(PHASE_D / "BORTEZOMIB_heme.tsv")
    score = compute_singscore(expr, sens, resg)[0].loc[keep].to_numpy()
    marginal = roc_auc_score(y, score)
    res["heme_marginal_auroc"] = float(marginal)

    # --- GATE: within molecular-subtype (TC_Class) stratum ---
    if tc_col:
        strata = phen.loc[keep, tc_col].astype(str).str.strip().to_numpy()
        within, wd = pooled_within_stratum_auroc(score, y, strata)
        res["gate"] = {"stratifier": tc_col, "marginal_auroc": float(marginal),
                       "within_stratum_auroc": within, "GAP": float(marginal - within),
                       "detail": wd,
                       "n_strata_used": len(wd.get("per_stratum", []))}

    # --- permutation null (marginal label shuffle; valid: 50% base rate) ---
    rng = np.random.default_rng(1337)
    perm = [roc_auc_score(rng.permutation(y), score) for _ in range(2000)]
    res["perm_p_marginal"] = float((np.sum(np.array(perm) >= marginal) + 1) / (len(perm) + 1))

    # --- paired bootstrap: signature vs MKI67 proliferation baseline ---
    mki67 = rank_pct["MKI67"].to_numpy() if "MKI67" in rank_pct.columns else None
    if mki67 is not None:
        auc_mki = roc_auc_score(y, mki67)
        rng2 = np.random.default_rng(42); diffs = []
        n = len(y)
        for _ in range(2000):
            idx = rng2.integers(0, n, n)
            if len(np.unique(y[idx])) < 2:
                continue
            diffs.append(roc_auc_score(y[idx], score[idx]) - roc_auc_score(y[idx], mki67[idx]))
        diffs = np.array(diffs)
        res["vs_mki67"] = {
            "sig_auroc": float(marginal), "mki67_auroc": float(auc_mki),
            "delta": float(marginal - auc_mki),
            "delta_ci": [float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))],
            "paired_boot_p_two_sided": float(2 * min((diffs <= 0).mean(), (diffs >= 0).mean())),
        }

    res["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT / "gse9782_robustness.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: res[k] for k in res if k not in ("started_utc", "finished_utc")}, indent=2))


if __name__ == "__main__":
    main()
