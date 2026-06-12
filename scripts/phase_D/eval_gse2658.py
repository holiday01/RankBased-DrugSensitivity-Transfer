"""Track D replication — GSE2658 (UAMS): bortezomib signatures vs survival.
Self-contained controlled design:
  TT3 (n=208, bortezomib-containing) = bortezomib test arm
  TT2 (n=351, NO bortezomib)         = control arm (prognostic-vs-predictive)
Survival: SURTIM (months) + SURIND (disease-related death 0/1), ragged-embedded.
Endpoint: Cox HR per +1 SD signature score (HR<1 = higher score -> better survival,
expected if signature marks bortezomib sensitivity). C-index. Plus signature x
protocol interaction. Baseline MKI67. All reported regardless of outcome.
"""
from __future__ import annotations
import sys, gzip, json, re
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from sklearn.metrics import roc_auc_score

REPO = Path(".")
sys.path.insert(0, str(REPO / "." / "scripts"))
from parse_gse_to_matrix import read_series_matrix, read_gpl_annot, collapse_probes_max_iqr
from score_external_cohort import compute_singscore, load_signature, auroc_with_bootstrap_ci

COHORT = REPO / "." / "external_data" / "GSE2658"
PHASE_D = REPO / "." / "." / "results" / "phase_D"
OUT = PHASE_D / "gse2658_eval"; OUT.mkdir(parents=True, exist_ok=True)


def parse_survival_and_protocol(meta):
    ids = meta["sample_ids"]; titles = meta.get("sample_titles", [""] * len(ids))
    char_lines = meta["sample_chars"]["characteristics"]
    rows = []
    for i, sid in enumerate(ids):
        surtim = surind = None
        for line in char_lines:
            v = line[i] if i < len(line) else ""
            mt = re.search(r"SURTIM=([0-9.]+)", v)
            mi = re.search(r"SURIND=(\d)", v)
            if mt:
                surtim = float(mt.group(1))
            if mi:
                surind = int(mi.group(1))
        prot = re.search(r"(TT\d)", titles[i] if i < len(titles) else "")
        rows.append({"patient_id": sid, "surtim": surtim, "surind": surind,
                     "protocol": prot.group(1) if prot else "NA"})
    return pd.DataFrame(rows).set_index("patient_id")


def cox_hr_per_sd(score, time, event):
    """Univariate Cox; score standardized to per-SD. Returns HR, CI, p, C-index."""
    df = pd.DataFrame({"t": time, "e": event, "z": (score - np.mean(score)) / np.std(score)})
    df = df.dropna()
    df = df[df["t"] > 0]
    if df["e"].sum() < 5 or len(df) < 20:
        return None
    cph = CoxPHFitter()
    cph.fit(df, duration_col="t", event_col="e")
    hr = float(np.exp(cph.params_["z"]))
    lo = float(np.exp(cph.confidence_intervals_.iloc[0, 0]))
    hi = float(np.exp(cph.confidence_intervals_.iloc[0, 1]))
    p = float(cph.summary.loc["z", "p"])
    c = float(cph.concordance_index_)
    return {"HR_per_SD": hr, "HR_ci": [lo, hi], "p": p, "c_index": c,
            "n": int(len(df)), "events": int(df["e"].sum())}


def main():
    res = {"started_utc": datetime.now(timezone.utc).isoformat(), "cohort": "GSE2658"}
    expr_raw, meta = read_series_matrix(COHORT / "series_matrix.txt.gz")
    if "sample_ids" not in meta:
        meta["sample_ids"] = list(expr_raw.columns)
    annot = read_gpl_annot(COHORT / "platform_annot.annot.gz")
    p2s = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    collapsed, diag = collapse_probes_max_iqr(expr_raw, p2s)
    res["expr_genes"] = diag["n_unique_symbols_after_collapse"]

    clin = parse_survival_and_protocol(meta)
    res["protocol_counts"] = clin["protocol"].value_counts().to_dict()
    res["n_with_survival"] = int(clin[["surtim", "surind"]].dropna().shape[0])

    expr = collapsed.rename(columns={c: c.upper() for c in collapsed.columns})
    common = [s for s in clin.index if s in expr.index]
    expr = expr.loc[common]; clin = clin.loc[common]
    rank_pct = expr.rank(axis=1, method="average", pct=True)

    # signature scores (whole-cohort singscore; deployment is per-sample rank)
    scores = {}
    for label, path in [("pan_cancer", PHASE_D / "BORTEZOMIB.tsv"),
                        ("heme_matched", PHASE_D / "BORTEZOMIB_heme.tsv")]:
        sens, resg, _, info = load_signature(path)
        sc, _ = compute_singscore(expr, sens, resg)
        scores[label] = sc.loc[common].to_numpy()

    res["by_arm"] = {}
    for arm in ["TT3", "TT2"]:  # TT3 = bortezomib, TT2 = control
        mask = (clin["protocol"] == arm).to_numpy()
        sub = clin[mask]
        block = {"n": int(mask.sum()), "events": int(sub["surind"].sum(skipna=True)),
                 "signatures": {}, "baselines": {}}
        t = sub["surtim"].to_numpy(); e = sub["surind"].to_numpy()
        for label in scores:
            cox = cox_hr_per_sd(scores[label][mask], t, e)
            block["signatures"][label] = cox
        # baselines
        for g in ["MKI67", "TNFRSF17"]:
            if g in rank_pct.columns:
                block["baselines"][g] = cox_hr_per_sd(rank_pct[g].to_numpy()[mask], t, e)
        res["by_arm"][arm] = block
        print(f"\n=== {arm} (n={block['n']}, events={block['events']}) ===")
        for label, c in block["signatures"].items():
            if c:
                print(f"  {label:12s} HR/SD {c['HR_per_SD']:.3f} [{c['HR_ci'][0]:.3f},{c['HR_ci'][1]:.3f}] p={c['p']:.3g} C={c['c_index']:.3f}")
        for g, c in block["baselines"].items():
            if c:
                print(f"  baseline {g:9s} HR/SD {c['HR_per_SD']:.3f} p={c['p']:.3g} C={c['c_index']:.3f}")

    res["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT / "gse2658_eval.json").write_text(json.dumps(res, indent=2))
    print("\nSaved ->", OUT / "gse2658_eval.json")
    print("NOTE direction: HR/SD < 1 = higher signature score -> better survival "
          "(expected if signature marks bortezomib sensitivity).")


if __name__ == "__main__":
    main()
