"""Track D — 2nd held-out: score BORTEZOMIB signatures on GSE68871 (VTD in
newly-diagnosed MM; Terragna 2016; GPL570; CD138+ purified plasma cells).

Myeloma response binarization is a researcher degree of freedom, so we report
ALL pre-specified thresholds transparently and flag the PRIMARY:
  PRIMARY     : >=VGPR (CR/nCR/VGPR) = responder vs (PR/SD/PD) = non  [deep response]
  sensitivity : CR/nCR vs rest ; and >=PR (mirrors GSE9782 'R') vs SD/PD
Higher signature score = more sensitive -> expected to associate with responder.
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
from score_external_cohort import compute_singscore, load_signature, auroc_with_bootstrap_ci

COHORT = REPO / "revise_bioadv" / "external_data" / "GSE68871"
PHASE_D = REPO / "revise_bioadv" / "resubmission_v2" / "results" / "phase_D"
OUT = PHASE_D / "gse68871_eval"; OUT.mkdir(parents=True, exist_ok=True)

THRESHOLDS = {
    "primary_geVGPR": {"resp": {"CR", "NCR", "VGPR"}},
    "sens_CRnCR":     {"resp": {"CR", "NCR"}},
    "sens_gePR":      {"resp": {"CR", "NCR", "VGPR", "PR"}},
}


def get_response_series(meta):
    ids = meta["sample_ids"]; lines = meta["sample_chars"]["characteristics"]
    out = {}
    for i, sid in enumerate(ids):
        for line in lines:
            v = line[i] if i < len(line) else ""
            if "response to vtd" in v.lower() and ":" in v:
                out[sid] = v.split(":", 1)[1].strip().upper()
    return pd.Series(out)


def main():
    res = {"started_utc": datetime.now(timezone.utc).isoformat(), "cohort": "GSE68871"}
    expr_raw, meta = read_series_matrix(COHORT / "series_matrix.txt.gz")
    annot = read_gpl_annot(COHORT / "platform_annot.annot.gz")
    p2s = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    collapsed, diag = collapse_probes_max_iqr(expr_raw, p2s)
    res["expr_genes"] = diag["n_unique_symbols_after_collapse"]
    if "sample_ids" not in meta:
        meta["sample_ids"] = list(expr_raw.columns)
    resp = get_response_series(meta)
    res["response_distribution"] = resp.value_counts().to_dict()
    print("response distribution:", res["response_distribution"])

    common = [s for s in resp.index if s in collapsed.index]
    expr = collapsed.loc[common].rename(columns={c: c.upper() for c in collapsed.columns})
    rcat = resp.loc[common]
    res["n_samples"] = len(common)
    rank_pct = expr.rank(axis=1, method="average", pct=True)

    # precompute signature scores once
    sigscores = {}
    for label, path in [("pan_cancer", PHASE_D / "BORTEZOMIB.tsv"),
                        ("heme_matched", PHASE_D / "BORTEZOMIB_heme.tsv")]:
        sens, resg, _, info = load_signature(path)
        sc, sd = compute_singscore(expr, sens, resg)
        sigscores[label] = (sc.loc[common].to_numpy(), sd)

    res["by_threshold"] = {}
    for tname, spec in THRESHOLDS.items():
        y = rcat.isin(spec["resp"]).astype(int).to_numpy()
        block = {"n_responder": int(y.sum()), "n_nonresponder": int((y == 0).sum()),
                 "responder_rate": float(y.mean()), "signatures": {}, "baselines": {}}
        if len(np.unique(y)) < 2:
            block["note"] = "degenerate (single class)"; res["by_threshold"][tname] = block; continue
        for label, (s, sd) in sigscores.items():
            auc, lo, hi = auroc_with_bootstrap_ci(s, y, 1000, 42)
            block["signatures"][label] = {"auroc": auc, "ci_low": lo, "ci_high": hi,
                                          "n_up_present": sd["n_up_genes_present_in_cohort"],
                                          "n_down_present": sd["n_down_genes_present_in_cohort"]}
        for g in ["MKI67", "TNFRSF17", "IRF4", "PSMB5"]:
            if g in rank_pct.columns:
                a, l, h = auroc_with_bootstrap_ci(rank_pct[g].to_numpy(), y, 1000, 42)
                block["baselines"][g] = {"auroc": a, "ci_low": l, "ci_high": h}
        res["by_threshold"][tname] = block
        print(f"\n[{tname}] responder rate {block['responder_rate']:.2f} "
              f"(n+{block['n_responder']}/n-{block['n_nonresponder']})")
        for label, d in block["signatures"].items():
            print(f"   {label:12s} AUROC {d['auroc']:.3f} [{d['ci_low']:.3f}, {d['ci_high']:.3f}]")
        for g, d in block["baselines"].items():
            print(f"   baseline {g:9s} {d['auroc']:.3f}")

    res["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT / "gse68871_eval.json").write_text(json.dumps(res, indent=2))
    print("\nSaved ->", OUT / "gse68871_eval.json")


if __name__ == "__main__":
    main()
