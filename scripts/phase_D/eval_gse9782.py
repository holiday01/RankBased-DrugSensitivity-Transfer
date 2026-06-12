"""Track D — score BORTEZOMIB signatures on GSE9782 (bortezomib MM, Mulligan 2007)
and measure response AUROC. NO inference — measured outcome decides D.

Signatures tested:
  1. BORTEZOMIB.tsv       (pan-cancer pool; framework verbatim)
  2. BORTEZOMIB_heme.tsv  (heme-lineage-matched; plasma-cell/myeloma biology)
Baselines (single-gene rank, within-sample percentile):
  PSMB5 (bortezomib's direct target), TNFRSF17/BCMA, IRF4, MKI67 (proliferation)

Label: PGx_Responder R=1 / NR=0; treatment filtered to PS341 (bortezomib arm).
Higher signature score = more 'sensitive' -> expected to associate with R.
"""
from __future__ import annotations
import sys, gzip, json, re
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

REPO = Path(".")
sys.path.insert(0, str(REPO / "." / "scripts"))
from parse_gse_to_matrix import read_series_matrix, read_gpl_annot, collapse_probes_max_iqr
from score_external_cohort import compute_singscore, load_signature, auroc_with_bootstrap_ci

COHORT = REPO / "." / "external_data" / "GSE9782"
PHASE_D = REPO / "." / "." / "results" / "phase_D"
OUT = PHASE_D / "gse9782_eval"
OUT.mkdir(parents=True, exist_ok=True)


def parse_ragged_phenotype(metadata: dict) -> pd.DataFrame:
    """GSE9782 characteristics are ragged 'key = value' (NOT positional, NOT ':')."""
    sample_ids = metadata.get("sample_ids", [])
    char_lines = metadata.get("sample_chars", {}).get("characteristics", [])
    rows = []
    for i, sid in enumerate(sample_ids):
        row = {"patient_id": sid}
        for line in char_lines:
            val = line[i] if i < len(line) else ""
            if "=" in val:
                k, v = val.split("=", 1)
                row[k.strip().lower().replace(" ", "_")] = v.strip()
        rows.append(row)
    return pd.DataFrame(rows)


def singlegene_auroc(expr_rank_pct: pd.DataFrame, gene: str, y: np.ndarray):
    g = gene.upper()
    if g not in expr_rank_pct.columns:
        return None
    s = expr_rank_pct[g].to_numpy()
    auc, lo, hi = auroc_with_bootstrap_ci(s, y, 1000, 42)
    return {"gene": g, "auroc": auc, "ci_low": lo, "ci_high": hi}


def main():
    res = {"started_utc": datetime.now(timezone.utc).isoformat(), "cohort": "GSE9782"}

    # --- expression: probe x sample -> sample x symbol ---
    expr_raw, meta = read_series_matrix(COHORT / "series_matrix.txt.gz")
    annot = read_gpl_annot(COHORT / "platform_annot.annot.gz")
    p2s = dict(zip(annot["probe_id"], annot["gene_symbol"]))
    collapsed, diag = collapse_probes_max_iqr(expr_raw, p2s)  # sample x symbol
    res["expr_diag"] = diag

    # --- phenotype ---
    if "sample_ids" not in meta:
        meta["sample_ids"] = list(expr_raw.columns)
    phen = parse_ragged_phenotype(meta)
    # locate treatment + responder cols
    treat_col = next((c for c in phen.columns if c.startswith("treatment")), None)
    resp_col = next((c for c in phen.columns if "pgx_responder" in c), None)
    respcat_col = next((c for c in phen.columns if "pgx_response" == c or c == "pgx_response"), None)
    res["phenotype_columns"] = list(phen.columns)
    res["treat_col"], res["resp_col"] = treat_col, resp_col
    phen = phen.set_index("patient_id")

    # filter: bortezomib arm (PS341) + evaluable responder
    treat = phen[treat_col].astype(str).str.upper()
    resp = phen[resp_col].astype(str).str.upper()
    keep = phen.index[(treat.str.contains("PS341")) & (resp.isin(["R", "NR"]))]
    res["n_total_samples"] = int(phen.shape[0])
    res["n_ps341_arm"] = int((treat.str.contains("PS341")).sum())
    res["n_kept_evaluable"] = int(len(keep))

    common = [s for s in keep if s in collapsed.index]
    expr = collapsed.loc[common]
    y = (resp.loc[common] == "R").astype(int).to_numpy()
    res["n_scored"] = int(len(common))
    res["n_responders"] = int(y.sum())
    res["n_nonresponders"] = int((y == 0).sum())
    res["response_rate"] = float(y.mean())

    # per-sample percentile rank matrix (for single-gene baselines)
    expr_u = expr.rename(columns={c: c.upper() for c in expr.columns})
    rank_pct = expr_u.rank(axis=1, method="average", pct=True)

    # --- signatures ---
    sig_results = {}
    for label, path in [("pan_cancer", PHASE_D / "BORTEZOMIB.tsv"),
                        ("heme_matched", PHASE_D / "BORTEZOMIB_heme.tsv")]:
        sens, resg, _, info = load_signature(path)
        try:
            scores, sdiag = compute_singscore(expr_u, sens, resg)
            s = scores.loc[common].to_numpy()
            auc, lo, hi = auroc_with_bootstrap_ci(s, y, 1000, 42)
            sig_results[label] = {"auroc": auc, "ci_low": lo, "ci_high": hi,
                                  "n_up_present": sdiag["n_up_genes_present_in_cohort"],
                                  "n_down_present": sdiag["n_down_genes_present_in_cohort"],
                                  "signature_sha256": info["signature_sha256"]}
            # save per-patient scores
            pd.DataFrame({"patient_id": common, "score": s, "responder": y}).to_csv(
                OUT / f"scores_{label}.tsv", sep="\t", index=False, float_format="%.10g")
            print(f"[{label:12s}] AUROC {auc:.3f} [{lo:.3f}, {hi:.3f}]  "
                  f"(up {sdiag['n_up_genes_present_in_cohort']}/{len(sens)}, "
                  f"down {sdiag['n_down_genes_present_in_cohort']}/{len(resg)})")
        except Exception as e:
            sig_results[label] = {"error": str(e)}
            print(f"[{label}] ERROR {e}")
    res["signatures"] = sig_results

    # --- single-gene baselines ---
    res["baselines"] = {}
    for g in ["PSMB5", "TNFRSF17", "IRF4", "IKZF3", "MKI67", "PSMB1", "XBP1"]:
        b = singlegene_auroc(rank_pct, g, y)
        if b:
            res["baselines"][g] = b
            print(f"  baseline {g:9s} AUROC {b['auroc']:.3f} [{b['ci_low']:.3f}, {b['ci_high']:.3f}]")

    res["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT / "gse9782_eval.json").write_text(json.dumps(res, indent=2))
    print("\nSaved ->", OUT / "gse9782_eval.json")


if __name__ == "__main__":
    main()
