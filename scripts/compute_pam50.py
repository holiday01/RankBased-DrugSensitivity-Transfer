"""Compute PAM50 basal score per patient.

Per pre-reg §3.1:
- Centroid source: genefu pam50$centroids (Parker 2009; exported to
  external_data/pam50_centroids.tsv with sha256 recorded).
- Algorithm: Spearman correlation of each patient's cohort-internal-rank
  expression with each of 5 PAM50 centroids; assign patient to highest-r class.
- Output: <COHORT>_pam50.tsv with columns
    patient_id, pam50_subtype, basal_score, her2_score, luma_score, lumb_score, normal_score
- If <40 of 50 PAM50 genes map to cohort symbols, log fallback and write
  pam50_subtype=NA + basal_score=NaN (caller must use ER-/HER2- indicator).
"""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--expression", required=True, type=Path,
                    help="TSV: rows=patient_id, cols=HGNC_symbol (cohort-internal expression).")
    ap.add_argument("--centroids", required=True, type=Path,
                    help="PAM50 centroid TSV: rows=gene, cols=Basal,Her2,LumA,LumB,Normal")
    ap.add_argument("--alias-table", default=None, type=Path,
                    help="HGNC alias TSV (optional, with symbol/alias_symbol columns).")
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--out-log", required=True, type=Path)
    args = ap.parse_args()

    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]

    cent = pd.read_csv(args.centroids, sep="\t", index_col=0)
    cent.index = [g.upper() for g in cent.index]
    # genefu standard column ordering
    cent = cent[["Basal", "Her2", "LumA", "LumB", "Normal"]]

    # Alias resolution
    alias_map = {}
    if args.alias_table and args.alias_table.exists():
        adf = pd.read_csv(args.alias_table, sep="\t")
        # Build a map from any prev/alias symbol -> approved symbol
        for _, row in adf.iterrows():
            approved = str(row.get("symbol", "")).upper().strip()
            if not approved:
                continue
            for col in ["alias_symbol", "prev_symbol"]:
                v = row.get(col)
                if pd.isna(v) or not v:
                    continue
                for tok in str(v).split("|"):
                    tok = tok.strip().upper()
                    if tok and tok != approved:
                        alias_map.setdefault(tok, approved)

    cohort_genes = set(expr.columns)
    pam50_genes = list(cent.index)
    missing = [g for g in pam50_genes if g not in cohort_genes]
    # Try alias resolution for missing
    resolved = {}
    for g in missing:
        if g in alias_map and alias_map[g] in cohort_genes:
            resolved[g] = alias_map[g]
        else:
            # Reverse lookup: if cohort has an alias that maps to PAM50 gene
            for alias, approved in alias_map.items():
                if approved == g and alias in cohort_genes:
                    resolved[g] = alias
                    break

    available_pam50 = []
    rename_pam50 = {}
    for g in pam50_genes:
        if g in cohort_genes:
            available_pam50.append(g)
        elif g in resolved:
            available_pam50.append(resolved[g])
            rename_pam50[resolved[g]] = g

    n_mapped = len(available_pam50)
    fallback = n_mapped < 40

    log = {
        "cohort": args.cohort,
        "n_pam50_genes_in_centroid": len(pam50_genes),
        "n_pam50_genes_in_cohort": n_mapped,
        "missing_pam50_genes": sorted(set(missing) - set(resolved.keys())),
        "resolved_via_alias": resolved,
        "fallback_triggered": fallback,
        "centroid_sha256_source": "external_data/pam50_centroids.tsv (genefu pam50 v2.42.0)",
    }

    if fallback:
        # Pre-reg §3.1 fallback: cannot compute PAM50; write NaN
        out = pd.DataFrame({"patient_id": expr.index})
        out["pam50_subtype"] = "NA_FALLBACK"
        for s in ["Basal", "Her2", "LumA", "LumB", "Normal"]:
            out[f"{s.lower()}_score"] = np.nan
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.out_tsv, sep="\t", index=False)
        with open(args.out_log, "w") as f:
            json.dump(log, f, indent=2)
        print(f"[{args.cohort}] FALLBACK: only {n_mapped}/50 PAM50 genes available; using ER-/HER2- indicator per pre-reg §3.1")
        return

    # Build expr_cohort restricted to PAM50 genes (cohort symbol space)
    expr_sub = expr[available_pam50].copy()
    # Re-label cohort symbols back to PAM50 canonical names (so centroid alignment works)
    expr_sub.columns = [rename_pam50.get(c, c) for c in expr_sub.columns]

    # Cohort-internal rank: within each gene across patients (genefu standard preprocessing)
    # genefu intrinsic.cluster.predict default: median-centered per gene + Spearman with centroid.
    # We use the Spearman with centroid directly, gene-wise centering by median.
    expr_sub = expr_sub.subtract(expr_sub.median(axis=0), axis=1)

    cent_sub = cent.loc[expr_sub.columns]

    # Per patient: Spearman correlation of patient's PAM50-gene vector with each centroid
    n_patients = expr_sub.shape[0]
    results = pd.DataFrame(index=expr_sub.index, columns=["Basal", "Her2", "LumA", "LumB", "Normal"], dtype=float)
    for pid in expr_sub.index:
        pvec = expr_sub.loc[pid].values
        for s in results.columns:
            cvec = cent_sub[s].values
            rho, _ = spearmanr(pvec, cvec)
            results.loc[pid, s] = rho

    # Assign subtype = argmax (tie-break: alphabetical class name, per pre-reg §3.1)
    def assign(row):
        m = row.max()
        cands = sorted([s for s in row.index if row[s] == m])
        return cands[0] if cands else "Unknown"
    subtype = results.apply(assign, axis=1)

    out = pd.DataFrame({"patient_id": expr.index})
    out["pam50_subtype"] = out["patient_id"].map(subtype.to_dict())
    out["basal_score"] = out["patient_id"].map(results["Basal"].to_dict())
    out["her2_score"] = out["patient_id"].map(results["Her2"].to_dict())
    out["luma_score"] = out["patient_id"].map(results["LumA"].to_dict())
    out["lumb_score"] = out["patient_id"].map(results["LumB"].to_dict())
    out["normal_score"] = out["patient_id"].map(results["Normal"].to_dict())

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False, float_format="%.6f")

    log["subtype_counts"] = out["pam50_subtype"].value_counts().to_dict()
    log["n_patients"] = int(n_patients)
    log["basal_score_mean"] = float(out["basal_score"].mean())
    log["basal_score_sd"] = float(out["basal_score"].std())
    with open(args.out_log, "w") as f:
        json.dump(log, f, indent=2)
    print(f"[{args.cohort}] PAM50: {n_mapped}/50 genes; subtypes={log['subtype_counts']}")


if __name__ == "__main__":
    main()
