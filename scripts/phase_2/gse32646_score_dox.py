"""Phase 2 / B3 — Apply frozen DOXORUBICIN signature to GSE32646 via singscore.

Quantile-normalize expression across samples first (per pre-reg deployment
recipe), then full-transcriptome rank-based singscore (matches Methods Eq. 1).

Output: phase_2/score_DOXORUBICIN_GSE32646.tsv (per sample)
        phase_2/score_DOXORUBICIN_GSE32646.json (diagnostics)
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
PHASE_2 = REPO / "." / "results" / "phase_2"
SIG_PATH = REPO / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all" / "DOXORUBICIN.tsv"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Column-wise quantile normalization (genes are rows in our internal view).

    Input here: rows = samples, cols = genes. We QN across samples (so each
    sample shares the same value distribution). Standard recipe: transpose,
    QN columns (= genes), transpose back.
    """
    # Work on samples-x-genes; QN across samples for each gene's rank vector
    # Standard approach: for matrix X (samples x genes), sort each sample's
    # gene-rank vector and replace with the row-wise mean of sorted rows.
    arr = df.values
    sorted_idx = np.argsort(arr, axis=1)
    sorted_vals = np.sort(arr, axis=1)
    mean_sorted = sorted_vals.mean(axis=0)
    qn = np.empty_like(arr, dtype=float)
    for i in range(arr.shape[0]):
        qn[i, sorted_idx[i]] = mean_sorted
    return pd.DataFrame(qn, index=df.index, columns=df.columns)


def load_signature(path: Path):
    df = pd.read_csv(path, sep="\t")
    sens = set(df[df["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(df[df["direction"] == "resistance"]["gene_symbol"].str.upper())
    info = {
        "signature_path": str(path),
        "signature_sha256": sha256_of(path),
        "n_sensitivity_genes": len(sens),
        "n_resistance_genes": len(res),
    }
    return sens, res, info


def singscore(expr: pd.DataFrame, up: set, down: set):
    expr.columns = [c.upper() for c in expr.columns]
    n_genes_total = expr.shape[1]
    if n_genes_total < 5000:
        raise ValueError(f"Only {n_genes_total} genes; refuse to score on subsetted matrix.")
    rank_matrix = expr.rank(axis=1, method="average", pct=True)
    up_present = sorted(up & set(expr.columns))
    down_present = sorted(down & set(expr.columns))
    if len(up_present) < 3 or len(down_present) < 3:
        raise ValueError(
            f"Too few sig genes present: up={len(up_present)}, down={len(down_present)}"
        )

    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    up_score = normalize(rank_matrix[up_present].mean(axis=1), len(up_present))
    down_score = normalize((1 - rank_matrix[down_present]).mean(axis=1), len(down_present))
    raw = up_score - down_score
    centered = raw - raw.median()
    diag = {
        "n_genes_in_expression": n_genes_total,
        "n_up_in_sig": len(up),
        "n_up_present": len(up_present),
        "n_down_in_sig": len(down),
        "n_down_present": len(down_present),
        "up_missing": sorted(up - set(expr.columns)),
        "down_missing": sorted(down - set(expr.columns)),
    }
    return centered, diag


def main():
    PHASE_2.mkdir(parents=True, exist_ok=True)
    expr_path = PHASE_2 / "GSE32646_expression.tsv"
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]

    log = {
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "expression_path": str(expr_path),
        "expression_sha256": sha256_of(expr_path),
    }

    print(f"[B3] expr shape: {expr.shape}; quantile-normalizing ...", flush=True)
    expr_qn = quantile_normalize(expr)
    log["quantile_normalized"] = True

    sens, res, sig_info = load_signature(SIG_PATH)
    log.update(sig_info)

    scores, diag = singscore(expr_qn, sens, res)
    log.update(diag)

    pcr_df = pd.read_csv(PHASE_2 / "GSE32646_pcr.tsv", sep="\t")
    pcr_df["patient_id"] = pcr_df["patient_id"].astype(str)
    common = sorted(set(scores.index.astype(str)) & set(pcr_df["patient_id"]))
    pcr_map = dict(zip(pcr_df["patient_id"], pcr_df["pcr"]))

    out = pd.DataFrame(
        {
            "patient_id": common,
            "sig_dox": scores.loc[common].values,
            "pcr": [pcr_map.get(p, np.nan) for p in common],
        }
    )
    out_path = PHASE_2 / "score_DOXORUBICIN_GSE32646.tsv"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.10g")
    log["score_output"] = str(out_path)
    log["n_patients_scored"] = int(len(common))
    log["sig_dox_mean"] = float(out["sig_dox"].mean())
    log["sig_dox_sd"] = float(out["sig_dox"].std())
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    log_path = PHASE_2 / "score_DOXORUBICIN_GSE32646.json"
    with open(log_path, "w") as f:
        json.dump(log, f, indent=2)

    print(f"[B3] Scored n={len(common)}; up_present={diag['n_up_present']}/{diag['n_up_in_sig']}, "
          f"down_present={diag['n_down_present']}/{diag['n_down_in_sig']}")
    print(f"    sig_dox mean={log['sig_dox_mean']:.4f} sd={log['sig_dox_sd']:.4f}")


if __name__ == "__main__":
    main()
