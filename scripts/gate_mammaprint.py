"""Phase 3 / C4 — Run the GATE diagnostic on the Mammaprint × TCGA-BRCA
worked example (§4.9 second worked example).

Primary analyses (PRIMARY worked example, full TCGA-BRCA cohort):
  - Mammaprint score vs 3-yr OS-or-recurrence outcome, stratified by PAM50.
  - Mammaprint score vs 5-yr OS-or-recurrence outcome, stratified by PAM50.

Sensitivity (chemo subset of TCGA-BRCA):
  - Same two outcome windows on the chemo-treated subset. Reported but flagged
    as under-powered (n_labeled<50).

Additional sensitivity (framework's Doxorubicin signature × TCGA-BRCA chemo cohort):
  - Score the frozen DOXORUBICIN.tsv signature on the chemo subset expression
    matrix; run GATE on the result. This is the §4.9 cross-check requested in
    the study spec.

Decision rule (locked):
  Layer-1 portability CONFIRMED if:
      GAP >= 0.10 AND within-stratum weighted mean AUROC ≤ 0.55.
  Otherwise CAVEATED.

Outputs (phase_3/):
  - gate_mammaprint_full_3yr.json
  - gate_mammaprint_full_5yr.json
  - gate_mammaprint_chemo_3yr.json
  - gate_mammaprint_chemo_5yr.json
  - gate_doxorubicin_chemo_3yr.json
  - gate_doxorubicin_chemo_5yr.json
  - gate_mammaprint.log.json   : aggregated headline summary
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(".")
PHASE3 = REPO / "." / "." / "results" / "phase_3"
SCRIPTS = REPO / "." / "." / "scripts"
sys.path.insert(0, str(SCRIPTS))
from gate_diagnostic import gate  # noqa: E402

DOX_SIGNATURE = REPO / "." / "frozen_signatures" / "quantile_30_oncokb_all" / "DOXORUBICIN.tsv"


def load_clinical_labels(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["patient_id"] = df["patient_id"].astype(str)
    return df


def merge_score_label(score_df: pd.DataFrame, clin_df: pd.DataFrame) -> pd.DataFrame:
    s = score_df.rename(columns={score_df.columns[1]: "score"})[["patient_id", "score"]]
    return s.merge(clin_df, on="patient_id", how="inner")


def run_gate(merged: pd.DataFrame, label_col: str) -> dict:
    sub = merged[merged[label_col] != -1].copy()
    if len(sub) < 10:
        return {"n": int(len(sub)), "verdict": "INSUFFICIENT_n"}
    y = sub[label_col].astype(int).to_numpy()
    s = sub["score"].astype(float).to_numpy()
    strat = sub["PAM50"].where(sub["PAM50"].notna()).astype(object).to_numpy()
    res = gate(s, y, strat, n_boot=1000, seed=42)
    # JSON-safe
    res["stratum_sizes"] = {str(k): int(v) for k, v in res["stratum_sizes"].items()}
    res["within_stratum_AUROC"] = {str(k): (None if (v is None or (isinstance(v, float) and not np.isfinite(v))) else float(v))
                                    for k, v in res["within_stratum_AUROC"].items()}
    res["marginal_CI"] = [float(v) for v in res["marginal_CI"]]
    return res


def decision_rule(gate_res: dict) -> str:
    """Layer 1 portability rule (locked, pre-reg §7.1 + §4.9 plan):
       CONFIRMED if GAP >= 0.10 AND weighted mean within-stratum AUROC ≤ 0.55.
       CAVEATED otherwise."""
    if "GAP" not in gate_res or gate_res["GAP"] is None:
        return "INSUFFICIENT"
    gap = float(gate_res["GAP"])
    wmean = gate_res.get("within_stratum_AUROC", {}).get("_weighted_mean")
    if wmean is None or not np.isfinite(wmean):
        return "INSUFFICIENT"
    if gap >= 0.10 and wmean <= 0.55:
        return "CONFIRMED_PORTABILITY"
    return "CAVEATED"


def score_doxorubicin(expr_path: Path, sig_path: Path) -> pd.DataFrame:
    """Score the frozen DOXORUBICIN signature on the expression matrix.
    Mirrors `score_external_cohort.compute_singscore` semantics."""
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    sig = pd.read_csv(sig_path, sep="\t")
    sens = set(sig[sig["direction"] == "sensitivity"]["gene_symbol"].str.upper())
    res = set(sig[sig["direction"] == "resistance"]["gene_symbol"].str.upper())

    n_genes_total = expr.shape[1]
    rank_matrix = expr.rank(axis=1, method="average", pct=True)
    up_present = sorted(sens & set(expr.columns))
    down_present = sorted(res & set(expr.columns))
    if len(up_present) < 3 or len(down_present) < 3:
        raise ValueError("Doxorubicin signature did not resolve enough genes")

    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    up_score = normalize(rank_matrix[up_present].mean(axis=1), len(up_present))
    down_score = normalize((1 - rank_matrix[down_present]).mean(axis=1), len(down_present))
    raw = up_score - down_score
    centered = raw - raw.median()
    return pd.DataFrame({"patient_id": centered.index, "dox_score": centered.values})


def main() -> int:
    PHASE3.mkdir(parents=True, exist_ok=True)

    mammaprint_path = PHASE3 / "mammaprint_scores.tsv"
    full_clin_path = PHASE3 / "tcga_brca_full_cohort.tsv"
    chemo_clin_path = PHASE3 / "tcga_brca_chemo_subset.tsv"
    expr_path = PHASE3 / "tcga_brca_expression.tsv"

    for p in (mammaprint_path, full_clin_path, chemo_clin_path, expr_path):
        if not p.exists():
            print(f"FATAL: required input {p} missing", file=sys.stderr)
            return 2

    score = pd.read_csv(mammaprint_path, sep="\t")
    full = load_clinical_labels(full_clin_path)
    chemo = load_clinical_labels(chemo_clin_path)

    merged_full = merge_score_label(score, full)
    merged_chemo = merge_score_label(score, chemo)

    headline: dict[str, dict] = {}

    # --- PRIMARY: Mammaprint on full TCGA-BRCA cohort ---
    for window in ["label_3yr", "label_5yr"]:
        key = f"mammaprint_full_{window}"
        res = run_gate(merged_full, window)
        res["decision"] = decision_rule(res)
        with open(PHASE3 / f"gate_{key}.json", "w") as f:
            json.dump(res, f, indent=2, default=str)
        headline[key] = res

    # --- SENSITIVITY: Mammaprint on chemo subset ---
    for window in ["label_3yr", "label_5yr"]:
        key = f"mammaprint_chemo_{window}"
        res = run_gate(merged_chemo, window)
        res["decision"] = decision_rule(res)
        with open(PHASE3 / f"gate_{key}.json", "w") as f:
            json.dump(res, f, indent=2, default=str)
        headline[key] = res

    # --- SENSITIVITY: Doxorubicin frozen signature on chemo subset ---
    if DOX_SIGNATURE.exists():
        dox_scores = score_doxorubicin(expr_path, DOX_SIGNATURE)
        for window in ["label_3yr", "label_5yr"]:
            merged_dox = merge_score_label(dox_scores, chemo)
            key = f"doxorubicin_chemo_{window}"
            res = run_gate(merged_dox, window)
            res["decision"] = decision_rule(res)
            with open(PHASE3 / f"gate_{key}.json", "w") as f:
                json.dump(res, f, indent=2, default=str)
            headline[key] = res
        # Also score Dox on the FULL cohort for completeness (cross-check that
        # marginal AUROC for a CCLE-derived chemo signature is near 0.5 on a
        # broad prognosis endpoint, which would corroborate the GATE story).
        for window in ["label_3yr", "label_5yr"]:
            merged_dox_full = merge_score_label(dox_scores, full)
            key = f"doxorubicin_full_{window}"
            res = run_gate(merged_dox_full, window)
            res["decision"] = decision_rule(res)
            with open(PHASE3 / f"gate_{key}.json", "w") as f:
                json.dump(res, f, indent=2, default=str)
            headline[key] = res
    else:
        print(f"WARN: {DOX_SIGNATURE} missing; skipping Dox sensitivity", file=sys.stderr)

    log = {
        "stage": "phase_3_C4_gate_mammaprint",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "mammaprint_scores": str(mammaprint_path),
            "full_cohort_labels": str(full_clin_path),
            "chemo_subset_labels": str(chemo_clin_path),
            "expression": str(expr_path),
            "doxorubicin_signature": str(DOX_SIGNATURE),
        },
        "decision_rule": "Portability CONFIRMED iff GAP >= 0.10 AND weighted-within-stratum AUROC <= 0.55",
        "results": headline,
    }
    with open(PHASE3 / "gate_mammaprint.log.json", "w") as f:
        json.dump(log, f, indent=2, default=str)

    # Short stdout summary.
    print("\n=== Phase 3 GATE summary ===")
    for key, res in headline.items():
        marg = res.get("marginal_AUROC", float("nan"))
        ci = res.get("marginal_CI", [None, None])
        wmean = res.get("within_stratum_AUROC", {}).get("_weighted_mean")
        gap = res.get("GAP")
        n = res.get("n_total")
        p = res.get("p_value")
        decision = res.get("decision")
        print(f"{key:36s}  n={n}  marg={marg:.3f} CI[{ci[0]:.3f},{ci[1]:.3f}]  "
              f"within_w_mean={('NA' if wmean is None else f'{wmean:.3f}')}  "
              f"GAP={('NA' if gap is None else f'{gap:.3f}')}  "
              f"p={('NA' if p is None else f'{p:.3f}')}  -> {decision}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
