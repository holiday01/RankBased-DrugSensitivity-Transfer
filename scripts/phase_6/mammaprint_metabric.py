"""Phase 6 / Sig 4 — Mammaprint 70-gene × METABRIC chemo subset.

Reuses the EXACT Mammaprint gene list from Phase 3 (mammaprint_score.py)
applied to the METABRIC expression matrix restricted to the chemo subset
(CHEMOTHERAPY=YES, with 5-y OS endpoint, n=399 — defined in Phase 4).

Pre-specified pipeline:
    sig_Mammaprint ~ P + H + B + H_PAM50basal  →  residual vs event5y (5-y OS).
3×3 verdict (same rule as Phase 4 metabric).
"""
from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _common import (
    sha256_of, load_alias_map, resolve_genes, compute_singscore,
    auroc_with_boot_ci, fit_residual, verdict_3x3,
)

REPO = Path(".")
PHASE_4M = REPO / "." / "results" / "phase_4" / "metabric"
PHASE_6 = REPO / "." / "results" / "phase_6"
EXPR_TSV = PHASE_4M / "METABRIC_expression.tsv"
AXES_TSV = PHASE_4M / "axes_METABRIC.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

# Copied verbatim from Phase 3 mammaprint_score.py to guarantee parity.
MAMMAPRINT = {
    "BBC3": +1, "CCNB2": +1, "CCNE2": +1, "CDC42BPA": +1, "CDCA7": +1,
    "CENPA": +1, "COL4A2": +1, "DCK": +1, "DIAPH3": +1, "DTL": +1,
    "EBF4": +1, "ECT2": +1, "ESM1": +1, "EXT1": +1, "FGF18": +1,
    "FLT1": +1, "GMPS": +1, "GPR180": +1, "GSTM3": +1, "HRASLS": +1,
    "IGFBP5": +1, "JHDM1D": +1, "NDC80": +1, "LIN9": +1, "MCM6": +1,
    "MELK": +1, "MMP9": +1, "MS4A7": +1, "MTDH": +1, "NMU": +1,
    "NUSAP1": +1, "ORC6": +1, "OXCT1": +1, "PALM2-AKAP2": +1, "PECR": +1,
    "PITRM1": +1, "PRC1": +1, "RAB6B": +1, "RECQL5": +1, "RFC4": +1,
    "RTN4RL1": +1, "RUNDC1": +1, "SLC2A3": +1, "STK32B": +1, "TGFB3": +1,
    "TSPYL5": +1, "UCHL5": +1, "WISP1": +1, "ZNF385B": +1, "PALMD": +1,
    "PRSS11": +1, "QSOX2": +1,
    "SCUBE2": -1, "ALDH4A1": -1, "AYTL2": -1, "FBXO31": -1,
    "EGLN1": -1, "ESR1": -1, "GPR126": -1, "HOXB13": -1,
    "IL17RB": -1, "NRIP1": -1, "SERF1A": -1, "TFAP2B": -1,
    "ZNF533": -1, "C16orf61": -1, "C9orf30": -1, "CDCA7L": -1,
    "CFFM4": -1, "PEX12": -1,
}
UP_POOR = sorted([g for g, w in MAMMAPRINT.items() if w == +1])
UP_GOOD = sorted([g for g, w in MAMMAPRINT.items() if w == -1])


def main():
    PHASE_6.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "Mammaprint 70-gene (Van't Veer 2002 / Glas 2006) — same gene list as Phase 3",
        "expression_path": str(EXPR_TSV),
        "cohort_filter": "METABRIC chemo subset (CHEMOTHERAPY=YES, with 5-y OS); defined in Phase 4",
        "decision_rule": "point<=0.55 AND ci_lo<=0.55 → CONFIRMED; >0.60 AND ci_lo>0.55 → FALSIFIED; else GRAY",
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    expr.index = expr.index.astype(str)
    log["expression_shape"] = list(expr.shape)
    alias_map = load_alias_map(ALIAS)

    up_present, _, up_missing = resolve_genes(UP_POOR, set(expr.columns), alias_map)
    down_present, _, down_missing = resolve_genes(UP_GOOD, set(expr.columns), alias_map)
    overlap = set(up_present) & set(down_present)
    if overlap:
        down_present = [g for g in down_present if g not in overlap]
    log["n_up_in_poor_resolved"] = len(up_present)
    log["n_up_in_good_resolved"] = len(down_present)
    log["up_missing"] = up_missing
    log["down_missing"] = down_missing

    score = compute_singscore(expr, up_present, down_present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_Mammaprint": score.values})
    score_df.to_csv(PHASE_6 / "score_Mammaprint_METABRIC.tsv", sep="\t", index=False,
                    float_format="%.6f")

    # Axes — restrict to chemo subset is already built into axes_METABRIC.tsv (412)
    axes = pd.read_csv(AXES_TSV, sep="\t")
    axes["sample_id"] = axes["sample_id"].astype(str)
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["event"] = pd.to_numeric(df["event5y"], errors="coerce")
    log["n_chemo_with_axes"] = int(len(axes))
    log["n_merged"] = int(len(df))
    log["n_with_event"] = int(df["event"].notna().sum())
    log["event5y_rate"] = float(df["event"].mean(skipna=True))

    sub_raw = df.dropna(subset=["sig_Mammaprint", "event"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_Mammaprint"].to_numpy(), sub_raw["event"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_Mammaprint",
                             axis_cols=["P", "H", "B", "H_PAM50basal"],
                             label_col="event")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_6 / "Mammaprint_metabric_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["event"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_Mammaprint ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_6 / "Mammaprint_metabric_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 6 / Sig 4 Mammaprint × METABRIC chemo (n={len(sub)}) =====")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
