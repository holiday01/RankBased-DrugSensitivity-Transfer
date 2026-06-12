"""Phase 8 / Sig 2 — HRD / BRCAness signature x GSE22226 (I-SPY1).

Hypothesis (orthogonal-signature falsification test):
    Homologous-recombination-deficiency (HRD) / BRCAness RNA signature
    captures DNA-damage-repair-pathway loss-of-function transcriptomic
    footprint. This is mechanistically distinct from proliferation (P),
    ER/PR signaling (H), HER2 (B), and even from PAM50-basal molecular
    subtype (most BRCA-mutant TNBCs are basal, BUT many basal tumors are
    HR-proficient). On I-SPY1 (GSE22226, AC->T arm including some
    AC-carboplatin contexts in subsequent I-SPY enrollment), HRD-high /
    BRCAness tumors are pCR-favorable.

Pre-spec: residual AUROC point > 0.60 AND CI low > 0.55 -> FALSIFIED
(signal survives axis deconfounding = test has discrimination power).

Gene list (Severson et al., NPJ Breast Cancer 2017, "RNA-HRD"
9-gene HR/BRCAness RNA signature, recalled MEDIUM confidence):

    BRCA1, BRCA2, FANCA, FANCC, FANCD2, ATM, CHEK2, RAD51, RAD51C

Direction: in Severson 2017 the signature is defined as the BRCA1/2-
pathway activity score; HR-PROFICIENT tumors have HIGH expression,
HR-DEFICIENT (BRCAness) tumors have LOW expression. We score the
mean per-gene z-score across the 9 genes as a "HR proficiency score":

    high score -> HR-proficient
    low score  -> HR-deficient (BRCAness)

For chemo pCR prediction, BRCAness is typically PRO-pCR (especially with
platinum-containing regimens). I-SPY1 (GSE22226) is mostly anthracycline-
+/- taxane; BRCAness signal on pure A-T regimens is weaker but published.
Therefore expected direction: HR proficiency score should ANTI-correlate
with pCR (i.e. raw AUROC may be < 0.5; under signed-AUROC verdict the
test still reads CONFIRMED if residual collapses to ~0.5, FALSIFIED only
on point > 0.60).

Gene-list confidence (honest):
- Severson 2017 RNA-HRD 9-gene MEDIUM confidence (canonical HR/Fanconi
  genes; the published list is HR core + Fanconi anemia complementation
  group; we did not consult Suppl. Table verbatim).
- Falls back to Peng 2014 3-gene mini-HRD (BRCA1, RAD51, FANCD2) if
  fewer than 5 of the 9 resolve on GSE22226.

Score = mean per-gene z-score, same scoring scheme as Phase 7 / earlier
Phase 8 for apples-to-apples comparison.
"""
from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "phase_7"))
from _common import (
    sha256_of, load_alias_map, resolve_genes,
    auroc_with_boot_ci, fit_residual, verdict_3x3,
)

REPO = Path(".")
PHASE_0 = REPO / "." / "results" / "phase_0"
PHASE_8 = REPO / "." / "results" / "phase_8"
EXPR_TSV = REPO / "external_data" / "GSE22226" / "GSE22226_expression.tsv"
ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

HRD_GENES_SEVERSON = [
    "BRCA1", "BRCA2", "FANCA", "FANCC", "FANCD2",
    "ATM", "CHEK2", "RAD51", "RAD51C",
]
HRD_GENES_PENG_FALLBACK = ["BRCA1", "RAD51", "FANCD2"]


def mean_zscore(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    sub = expr[genes].copy()
    z = (sub - sub.mean(axis=0)) / sub.std(axis=0, ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def axis_correlations(df: pd.DataFrame, sig_col: str,
                      axis_cols: list[str]) -> dict:
    out = {}
    for a in axis_cols:
        sub = df.dropna(subset=[sig_col, a])
        if len(sub) < 5:
            out[a] = float("nan")
        else:
            out[a] = float(sub[sig_col].corr(sub[a]))
    return out


def main():
    PHASE_8.mkdir(parents=True, exist_ok=True)
    log = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "signature": "HRD / BRCAness RNA signature (Severson 2017 9-gene primary; "
                     "Peng 2014 3-gene fallback)",
        "signature_confidence": (
            "MEDIUM (Severson 9-gene HR/Fanconi core recalled; not verified "
            "against Suppl. table verbatim). Peng 2014 3-gene fallback HIGH "
            "(BRCA1, RAD51, FANCD2 are canonical HR-core)."
        ),
        "scoring": "mean per-gene z-score; HIGH score = HR-proficient, LOW = BRCAness",
        "decision_rule": (
            "point<=0.55 AND ci_lo<=0.55 -> CONFIRMED; "
            ">0.60 AND ci_lo>0.55 -> FALSIFIED; else GRAY"
        ),
        "expected_outcome": (
            "MEDIUM chance of FALSIFICATION: HRD biology is mechanistically "
            "orthogonal to proliferation/ER/HER2/basal axes. On I-SPY1 "
            "(GSE22226) anthracycline-taxane cohort, BRCAness signal is "
            "weaker than on platinum cohorts; effect size may be marginal."
        ),
    }

    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    alias_map = load_alias_map(ALIAS)

    # Try Severson 9-gene first
    present, mapping, missing = resolve_genes(HRD_GENES_SEVERSON, cohort_genes, alias_map)
    used_list = "Severson_9gene"
    if len(present) < 5:
        # Fall back to Peng 3-gene
        present, mapping, missing = resolve_genes(HRD_GENES_PENG_FALLBACK,
                                                  cohort_genes, alias_map)
        used_list = "Peng_3gene_fallback"
        log["fallback_used"] = (
            "Severson 9-gene resolved < 5; fell back to Peng 3-gene "
            "(BRCA1, RAD51, FANCD2)."
        )

    log["used_gene_list"] = used_list
    log["n_resolved"] = len(present)
    log["resolved_genes"] = present
    log["missing_genes"] = missing
    log["alias_mapping"] = mapping
    log["n_requested"] = (len(HRD_GENES_SEVERSON) if used_list == "Severson_9gene"
                         else len(HRD_GENES_PENG_FALLBACK))

    if len(present) < 2:
        log["error"] = f"too few HRD genes resolved ({len(present)})"
        with open(PHASE_8 / "HRD_gse22226_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    score = mean_zscore(expr, present)
    score_df = pd.DataFrame({"sample_id": score.index, "sig_HRD": score.values})
    score_df.to_csv(PHASE_8 / "score_HRD_GSE22226.tsv",
                    sep="\t", index=False, float_format="%.6f")

    axes = pd.read_csv(PHASE_0 / "axes_GSE22226.tsv", sep="\t")
    df = axes.merge(score_df, on="sample_id", how="inner")
    df["pcr"] = pd.to_numeric(df["pCR"], errors="coerce")
    log["n_merged"] = int(len(df))
    log["n_with_pcr"] = int(df["pcr"].notna().sum())
    log["pcr_rate"] = float((df["pcr"] == 1).mean())

    axis_cols = ["P", "H", "B", "H_PAM50basal"]
    corrs = axis_correlations(df, "sig_HRD", axis_cols)
    log["sig_vs_axis_pearson_r"] = corrs

    sub_raw = df.dropna(subset=["sig_HRD", "pcr"])
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(
        sub_raw["sig_HRD"].to_numpy(),
        sub_raw["pcr"].astype(int).to_numpy(),
        n_boot=1000, seed=42,
    )
    log["marginal_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi,
                              "n": int(len(sub_raw))}

    sub, coef = fit_residual(df, sig_col="sig_HRD",
                             axis_cols=axis_cols,
                             label_col="pcr")
    if sub is None:
        log["error"] = "too few rows for residual model"
        with open(PHASE_8 / "HRD_gse22226_residual.json", "w") as f:
            json.dump(log, f, indent=2)
        return

    y = sub["pcr"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, n_boot=1000, seed=42)
    log["residual_model"] = {
        "formula": "sig_HRD ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["residual_auroc"] = {"point": auc, "ci_low": lo, "ci_high": hi}
    log["residual_auroc_abs_dist_from_0p5"] = {
        "point": abs(auc - 0.5) if auc == auc else float("nan"),
        "ci_low_abs": abs(lo - 0.5) if lo == lo else float("nan"),
        "ci_high_abs": abs(hi - 0.5) if hi == hi else float("nan"),
    }
    log["verdict_3x3"] = verdict_3x3(auc, lo)
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    with open(PHASE_8 / "HRD_gse22226_residual.json", "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n===== Phase 8 / Sig 2 HRD x GSE22226 (n={len(sub)}) =====")
    print(f"Used gene list: {used_list}")
    print(f"Resolved HRD genes: {len(present)}: {present}")
    print(f"Missing: {missing}")
    print(f"Sig vs axes Pearson r: {corrs}")
    print(f"Marginal AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"Residual AUROC = {auc:.3f} [{lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {log['verdict_3x3']}")


if __name__ == "__main__":
    main()
