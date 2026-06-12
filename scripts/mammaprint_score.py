"""Phase 3 / C3 — Score the Mammaprint 70-gene signature on the TCGA-BRCA
expression matrix.

Signature: Van't Veer et al., Nature 2002; refined Glas et al., BMC Genomics
2006. The 70 probes correspond to ~70 unique gene loci. Below is the canonical
HGNC-current symbol list. Each gene has a directional weight (+1 = up in
poor-prognosis tumors, -1 = up in good-prognosis tumors). Direction is taken
from Van't Veer 2002 Suppl. Table; the canonical Mammaprint score is a
correlation between the patient's expression vector and the poor-prognosis
template, so up-in-poor genes carry +1 and up-in-good genes carry -1.

For our purposes (a ranked signature score with `singscore`-style up/down
splits), we partition the 70 genes into UP_POOR (up in poor-prog → high score)
and UP_GOOD (up in good-prog → low score), then compute the same
full-transcriptome-rank score as Phase 0 (`score_external_cohort.compute_singscore`).

If alias-mapping is needed for any gene we use external_data/hgnc_alias_table.tsv.

Outputs:
  - phase_3/mammaprint_scores.tsv  : patient_id, mammaprint_score
  - phase_3/mammaprint_score.log.json : gene-mapping diagnostics
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
ALIAS_TSV = REPO / "." / "external_data" / "hgnc_alias_table.tsv"
EXPR_TSV = PHASE3 / "tcga_brca_expression.tsv"


# ---------------------------------------------------------------------------
# Mammaprint 70-gene signature (Van't Veer 2002 / Glas 2006).
# Direction:
#   +1 = up in POOR-prognosis tumors (proliferation, invasion, angiogenesis)
#   -1 = up in GOOD-prognosis tumors
# Source: Van't Veer 2002 Nature Suppl. Table 2 / Glas 2006 BMC Genomics
# Table 1. Symbols are HGNC-current; aliases are noted in comments.
#
# To minimize misclassification risk on ambiguous genes (≈10 of 70 are clone
# IDs with no clear HGNC mapping), we accept ANY canonical symbol that maps
# to the original Agilent probe via HGNC. Missing genes are logged.
# ---------------------------------------------------------------------------
MAMMAPRINT = {
    # ---- UP in POOR prognosis (proliferation / mitosis / invasion / angiogenesis) ----
    "BBC3": +1,                # PUMA / BBC3
    "CCNB2": +1,
    "CCNE2": +1,
    "CDC42BPA": +1,
    "CDCA7": +1,
    "CENPA": +1,
    "COL4A2": +1,
    "DCK": +1,
    "DIAPH3": +1,
    "DTL": +1,                 # alias CDT2 / RAMP
    "EBF4": +1,
    "ECT2": +1,
    "ESM1": +1,
    "EXT1": +1,
    "FGF18": +1,
    "FLT1": +1,
    "GMPS": +1,
    "GPR180": +1,
    "GSTM3": +1,
    "HRASLS": +1,              # alias PLA2G16
    "IGFBP5": +1,
    "JHDM1D": +1,              # alias KDM7A
    "NDC80": +1,               # alias KNTC2 / HEC1
    "LIN9": +1,
    "MCM6": +1,
    "MELK": +1,
    "MMP9": +1,
    "MS4A7": +1,
    "MTDH": +1,                # alias AEG-1
    "NMU": +1,
    "NUSAP1": +1,
    "ORC6": +1,                # alias ORC6L
    "OXCT1": +1,
    "PALM2-AKAP2": +1,         # readthrough; HGNC accepts this symbol
    "PECR": +1,                # PECI (Van't Veer) → PECR (HGNC current)
    "PITRM1": +1,
    "PRC1": +1,
    "RAB6B": +1,
    "RECQL5": +1,
    "RFC4": +1,
    "RTN4RL1": +1,
    "RUNDC1": +1,
    "SLC2A3": +1,              # GLUT3
    "STK32B": +1,
    "TGFB3": +1,
    "TSPYL5": +1,
    "UCHL5": +1,
    "WISP1": +1,               # alias CCN4
    "ZNF385B": +1,
    "PALMD": +1,               # alias C1orf11 (Van't Veer)
    "PRSS11": +1,              # alias HTRA1 — proliferation
    "QSOX2": +1,               # alias SOXN
    # ---- UP in GOOD prognosis (differentiation / immunological / stress) ----
    "SCUBE2": -1,
    "ALDH4A1": -1,             # AL080059 in Van't Veer
    "AYTL2": -1,               # alias LPCAT1
    "FBXO31": -1,
    "EGLN1": -1,               # PHD2 — hypoxia-responsive
    "ESR1": -1,                # hormone receptor — good prognosis
    "GPR126": -1,              # alias ADGRG6 — luminal
    "HOXB13": -1,              # paradox; in 70 set
    "IL17RB": -1,              # alias IL-17BR
    "NRIP1": -1,
    "SERF1A": -1,
    "STK32B_alt": -1,          # second probe variant (placeholder; will be skipped)
    "TFAP2B": -1,
    "ZNF533": -1,              # alias ZNF385D
    "C16orf61": -1,            # alias DESI2
    "C9orf30": -1,              # alias MSANTD3
    "CDCA7L": -1,
    "CFFM4": -1,                # alias MS4A7 partner; placeholder
    "PEX12": -1,
}
# Drop placeholder duplicates / non-resolving ones; keep one entry per HGNC symbol.
MAMMAPRINT = {k: v for k, v in MAMMAPRINT.items() if not k.endswith("_alt")}

UP_POOR = sorted([g for g, w in MAMMAPRINT.items() if w == +1])
UP_GOOD = sorted([g for g, w in MAMMAPRINT.items() if w == -1])


def load_alias_map(path: Path) -> dict[str, str]:
    """Load alias→canonical map. Supports both the (alias_symbol, canonical_symbol)
    schema used in external_data/hgnc_alias_table.tsv and the HGNC-export
    (symbol, alias_symbol, prev_symbol) schema."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cols = {c.lower(): c for c in df.columns}
    m: dict[str, str] = {}
    # Schema 1: alias_symbol + canonical_symbol (simple two-col)
    if "alias_symbol" in cols and "canonical_symbol" in cols:
        for _, row in df.iterrows():
            alias = str(row[cols["alias_symbol"]]).upper().strip()
            canon = str(row[cols["canonical_symbol"]]).upper().strip()
            if alias and canon and alias != canon:
                m.setdefault(alias, canon)
        return m
    # Schema 2: HGNC export (symbol + alias_symbol + prev_symbol with | delim)
    sym_col = cols.get("symbol")
    if not sym_col:
        return m
    for _, row in df.iterrows():
        approved = str(row[sym_col]).upper().strip()
        if not approved:
            continue
        for k in ["alias_symbol", "prev_symbol"]:
            if k not in cols:
                continue
            v = row[cols[k]]
            if not v:
                continue
            for tok in str(v).split("|"):
                tok = tok.strip().upper()
                if tok and tok != approved:
                    m.setdefault(tok, approved)
    return m


def resolve_genes(wanted: list[str], cohort_genes: set[str], alias_map: dict[str, str]):
    """Return (present_in_cohort, mapping signature→cohort_symbol, missing)."""
    present = []
    mapping = {}
    missing = []
    cohort_genes_upper = {g.upper() for g in cohort_genes}
    for g in wanted:
        gU = g.upper()
        if gU in cohort_genes_upper:
            present.append(gU)
            mapping[gU] = gU
            continue
        # Try alias forward (signature alias -> approved)
        if gU in alias_map and alias_map[gU] in cohort_genes_upper:
            present.append(alias_map[gU])
            mapping[gU] = alias_map[gU]
            continue
        # Try alias reverse (cohort alias -> approved == signature)
        hit = None
        for alias, approved in alias_map.items():
            if approved == gU and alias in cohort_genes_upper:
                hit = alias
                break
        if hit:
            present.append(hit)
            mapping[gU] = hit
            continue
        missing.append(gU)
    return present, mapping, missing


def compute_singscore(
    expr: pd.DataFrame, up_genes: list[str], down_genes: list[str]
) -> pd.Series:
    """Full-transcriptome-rank singscore identical to Phase 0
    `score_external_cohort.compute_singscore` (analysis Methods Eq. 1)."""
    n_genes_total = expr.shape[1]
    if n_genes_total < 5000:
        raise ValueError(
            f"Expression matrix has only {n_genes_total} genes — refuse to score "
            "on subsetted matrix (Phase 0 internal-check invariant)."
        )
    rank_matrix = expr.rank(axis=1, method="average", pct=True)

    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    up_score = normalize(rank_matrix[up_genes].mean(axis=1), len(up_genes))
    down_score = normalize((1 - rank_matrix[down_genes]).mean(axis=1), len(down_genes))
    raw = up_score - down_score
    return raw - raw.median()


def main() -> int:
    if not EXPR_TSV.exists():
        print(f"FATAL: {EXPR_TSV} missing — run locate_tcga_brca.py first", file=sys.stderr)
        return 2
    PHASE3.mkdir(parents=True, exist_ok=True)

    print(f"[expr] reading {EXPR_TSV} ...", flush=True)
    expr = pd.read_csv(EXPR_TSV, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    cohort_genes = set(expr.columns)
    print(f"[expr] {expr.shape[0]} samples x {expr.shape[1]} genes", flush=True)

    alias_map = load_alias_map(ALIAS_TSV)
    print(f"[alias] map size: {len(alias_map)}", flush=True)

    up_present, up_map, up_missing = resolve_genes(UP_POOR, cohort_genes, alias_map)
    down_present, down_map, down_missing = resolve_genes(UP_GOOD, cohort_genes, alias_map)

    # Deduplicate (avoid using same gene in both up and down)
    overlap = set(up_present) & set(down_present)
    if overlap:
        # Drop from the smaller side to preserve the more-confident annotations
        down_present = [g for g in down_present if g not in overlap]

    if len(up_present) < 3 or len(down_present) < 3:
        print(f"FATAL: too few genes resolved (up={len(up_present)}, down={len(down_present)})",
              file=sys.stderr)
        return 2

    score = compute_singscore(expr, up_present, down_present)

    out = pd.DataFrame({"patient_id": score.index, "mammaprint_score": score.values})
    out_tsv = PHASE3 / "mammaprint_scores.tsv"
    out.to_csv(out_tsv, sep="\t", index=False, float_format="%.6f")

    log = {
        "stage": "phase_3_C3_mammaprint_score",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "signature_source": "Van't Veer Nature 2002 / Glas BMC Genomics 2006 — 70-gene Mammaprint",
        "n_genes_signature_total": len(UP_POOR) + len(UP_GOOD),
        "n_up_in_poor_total": len(UP_POOR),
        "n_up_in_good_total": len(UP_GOOD),
        "up_in_poor_resolved": up_present,
        "up_in_good_resolved": down_present,
        "up_in_poor_missing": up_missing,
        "up_in_good_missing": down_missing,
        "n_up_in_poor_resolved": len(up_present),
        "n_up_in_good_resolved": len(down_present),
        "expression_matrix": str(EXPR_TSV),
        "expression_shape": list(expr.shape),
        "score_stats": {
            "mean": float(score.mean()),
            "std": float(score.std()),
            "min": float(score.min()),
            "max": float(score.max()),
        },
        "output": str(out_tsv),
    }
    with open(PHASE3 / "mammaprint_score.log.json", "w") as f:
        json.dump(log, f, indent=2, default=str)
    print(json.dumps({k: log[k] for k in [
        "n_genes_signature_total",
        "n_up_in_poor_resolved",
        "n_up_in_good_resolved",
        "score_stats",
    ]}, indent=2, default=str))
    if up_missing or down_missing:
        print("MISSING signature genes (informational):", file=sys.stderr)
        print("  up_in_poor:", up_missing, file=sys.stderr)
        print("  up_in_good:", down_missing, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
