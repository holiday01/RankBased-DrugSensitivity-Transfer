"""Phase 0 — Compute axes P (proliferation), H (hormone-receptor), B (ERBB2)
on 3 external BC cohorts (GSE16446 / GSE25066 / GSE22226).

P axis = mean z-score across:
  HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT, HALLMARK_MITOTIC_SPINDLE,
  HALLMARK_MYC_TARGETS_V1, Whitfield_G2M (if available)
H axis        = mean z-score of (ESR1 + PGR)  with HGNC alias resolution
H_ESR1only    = ESR1 z-score
H_PAM50basal  = basal_score column from existing PAM50 output
B axis        = ERBB2 z-score

Singscore convention here is the same per-sample percentile-rank deployment
used in score_external_cohort.py (Methods Eq. 1) BUT for a gene set the
"down" set is empty, so the score reduces to:
   score = 2*(rank_pct.mean - mn)/(mx - mn) - 1
i.e., the unsigned singscore. We z-score across samples to make P comparable
across the 4–5 component sets before averaging.

Output: axes_<cohort>.tsv with columns
  sample_id, P, H, H_ESR1only, H_PAM50basal, B, pCR
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


REPO = Path(".")
COHORTS = ["GSE16446", "GSE25066", "GSE22226"]
HALLMARK_SETS = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_MYC_TARGETS_V1",
]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_gmt(path: Path) -> dict:
    out: dict[str, list[str]] = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g.strip().upper() for g in parts[2:] if g.strip()]
            out[name] = genes
    return out


def load_alias(path: Path) -> dict:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["alias_symbol"].astype(str).str.upper(),
                    df["canonical_symbol"].astype(str).str.upper()))


def resolve_gene(symbol: str, expr_cols: set, alias: dict) -> str | None:
    s = symbol.upper()
    if s in expr_cols:
        return s
    a = alias.get(s)
    if a and a in expr_cols:
        return a
    return None


def singscore_unsigned(rank_pct: pd.DataFrame, gene_set: list[str]) -> pd.Series:
    """Per-sample unsigned singscore for an "up only" gene set.
    rank_pct: DataFrame samples x genes (already per-sample percentile rank).
    """
    n_genes_total = rank_pct.shape[1]
    present = [g for g in gene_set if g in rank_pct.columns]
    if len(present) < 3:
        return pd.Series(np.nan, index=rank_pct.index)
    n_set = len(present)
    mn = (1 + n_set) / (2 * n_genes_total)
    mx = 1 - mn
    raw = rank_pct[present].mean(axis=1)
    norm = 2 * (raw - mn) / (mx - mn) - 1
    return norm


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (s - mu) / sd


def compute_axes_for_cohort(
    cohort: str,
    hallmark: dict,
    whitfield: list | None,
    alias: dict,
    log: dict,
) -> pd.DataFrame:
    cohort_dir = REPO / "external_data" / cohort
    expr_path = cohort_dir / f"{cohort}_expression.tsv"
    pam50_path = cohort_dir / f"{cohort}_pam50.tsv"
    pcr_path = cohort_dir / f"{cohort}_pcr.tsv"

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]

    n_samples_expr = expr.shape[0]
    n_genes = expr.shape[1]
    if n_genes < 5000:
        raise RuntimeError(f"{cohort}: only {n_genes} genes; refuse to score on subsetted matrix.")

    # Per-sample percentile ranks for unsigned singscore
    rank_pct = expr.rank(axis=1, method="average", pct=True)

    # --- P axis component sets ------------------------------------------------
    comp_z_cols: dict[str, pd.Series] = {}
    p_set_diag: dict[str, dict] = {}
    for name in HALLMARK_SETS:
        gs = hallmark.get(name, [])
        gs_resolved = []
        missing = []
        expr_cols_set = set(expr.columns)
        for g in gs:
            r = resolve_gene(g, expr_cols_set, alias)
            if r is None:
                missing.append(g)
            else:
                gs_resolved.append(r)
        gs_resolved = sorted(set(gs_resolved))
        score = singscore_unsigned(rank_pct, gs_resolved)
        comp_z_cols[name] = zscore(score)
        p_set_diag[name] = {
            "n_genes_in_set": len(gs),
            "n_resolved_present": len(gs_resolved),
            "n_missing": len(missing),
        }

    if whitfield is not None:
        gs_resolved = []
        missing = []
        expr_cols_set = set(expr.columns)
        for g in whitfield:
            r = resolve_gene(g, expr_cols_set, alias)
            if r is None:
                missing.append(g)
            else:
                gs_resolved.append(r)
        gs_resolved = sorted(set(gs_resolved))
        score = singscore_unsigned(rank_pct, gs_resolved)
        comp_z_cols["Whitfield_G2M"] = zscore(score)
        p_set_diag["Whitfield_G2M"] = {
            "n_genes_in_set": len(whitfield),
            "n_resolved_present": len(gs_resolved),
            "n_missing": len(missing),
        }

    p_df = pd.DataFrame(comp_z_cols)
    P = p_df.mean(axis=1)

    # --- H axis: z-score ESR1 + PGR using expression values directly ----------
    expr_cols_set = set(expr.columns)
    esr1 = resolve_gene("ESR1", expr_cols_set, alias)
    pgr = resolve_gene("PGR", expr_cols_set, alias)
    erbb2 = resolve_gene("ERBB2", expr_cols_set, alias)

    H_components = []
    h_diag = {"ESR1_symbol": esr1, "PGR_symbol": pgr, "ERBB2_symbol": erbb2}
    if esr1 is not None:
        H_ESR1_z = zscore(expr[esr1])
        H_components.append(H_ESR1_z)
    else:
        H_ESR1_z = pd.Series(np.nan, index=expr.index)
    if pgr is not None:
        H_PGR_z = zscore(expr[pgr])
        H_components.append(H_PGR_z)
    if H_components:
        H = pd.concat(H_components, axis=1).mean(axis=1)
    else:
        H = pd.Series(np.nan, index=expr.index)

    # --- H_PAM50basal: basal_score from PAM50 output --------------------------
    if pam50_path.exists():
        pam50 = pd.read_csv(pam50_path, sep="\t")
        pam50 = pam50.set_index("patient_id")
        H_PAM50basal = pam50["basal_score"].reindex(expr.index)
    else:
        H_PAM50basal = pd.Series(np.nan, index=expr.index)
        h_diag["pam50_missing"] = True

    # --- B axis ---------------------------------------------------------------
    if erbb2 is not None:
        B = zscore(expr[erbb2])
    else:
        B = pd.Series(np.nan, index=expr.index)

    # --- pCR -----------------------------------------------------------------
    pcr_df = pd.read_csv(pcr_path, sep="\t")
    pcr_df["patient_id"] = pcr_df["patient_id"].astype(str)
    pcr_df = pcr_df.set_index("patient_id")
    pCR = pd.to_numeric(pcr_df["pcr"], errors="coerce").reindex(expr.index)

    out = pd.DataFrame({
        "sample_id": expr.index,
        "P": P.values,
        "H": H.values,
        "H_ESR1only": H_ESR1_z.values,
        "H_PAM50basal": H_PAM50basal.values,
        "B": B.values,
        "pCR": pCR.values,
    })

    log[cohort] = {
        "n_samples_expression": int(n_samples_expr),
        "n_genes_expression": int(n_genes),
        "p_set_diagnostics": p_set_diag,
        "h_axis_diagnostics": h_diag,
        "n_samples_with_pCR": int(out["pCR"].notna().sum()),
        "n_samples_with_P": int(out["P"].notna().sum()),
        "n_samples_with_H": int(out["H"].notna().sum()),
        "n_samples_with_B": int(out["B"].notna().sum()),
        "n_samples_with_PAM50basal": int(out["H_PAM50basal"].notna().sum()),
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(REPO / "results/phase_0"))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    script_path = Path(__file__).resolve()
    script_sha = sha256_of(script_path)

    log_path = out_dir / "compute_axes.log"
    log_entries: dict = {
        "script": str(script_path),
        "script_sha256": script_sha,
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "cohorts": {},
    }

    hallmark_path = REPO / "h.all.v2025.1.Hs.symbols.gmt"
    hallmark = load_gmt(hallmark_path)
    log_entries["hallmark_sha256"] = sha256_of(hallmark_path)

    # Whitfield_G2M — search for any file containing "Whitfield"
    whitfield = None
    for p in REPO.rglob("*"):
        if p.is_file() and "whitfield" in p.name.lower():
            try:
                if p.suffix.lower() in (".gmt",):
                    sets = load_gmt(p)
                    for k, v in sets.items():
                        if "G2M" in k.upper() or "G2_M" in k.upper():
                            whitfield = v
                            log_entries["whitfield_source"] = str(p)
                            log_entries["whitfield_set_name"] = k
                            break
                else:
                    df = pd.read_csv(p, sep="\t", header=None)
                    whitfield = df.iloc[:, 0].astype(str).str.upper().tolist()
                    log_entries["whitfield_source"] = str(p)
            except Exception:
                continue
            if whitfield is not None:
                break
    if whitfield is None:
        log_entries["whitfield_available"] = False
        log_entries["whitfield_note"] = "Whitfield_G2M not found; falling back to HALLMARK_G2M_CHECKPOINT only (3 + 1 = 4 component sets)."
        print("[WARN] Whitfield_G2M unavailable; falling back to 4 Hallmark sets.")
    else:
        log_entries["whitfield_available"] = True
        log_entries["whitfield_n_genes"] = len(whitfield)

    alias = load_alias(REPO / "external_data" / "hgnc_alias_table.tsv")

    for cohort in COHORTS:
        cohort_log: dict = {}
        try:
            axes = compute_axes_for_cohort(cohort, hallmark, whitfield, alias, cohort_log)
            out_path = out_dir / f"axes_{cohort}.tsv"
            axes.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
            cohort_log[cohort]["output"] = str(out_path)
            cohort_log[cohort]["status"] = "ok"
            print(f"[{cohort}] axes -> {out_path}  (n={len(axes)})")
        except Exception as e:
            cohort_log[cohort] = {"status": "fail", "error": str(e)}
            print(f"[{cohort}] FAIL: {e}", file=sys.stderr)
        log_entries["cohorts"].update(cohort_log)

    log_entries["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    with open(log_path, "w") as f:
        json.dump(log_entries, f, indent=2)
    print(f"Log written: {log_path}")


if __name__ == "__main__":
    main()
