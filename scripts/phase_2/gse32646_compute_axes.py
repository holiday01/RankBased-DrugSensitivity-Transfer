"""Phase 2 / B2 — Compute axes (P, H, H_ESR1only, H_PAM50basal, B) on GSE32646
plus PAM50 calls. Mirrors Phase 0 compute_axes.py logic but reads inputs from
phase_2 outputs (not external_data/<cohort>/).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(".")
PHASE_2 = REPO / "." / "results" / "phase_2"
COHORT = "GSE32646"
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
    return dict(
        zip(df["alias_symbol"].astype(str).str.upper(), df["canonical_symbol"].astype(str).str.upper())
    )


def resolve_gene(symbol: str, expr_cols: set, alias: dict):
    s = symbol.upper()
    if s in expr_cols:
        return s
    a = alias.get(s)
    if a and a in expr_cols:
        return a
    return None


def singscore_unsigned(rank_pct: pd.DataFrame, gene_set: list) -> pd.Series:
    n_genes_total = rank_pct.shape[1]
    present = [g for g in gene_set if g in rank_pct.columns]
    if len(present) < 3:
        return pd.Series(np.nan, index=rank_pct.index)
    n_set = len(present)
    mn = (1 + n_set) / (2 * n_genes_total)
    mx = 1 - mn
    raw = rank_pct[present].mean(axis=1)
    return 2 * (raw - mn) / (mx - mn) - 1


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=s.index)
    return (s - mu) / sd


def compute_pam50_inline(expr: pd.DataFrame, log: dict) -> pd.DataFrame:
    """Compute PAM50 directly here (avoids re-shelling out)."""
    centroid_path = REPO / "external_data" / "pam50_centroids.tsv"
    cent = pd.read_csv(centroid_path, sep="\t", index_col=0)
    cent.index = [g.upper() for g in cent.index]
    cent = cent[["Basal", "Her2", "LumA", "LumB", "Normal"]]
    log["pam50_centroid_sha256"] = sha256_of(centroid_path)

    # alias resolution
    alias_path = REPO / "external_data" / "hgnc_alias_table.tsv"
    alias_map = {}
    if alias_path.exists():
        adf = pd.read_csv(alias_path, sep="\t")
        for _, row in adf.iterrows():
            approved = str(row.get("symbol", "")).upper().strip() if "symbol" in adf.columns else ""
            if not approved and "canonical_symbol" in adf.columns:
                approved = str(row.get("canonical_symbol", "")).upper().strip()
            if not approved:
                continue
            for col in ["alias_symbol", "prev_symbol"]:
                if col not in adf.columns:
                    continue
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
    resolved = {}
    for g in missing:
        if g in alias_map and alias_map[g] in cohort_genes:
            resolved[g] = alias_map[g]
        else:
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
    log["pam50_n_mapped"] = n_mapped
    fallback = n_mapped < 40

    if fallback:
        log["pam50_fallback_triggered"] = True
        out = pd.DataFrame({"patient_id": expr.index})
        out["pam50_subtype"] = "NA_FALLBACK"
        for s in ["Basal", "Her2", "LumA", "LumB", "Normal"]:
            out[f"{s.lower()}_score"] = np.nan
        return out

    expr_sub = expr[available_pam50].copy()
    expr_sub.columns = [rename_pam50.get(c, c) for c in expr_sub.columns]
    expr_sub = expr_sub.subtract(expr_sub.median(axis=0), axis=1)
    cent_sub = cent.loc[expr_sub.columns]

    results = pd.DataFrame(
        index=expr_sub.index,
        columns=["Basal", "Her2", "LumA", "LumB", "Normal"],
        dtype=float,
    )
    for pid in expr_sub.index:
        pvec = expr_sub.loc[pid].values
        for s in results.columns:
            cvec = cent_sub[s].values
            rho, _ = spearmanr(pvec, cvec)
            results.loc[pid, s] = rho

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
    log["pam50_subtype_counts"] = out["pam50_subtype"].value_counts().to_dict()
    return out


def main():
    PHASE_2.mkdir(parents=True, exist_ok=True)
    log: dict = {
        "cohort": COHORT,
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
    }

    expr_path = PHASE_2 / "GSE32646_expression.tsv"
    pcr_path = PHASE_2 / "GSE32646_pcr.tsv"
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    log["expression_path"] = str(expr_path)
    log["expression_sha256"] = sha256_of(expr_path)
    log["n_samples_expr"] = int(expr.shape[0])
    log["n_genes_expr"] = int(expr.shape[1])

    hallmark_path = REPO / "h.all.v2025.1.Hs.symbols.gmt"
    hallmark = load_gmt(hallmark_path)
    log["hallmark_sha256"] = sha256_of(hallmark_path)

    alias = load_alias(REPO / "external_data" / "hgnc_alias_table.tsv")

    # --- PAM50 (compute first; H_PAM50basal needs it) -------------------------
    pam50 = compute_pam50_inline(expr, log)
    pam50_out = PHASE_2 / "GSE32646_pam50.tsv"
    pam50.to_csv(pam50_out, sep="\t", index=False, float_format="%.6f")
    log["pam50_output"] = str(pam50_out)

    # --- Per-sample percentile ranks for singscore ----------------------------
    rank_pct = expr.rank(axis=1, method="average", pct=True)

    # --- P axis components ----------------------------------------------------
    comp_z_cols = {}
    p_set_diag = {}
    expr_cols_set = set(expr.columns)
    for name in HALLMARK_SETS:
        gs = hallmark.get(name, [])
        gs_resolved = []
        missing = []
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
    log["p_set_diagnostics"] = p_set_diag

    # Whitfield_G2M (optional fallback)
    whitfield = None
    for p in REPO.rglob("*"):
        if p.is_file() and "whitfield" in p.name.lower():
            try:
                if p.suffix.lower() == ".gmt":
                    sets = load_gmt(p)
                    for k, v in sets.items():
                        if "G2M" in k.upper() or "G2_M" in k.upper():
                            whitfield = v
                            log["whitfield_source"] = str(p)
                            log["whitfield_set_name"] = k
                            break
                else:
                    df = pd.read_csv(p, sep="\t", header=None)
                    whitfield = df.iloc[:, 0].astype(str).str.upper().tolist()
                    log["whitfield_source"] = str(p)
            except Exception:
                continue
            if whitfield is not None:
                break
    if whitfield is not None:
        gs_resolved = sorted({resolve_gene(g, expr_cols_set, alias) for g in whitfield} - {None})
        score = singscore_unsigned(rank_pct, gs_resolved)
        comp_z_cols["Whitfield_G2M"] = zscore(score)
        log["whitfield_available"] = True
    else:
        log["whitfield_available"] = False
        print("[WARN] Whitfield_G2M unavailable; using 4 Hallmark sets.")

    p_df = pd.DataFrame(comp_z_cols)
    P = p_df.mean(axis=1)

    # --- H axis ---------------------------------------------------------------
    esr1 = resolve_gene("ESR1", expr_cols_set, alias)
    pgr = resolve_gene("PGR", expr_cols_set, alias)
    erbb2 = resolve_gene("ERBB2", expr_cols_set, alias)
    log["ESR1_symbol"] = esr1
    log["PGR_symbol"] = pgr
    log["ERBB2_symbol"] = erbb2

    H_ESR1_z = zscore(expr[esr1]) if esr1 else pd.Series(np.nan, index=expr.index)
    H_PGR_z = zscore(expr[pgr]) if pgr else pd.Series(np.nan, index=expr.index)
    if esr1 and pgr:
        H = pd.concat([H_ESR1_z, H_PGR_z], axis=1).mean(axis=1)
    elif esr1:
        H = H_ESR1_z
    else:
        H = pd.Series(np.nan, index=expr.index)

    # --- H_PAM50basal ---------------------------------------------------------
    pam50_indexed = pam50.set_index("patient_id")
    H_PAM50basal = pam50_indexed["basal_score"].reindex(expr.index)

    # --- B axis ---------------------------------------------------------------
    B = zscore(expr[erbb2]) if erbb2 else pd.Series(np.nan, index=expr.index)

    # --- pCR ------------------------------------------------------------------
    pcr_df = pd.read_csv(pcr_path, sep="\t")
    pcr_df["patient_id"] = pcr_df["patient_id"].astype(str)
    pcr_df = pcr_df.set_index("patient_id")
    pCR = pd.to_numeric(pcr_df["pcr"], errors="coerce").reindex(expr.index)

    out = pd.DataFrame(
        {
            "sample_id": expr.index,
            "P": P.values,
            "H": H.values,
            "H_ESR1only": H_ESR1_z.values,
            "H_PAM50basal": H_PAM50basal.values,
            "B": B.values,
            "pCR": pCR.values,
        }
    )
    out_path = PHASE_2 / "axes_GSE32646.tsv"
    out.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    log["axes_output"] = str(out_path)
    log["n_samples_with_pCR"] = int(out["pCR"].notna().sum())
    log["n_samples_with_P"] = int(out["P"].notna().sum())
    log["n_samples_with_H"] = int(out["H"].notna().sum())
    log["n_samples_with_B"] = int(out["B"].notna().sum())
    log["n_samples_with_PAM50basal"] = int(out["H_PAM50basal"].notna().sum())

    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"
    log_out = PHASE_2 / "compute_axes_GSE32646.log"
    with open(log_out, "w") as f:
        json.dump(log, f, indent=2)
    print(f"[GSE32646] axes -> {out_path}  (n={len(out)})")
    print(f"  PAM50 subtype counts: {log.get('pam50_subtype_counts')}")
    print(f"  axes nonNaN: P={log['n_samples_with_P']}, H={log['n_samples_with_H']}, "
          f"B={log['n_samples_with_B']}, PAM50basal={log['n_samples_with_PAM50basal']}, pCR={log['n_samples_with_pCR']}")


if __name__ == "__main__":
    main()
