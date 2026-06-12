"""Phase 4 / Critical Gaps — Task A.

GSE32646: 2-gene MKI67-ESR1 signed z-score "oracle" vs framework Dox signature.

Convention (matches Phase 1 baseline_expansion.py):
  oracle = z(MKI67) - z(ESR1)
  (signed sum where ESR1 carries negative sign — high prolif, low ER = positive)

Steps:
  1. Load GSE32646 expression
  2. Extract MKI67 and ESR1 (canonical; HGNC alias fallback)
  3. Compute oracle
  4. AUROC oracle vs pCR + bootstrap 95% CI (1000x, seed 42)
  5. Re-extract framework AUROC + CI from Phase 2 score_DOXORUBICIN_GSE32646.tsv
  6. Paired DeLong test (DeLong-Sun-Xu 2014 / Sun-Xu fast 2014) — proper
  7. Output JSON + TSV + headline

Outputs:
  results/phase_4/critical_gaps/oracle_gse32646.json
  results/phase_4/critical_gaps/oracle_gse32646.tsv
"""
from __future__ import annotations
import datetime as _dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

REPO = Path(".")
PHASE2 = REPO / "." / "results" / "phase_2"
OUT_DIR = REPO / "." / "results" / "phase_4" / "critical_gaps"
HGNC_ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_gene(symbol: str, columns: list[str], alias_map: dict) -> str | None:
    cols_upper = {c.upper(): c for c in columns}
    if symbol.upper() in cols_upper:
        return cols_upper[symbol.upper()]
    # Try alias
    if symbol.upper() in alias_map:
        for alt in alias_map[symbol.upper()]:
            if alt.upper() in cols_upper:
                return cols_upper[alt.upper()]
    return None


def load_alias_map() -> dict:
    """Return: alias_symbol_upper -> [canonical_symbols]."""
    if not HGNC_ALIAS.exists():
        return {}
    df = pd.read_csv(HGNC_ALIAS, sep="\t", dtype=str).fillna("")
    out: dict[str, list[str]] = {}
    # Heuristic: assume columns include 'symbol' and 'alias_symbol' / 'prev_symbol'
    sym_col = None
    for cand in ("symbol", "Approved symbol", "approved_symbol"):
        if cand in df.columns:
            sym_col = cand
            break
    if sym_col is None:
        return {}
    alias_cols = [c for c in df.columns if "alias" in c.lower() or "prev" in c.lower() or "previous" in c.lower()]
    for _, row in df.iterrows():
        canonical = row[sym_col].strip().upper()
        if not canonical:
            continue
        # canonical -> canonical (self)
        out.setdefault(canonical, []).append(canonical)
        # any aliases listed map to canonical
        for ac in alias_cols:
            field = row[ac]
            for piece in str(field).replace("|", ",").split(","):
                p = piece.strip().upper()
                if p:
                    out.setdefault(p, []).append(canonical)
    return out


def boot_auroc(y: np.ndarray, s: np.ndarray, n_boot: int = 1000, seed: int = 42):
    rng = np.random.default_rng(seed)
    if len(np.unique(y)) < 2:
        return float("nan"), float("nan"), float("nan")
    obs = float(roc_auc_score(y, s))
    boot = []
    n = len(y)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(y[idx])) < 2:
            continue
        boot.append(roc_auc_score(y[idx], s[idx]))
    if len(boot) < 10:
        return obs, float("nan"), float("nan")
    return obs, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


# -------- DeLong (Sun & Xu 2014 fast implementation) --------
def _midrank(x: np.ndarray) -> np.ndarray:
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=float)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1  # 1-based midrank
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def delong_var_cov(scores_list: list[np.ndarray], y: np.ndarray):
    """Return AUC array (k,), covariance matrix (k,k) for k classifiers on same y."""
    pos = (y == 1)
    neg = ~pos
    m = int(pos.sum())
    n = int(neg.sum())
    k = len(scores_list)
    V10 = np.zeros((k, m), dtype=float)
    V01 = np.zeros((k, n), dtype=float)
    aucs = np.zeros(k, dtype=float)
    for r, sc in enumerate(scores_list):
        X = sc[pos].astype(float)
        Y = sc[neg].astype(float)
        TX = _midrank(X)
        TY = _midrank(Y)
        TZ = _midrank(np.concatenate([X, Y]))
        TZX = TZ[:m]
        TZY = TZ[m:]
        aucs[r] = (TZX.sum() / m - (m + 1) / 2.0) / n
        V10[r, :] = 1.0 - (TZX - TX) / n
        V01[r, :] = (TZY - TY) / m
    S10 = np.cov(V10, ddof=1)
    S01 = np.cov(V01, ddof=1)
    if k == 1:
        S10 = np.array([[np.var(V10[0], ddof=1)]])
        S01 = np.array([[np.var(V01[0], ddof=1)]])
    S = S10 / m + S01 / n
    return aucs, S


def delong_paired_pvalue(y: np.ndarray, s_a: np.ndarray, s_b: np.ndarray):
    """Two-sided paired DeLong test. Returns (AUC_a, AUC_b, delta, z, p_two)."""
    aucs, S = delong_var_cov([s_a, s_b], y)
    delta = aucs[0] - aucs[1]
    var_delta = float(S[0, 0] + S[1, 1] - 2.0 * S[0, 1])
    if var_delta <= 0 or not np.isfinite(var_delta):
        return float(aucs[0]), float(aucs[1]), float(delta), float("nan"), float("nan")
    z = delta / np.sqrt(var_delta)
    p_two = 2.0 * (1.0 - norm.cdf(abs(z)))
    return float(aucs[0]), float(aucs[1]), float(delta), float(z), float(p_two)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    log = {
        "script": str(me),
        "script_sha256": sha256_of(me),
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "convention": {
            "oracle_formula": "z(MKI67) - z(ESR1)  # signed sum, high prolif + low ER = positive",
            "framework_score": "sig_dox column from Phase 2 score_DOXORUBICIN_GSE32646.tsv",
            "delong": "Sun & Xu 2014 fast algorithm (paired two-sided)",
            "bootstrap": "1000x, seed 42, percentile 95% CI",
        },
    }

    # --- Load expression ---
    expr_path = PHASE2 / "GSE32646_expression.tsv"
    pcr_path = PHASE2 / "GSE32646_pcr.tsv"
    fw_path = PHASE2 / "score_DOXORUBICIN_GSE32646.tsv"
    log["expression_sha256"] = sha256_of(expr_path)
    log["pcr_sha256"] = sha256_of(pcr_path)
    log["framework_score_sha256"] = sha256_of(fw_path)

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    pcr = pd.read_csv(pcr_path, sep="\t")
    fw = pd.read_csv(fw_path, sep="\t")

    alias_map = load_alias_map()
    log["alias_table_loaded"] = bool(alias_map)
    log["alias_table_size"] = len(alias_map)

    mki67_col = resolve_gene("MKI67", list(expr.columns), alias_map)
    esr1_col = resolve_gene("ESR1", list(expr.columns), alias_map)
    if mki67_col is None or esr1_col is None:
        raise RuntimeError(f"MKI67 or ESR1 not resolvable; got {mki67_col=} {esr1_col=}")
    log["MKI67_column_used"] = mki67_col
    log["ESR1_column_used"] = esr1_col

    # --- Align patients across expr + pcr + fw ---
    common = (set(expr.index.astype(str)) &
              set(pcr["patient_id"].astype(str)) &
              set(fw["patient_id"].astype(str)))
    common = sorted(common)
    log["n_common"] = len(common)

    expr_c = expr.loc[common]
    pcr_c = pcr.set_index("patient_id").loc[common]["pcr"].astype(int).values
    framework_c = fw.set_index("patient_id").loc[common]["sig_dox"].astype(float).values

    # --- Build 2-gene oracle (z-score uses cohort common-sample mean/sd) ---
    mki = expr_c[mki67_col].astype(float)
    esr = expr_c[esr1_col].astype(float)
    zmki = (mki - mki.mean()) / (mki.std(ddof=0) + 1e-12)
    zesr = (esr - esr.mean()) / (esr.std(ddof=0) + 1e-12)
    oracle = (zmki - zesr).values  # signed sum: high MKI67 - high ESR1

    # Mask non-finite
    mask = np.isfinite(oracle) & np.isfinite(framework_c) & np.isfinite(pcr_c)
    y = pcr_c[mask]
    o = oracle[mask]
    f = framework_c[mask]
    log["n_used"] = int(mask.sum())
    log["n_pos"] = int((y == 1).sum())
    log["n_neg"] = int((y == 0).sum())

    # --- AUROC + bootstrap CI for each ---
    auc_o, lo_o, hi_o = boot_auroc(y, o, 1000, 42)
    auc_f, lo_f, hi_f = boot_auroc(y, f, 1000, 42)

    # --- Bootstrap delta CI ---
    rng = np.random.default_rng(42)
    diffs = []
    for _ in range(1000):
        idx = rng.integers(0, len(y), size=len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        diffs.append(roc_auc_score(y[idx], o[idx]) - roc_auc_score(y[idx], f[idx]))
    diffs = np.array(diffs)
    delta_obs = auc_o - auc_f
    delta_lo, delta_hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))

    # --- Paired DeLong (proper) ---
    auc_o_dl, auc_f_dl, delta_dl, z_dl, p_dl = delong_paired_pvalue(y, o, f)

    verdict = (
        "2-gene > Framework" if (delta_obs > 0 and p_dl < 0.05) else
        "2-gene matches Framework (NS)" if p_dl >= 0.05 else
        "2-gene < Framework"
    )

    result = {
        "cohort": "GSE32646",
        "n": int(mask.sum()),
        "n_pos": int((y == 1).sum()),
        "n_neg": int((y == 0).sum()),
        "two_gene": {
            "formula": "z(MKI67) - z(ESR1)",
            "MKI67_column": mki67_col,
            "ESR1_column": esr1_col,
            "AUROC": auc_o,
            "AUROC_ci_95": [lo_o, hi_o],
        },
        "framework": {
            "source": "Phase 2 score_DOXORUBICIN_GSE32646.tsv [sig_dox]",
            "AUROC": auc_f,
            "AUROC_ci_95": [lo_f, hi_f],
        },
        "delta_two_gene_minus_framework": {
            "value": float(delta_obs),
            "ci_95_bootstrap": [delta_lo, delta_hi],
        },
        "paired_delong": {
            "AUROC_two_gene": auc_o_dl,
            "AUROC_framework": auc_f_dl,
            "delta": delta_dl,
            "z": z_dl,
            "p_value_two_sided": p_dl,
        },
        "verdict": verdict,
    }

    log["result"] = result
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    # --- Write JSON ---
    with open(OUT_DIR / "oracle_gse32646.json", "w") as fh:
        json.dump(log, fh, indent=2)

    # --- Write TSV (1 row) ---
    tsv = pd.DataFrame([{
        "cohort": "GSE32646",
        "n": int(mask.sum()),
        "two_gene_AUROC": auc_o,
        "two_gene_CI_low": lo_o,
        "two_gene_CI_high": hi_o,
        "framework_AUROC": auc_f,
        "framework_CI_low": lo_f,
        "framework_CI_high": hi_f,
        "delta_2gene_minus_framework": delta_obs,
        "delta_CI_low": delta_lo,
        "delta_CI_high": delta_hi,
        "paired_DeLong_z": z_dl,
        "paired_DeLong_p": p_dl,
        "verdict": verdict,
    }])
    tsv.to_csv(OUT_DIR / "oracle_gse32646.tsv", sep="\t", index=False)

    # Headline print
    print("[TASK A] GSE32646:")
    print(f"  2-gene AUROC = {auc_o:.4f} [{lo_o:.4f}, {hi_o:.4f}]")
    print(f"  Framework AUROC = {auc_f:.4f} [{lo_f:.4f}, {hi_f:.4f}]")
    print(f"  Delta = {delta_obs:+.4f} [{delta_lo:+.4f}, {delta_hi:+.4f}]")
    print(f"  Paired DeLong p = {p_dl:.4g} (z = {z_dl:.3f})")
    print(f"  Verdict: {verdict}")


if __name__ == "__main__":
    main()
