"""Phase 4 / METABRIC — second held-out cohort confirmation.

Pipeline (single script, mirrors Phase 2 GSE32646 logic):
  M2  extract chemo subset with 5-year OS endpoint
  M3  compute axes (P, H, H_ESR1only, H_PAM50basal, B) — uses METABRIC's
      built-in CLAUDIN_SUBTYPE for PAM50 (Basal indicator → H_PAM50basal)
  M4  apply frozen DOXORUBICIN signature via quantile-norm + singscore
  M5  compute 2-gene MKI67-ESR1 oracle (z(MKI67) - z(ESR1))
  M6  residual confirmatory test:
        sig_dox ~ β0 + β1·P + β2·H + β3·B + β4·H_PAM50basal + ε
        → residual AUROC vs 5y-death, 1000× bootstrap CI, seed 42
        3×3 verdict (point ≤ 0.55 AND ci_lo ≤ 0.55 → H-b CONFIRMED)
      Supplementary variants (Bonferroni × 3):
        - H_ESR1only in place of H
        - H_PAM50basal as H replacement
        - spline residual (pygam) — skipped if pygam absent

Outputs to results/phase_4/metabric/:
  axes_METABRIC.tsv
  score_DOXORUBICIN_METABRIC.tsv
  oracle_METABRIC.tsv
  metabric_residual_test.json
  metabric_residual_summary.tsv
  PHASE_4_METABRIC_LOG.md
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

REPO = Path(".")
PHASE_4 = REPO / "." / "results" / "phase_4" / "metabric"
SIG_PATH = REPO / "frozen_signatures" / "solid_only" / "quantile_30_oncokb_all" / "DOXORUBICIN.tsv"
HALLMARK_PATH = REPO / "h.all.v2025.1.Hs.symbols.gmt"
HGNC_ALIAS = REPO / "external_data" / "hgnc_alias_table.tsv"

HALLMARK_SETS = [
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_MYC_TARGETS_V1",
]


# ===== utilities =====================================================
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
            out[parts[0]] = [g.strip().upper() for g in parts[2:] if g.strip()]
    return out


def load_alias_simple(path: Path) -> dict:
    """alias_symbol_upper -> canonical_symbol_upper."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "alias_symbol" in df.columns and "canonical_symbol" in df.columns:
        return dict(zip(df["alias_symbol"].str.upper(), df["canonical_symbol"].str.upper()))
    return {}


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


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Rows = samples, cols = genes. Standard per-sample QN to a shared dist."""
    arr = df.values.astype(float)
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
    return sens, res


def singscore_signed(expr: pd.DataFrame, up: set, down: set):
    expr.columns = [c.upper() for c in expr.columns]
    n_genes_total = expr.shape[1]
    if n_genes_total < 5000:
        raise ValueError(f"Only {n_genes_total} genes; refuse to score.")
    rank_matrix = expr.rank(axis=1, method="average", pct=True)
    up_present = sorted(up & set(expr.columns))
    down_present = sorted(down & set(expr.columns))
    if len(up_present) < 3 or len(down_present) < 3:
        raise ValueError(f"Too few sig genes: up={len(up_present)}, down={len(down_present)}")

    def normalize(score_mean, n_set):
        mn = (1 + n_set) / (2 * n_genes_total)
        mx = 1 - mn
        return 2 * (score_mean - mn) / (mx - mn) - 1

    up_score = normalize(rank_matrix[up_present].mean(axis=1), len(up_present))
    down_score = normalize((1 - rank_matrix[down_present]).mean(axis=1), len(down_present))
    raw = up_score - down_score
    centered = raw - raw.median()
    return centered, {
        "n_genes_in_expression": n_genes_total,
        "n_up_in_sig": len(up),
        "n_up_present": len(up_present),
        "n_down_in_sig": len(down),
        "n_down_present": len(down_present),
    }


def auroc_with_boot_ci(scores, labels, n_boot=1000, seed=42):
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    if len(np.unique(labels)) < 2:
        return float("nan"), float("nan"), float("nan"), []
    obs = float(roc_auc_score(labels, scores))
    rng = np.random.default_rng(seed)
    n = len(labels)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        if len(np.unique(labels[idx])) < 2:
            continue
        boot.append(roc_auc_score(labels[idx], scores[idx]))
    if len(boot) < 10:
        return obs, float("nan"), float("nan"), boot
    return obs, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)), boot


def verdict_3x3(point: float, ci_lo: float):
    if not (math.isfinite(point) and math.isfinite(ci_lo)):
        return "INDETERMINATE_NAN"
    if point <= 0.55 and ci_lo <= 0.55:
        return "H-b CONFIRMED (no residual signal after deconfounding)"
    if point > 0.60 and ci_lo > 0.55:
        return "H-b FALSIFIED (residual signal survives)"
    return "GRAY ZONE"


# ===== DeLong (Sun & Xu 2014) =====================================
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
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, dtype=float)
    T2[J] = T
    return T2


def delong_paired_pvalue(y: np.ndarray, s_a: np.ndarray, s_b: np.ndarray):
    pos = (y == 1)
    neg = ~pos
    m = int(pos.sum())
    n = int(neg.sum())
    if m == 0 or n == 0:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
    scores_list = [s_a, s_b]
    k = 2
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
    S = S10 / m + S01 / n
    delta = float(aucs[0] - aucs[1])
    var_delta = float(S[0, 0] + S[1, 1] - 2.0 * S[0, 1])
    if var_delta <= 0 or not np.isfinite(var_delta):
        return float(aucs[0]), float(aucs[1]), delta, float("nan"), float("nan")
    z = delta / np.sqrt(var_delta)
    p_two = 2.0 * (1.0 - norm.cdf(abs(z)))
    return float(aucs[0]), float(aucs[1]), delta, float(z), float(p_two)


# ===== fit residual ===============================================
def fit_residual(df: pd.DataFrame, h_col: str = "H", include_pam50: bool = True):
    if include_pam50 and h_col != "H_PAM50basal":
        cols = ["P", h_col, "B", "H_PAM50basal"]
    else:
        cols = ["P", h_col, "B"]
    needed = ["sig_dox", "event5y"] + cols
    sub = df.dropna(subset=needed).copy()
    if len(sub) < 20:
        return None, sub, None
    X = sub[cols].values
    y = sub["sig_dox"].values
    lr = LinearRegression()
    lr.fit(X, y)
    sub["residual"] = y - lr.predict(X)
    coef = {c: float(v) for c, v in zip(cols, lr.coef_)}
    coef["intercept"] = float(lr.intercept_)
    return sub, sub, coef


# ===== main pipeline ===============================================
def main():
    PHASE_4.mkdir(parents=True, exist_ok=True)
    log = {
        "cohort": "METABRIC",
        "started_utc": _dt.datetime.utcnow().isoformat() + "Z",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_of(Path(__file__).resolve()),
        "endpoint": "5-year all-cause death (OS_STATUS=1:DECEASED AND OS_MONTHS<=60); censored before 60mo dropped",
        "chemo_filter": "CHEMOTHERAPY == 'YES'",
        "decision_rule": {
            "confirmed": "point ≤ 0.55 AND ci_lo ≤ 0.55",
            "gray_zone": "0.55 < point ≤ 0.60 OR ci crosses 0.55",
            "falsified": "point > 0.60 AND ci_lo > 0.55",
        },
    }
    print("[M] Loading inputs...", flush=True)
    clin_path = PHASE_4 / "METABRIC_clinical.tsv"
    expr_path = PHASE_4 / "METABRIC_expression.tsv"
    log["clinical_sha256"] = sha256_of(clin_path)
    log["expression_sha256"] = sha256_of(expr_path)
    log["signature_sha256"] = sha256_of(SIG_PATH)
    log["hallmark_sha256"] = sha256_of(HALLMARK_PATH)

    clin = pd.read_csv(clin_path, sep="\t", low_memory=False)
    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.columns = [c.upper() for c in expr.columns]
    log["expression_shape_raw"] = list(expr.shape)
    log["clinical_shape_raw"] = list(clin.shape)

    # ===== M2: chemo subset + 5-year endpoint =====
    chemo = clin[clin["CHEMOTHERAPY"] == "YES"].copy()
    chemo["sample_id"] = chemo["SAMPLE_ID"].astype(str)

    def to_event5y(row):
        deceased = str(row["OS_STATUS"]).startswith("1")
        try:
            months = float(row["OS_MONTHS"])
        except Exception:
            return np.nan
        if not np.isfinite(months):
            return np.nan
        if deceased and months <= 60:
            return 1
        if months >= 60:
            return 0
        return np.nan

    chemo["event5y"] = chemo.apply(to_event5y, axis=1)
    log["n_chemo_total"] = int(len(chemo))
    log["n_chemo_with_event5y"] = int(chemo["event5y"].notna().sum())
    log["n_chemo_event_pos"] = int((chemo["event5y"] == 1).sum())
    log["n_chemo_event_neg"] = int((chemo["event5y"] == 0).sum())
    log["chemo_event5y_rate"] = float((chemo["event5y"] == 1).mean())

    chemo_samples = sorted(set(chemo["sample_id"]) & set(expr.index.astype(str)))
    log["n_chemo_with_expr"] = len(chemo_samples)
    expr_chemo = expr.loc[chemo_samples].copy()
    chemo_with_event = chemo[chemo["sample_id"].isin(chemo_samples) & chemo["event5y"].notna()].copy()
    chemo_with_event = chemo_with_event.set_index("sample_id")
    log["n_chemo_with_expr_and_event5y"] = int(len(chemo_with_event))

    print(f"[M2] chemo subset: total={log['n_chemo_total']}, with event5y={log['n_chemo_with_event5y']}, "
          f"in expr={log['n_chemo_with_expr']}, both={log['n_chemo_with_expr_and_event5y']}", flush=True)

    # ===== M3: axes (P, H, H_ESR1only, H_PAM50basal, B) =====
    alias = load_alias_simple(HGNC_ALIAS)
    hallmark = load_gmt(HALLMARK_PATH)

    # Use the chemo-only expression matrix for ranking (consistent with cohort definition)
    rank_pct = expr_chemo.rank(axis=1, method="average", pct=True)

    expr_cols_set = set(expr_chemo.columns)
    p_set_diag = {}
    comp_z = {}
    for name in HALLMARK_SETS:
        gs = hallmark.get(name, [])
        resolved = sorted({resolve_gene(g, expr_cols_set, alias) for g in gs} - {None})
        score = singscore_unsigned(rank_pct, resolved)
        comp_z[name] = zscore(score)
        p_set_diag[name] = {"n_genes_in_set": len(gs), "n_resolved_present": len(resolved)}
    log["p_set_diagnostics"] = p_set_diag
    log["whitfield_available"] = False
    P = pd.DataFrame(comp_z).mean(axis=1)

    esr1 = resolve_gene("ESR1", expr_cols_set, alias)
    pgr = resolve_gene("PGR", expr_cols_set, alias)
    erbb2 = resolve_gene("ERBB2", expr_cols_set, alias)
    mki67 = resolve_gene("MKI67", expr_cols_set, alias)
    log.update({"ESR1": esr1, "PGR": pgr, "ERBB2": erbb2, "MKI67": mki67})

    H_ESR1_z = zscore(expr_chemo[esr1]) if esr1 else pd.Series(np.nan, index=expr_chemo.index)
    H_PGR_z = zscore(expr_chemo[pgr]) if pgr else pd.Series(np.nan, index=expr_chemo.index)
    H = pd.concat([H_ESR1_z, H_PGR_z], axis=1).mean(axis=1) if (esr1 and pgr) else H_ESR1_z
    B = zscore(expr_chemo[erbb2]) if erbb2 else pd.Series(np.nan, index=expr_chemo.index)

    # H_PAM50basal from CLAUDIN_SUBTYPE (METABRIC built-in PAM50+claudin call)
    subtype_map = chemo.set_index("sample_id")["CLAUDIN_SUBTYPE"].to_dict()
    H_PAM50basal_raw = pd.Series(
        {s: (1.0 if subtype_map.get(s) == "Basal" else 0.0)
         for s in expr_chemo.index},
        index=expr_chemo.index,
    )
    # z-score the indicator so it lives on the same scale as other axes
    H_PAM50basal = zscore(H_PAM50basal_raw)
    log["pam50_source"] = "METABRIC CLAUDIN_SUBTYPE (built-in); H_PAM50basal = z(I[subtype==Basal])"
    log["pam50_subtype_counts"] = chemo["CLAUDIN_SUBTYPE"].value_counts(dropna=False).to_dict()

    event_map = chemo.set_index("sample_id")["event5y"].to_dict()
    axes = pd.DataFrame({
        "sample_id": expr_chemo.index,
        "P": P.values,
        "H": H.values,
        "H_ESR1only": H_ESR1_z.values,
        "H_PAM50basal": H_PAM50basal.values,
        "B": B.values,
        "event5y": [event_map.get(s, np.nan) for s in expr_chemo.index],
    })
    axes_path = PHASE_4 / "axes_METABRIC.tsv"
    axes.to_csv(axes_path, sep="\t", index=False, float_format="%.6g")
    log["axes_output_sha256"] = sha256_of(axes_path)

    # ===== M4: frozen DOXORUBICIN signature with QN + singscore =====
    print("[M4] quantile-normalizing chemo expression matrix...", flush=True)
    expr_qn = quantile_normalize(expr_chemo)
    sens, res = load_signature(SIG_PATH)
    scores, diag = singscore_signed(expr_qn, sens, res)
    log["dox_sig_diagnostics"] = diag

    score_df = pd.DataFrame({
        "sample_id": expr_chemo.index,
        "sig_dox": scores.values,
        "event5y": [event_map.get(s, np.nan) for s in expr_chemo.index],
    })
    score_path = PHASE_4 / "score_DOXORUBICIN_METABRIC.tsv"
    score_df.to_csv(score_path, sep="\t", index=False, float_format="%.10g")
    log["score_output_sha256"] = sha256_of(score_path)

    # ===== M5: 2-gene MKI67-ESR1 oracle =====
    if mki67 is None or esr1 is None:
        raise RuntimeError(f"MKI67={mki67}, ESR1={esr1} not resolvable")
    zmki = zscore(expr_chemo[mki67])
    zesr = zscore(expr_chemo[esr1])
    oracle = zmki - zesr
    oracle_df = pd.DataFrame({
        "sample_id": expr_chemo.index,
        "z_MKI67": zmki.values,
        "z_ESR1": zesr.values,
        "oracle_2gene": oracle.values,
        "event5y": [event_map.get(s, np.nan) for s in expr_chemo.index],
    })
    oracle_path = PHASE_4 / "oracle_METABRIC.tsv"
    oracle_df.to_csv(oracle_path, sep="\t", index=False, float_format="%.6g")
    log["oracle_output_sha256"] = sha256_of(oracle_path)

    # ===== M6: residual confirmatory test =====
    df = axes.merge(score_df[["sample_id", "sig_dox"]], on="sample_id", how="inner")
    df = df.merge(oracle_df[["sample_id", "oracle_2gene"]], on="sample_id", how="inner")
    df["event5y"] = pd.to_numeric(df["event5y"], errors="coerce")
    log["n_merged"] = int(len(df))

    sub, _, coef = fit_residual(df, h_col="H", include_pam50=True)
    if sub is None:
        raise RuntimeError("Too few rows for primary residual model")
    y = sub["event5y"].astype(int).to_numpy()
    r = sub["residual"].to_numpy()
    auc, lo, hi, _ = auroc_with_boot_ci(r, y, 1000, 42)
    primary_verdict = verdict_3x3(auc, lo)
    log["primary_model"] = {
        "formula": "sig_dox ~ P + H + B + H_PAM50basal",
        "n_used": int(len(sub)),
        "coefficients": coef,
    }
    log["primary_residual_auroc"] = {
        "point": auc, "ci_low": lo, "ci_high": hi, "n_bootstrap": 1000, "seed": 42,
    }
    log["primary_verdict"] = primary_verdict

    # Raw sig_dox AUROC
    raw_auc, raw_lo, raw_hi, _ = auroc_with_boot_ci(sub["sig_dox"].to_numpy(), y, 1000, 42)
    log["raw_sig_dox_auroc"] = {"point": raw_auc, "ci_low": raw_lo, "ci_high": raw_hi}

    # 2-gene oracle vs framework
    auc_o, lo_o, hi_o, _ = auroc_with_boot_ci(sub["oracle_2gene"].to_numpy(), y, 1000, 42)
    auc_f, lo_f, hi_f, _ = auroc_with_boot_ci(sub["sig_dox"].to_numpy(), y, 1000, 42)
    rng = np.random.default_rng(42)
    diffs = []
    o_arr = sub["oracle_2gene"].to_numpy()
    f_arr = sub["sig_dox"].to_numpy()
    for _ in range(1000):
        idx = rng.integers(0, len(y), size=len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        diffs.append(roc_auc_score(y[idx], o_arr[idx]) - roc_auc_score(y[idx], f_arr[idx]))
    diffs = np.array(diffs)
    delta_obs = auc_o - auc_f
    delta_lo, delta_hi = float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))
    auc_o_dl, auc_f_dl, delta_dl, z_dl, p_dl = delong_paired_pvalue(y, o_arr, f_arr)
    log["two_gene_vs_framework"] = {
        "n": int(len(sub)),
        "two_gene_AUROC": auc_o,
        "two_gene_CI_95": [lo_o, hi_o],
        "framework_AUROC": auc_f,
        "framework_CI_95": [lo_f, hi_f],
        "delta_two_gene_minus_framework": delta_obs,
        "delta_CI_95_bootstrap": [delta_lo, delta_hi],
        "paired_DeLong_z": z_dl,
        "paired_DeLong_p": p_dl,
    }

    # ===== supplementary variants (Bonf×3) =====
    supp = {}
    for variant_name, h_arg, include_pam in [
        ("H_ESR1only", "H_ESR1only", True),
        ("H_PAM50basal_only", "H_PAM50basal", True),
    ]:
        sub_v, _, coef_v = fit_residual(df, h_col=h_arg, include_pam50=include_pam)
        if sub_v is None or len(sub_v) < 20:
            supp[variant_name] = {"status": "skipped_n<20"}
            continue
        y_v = sub_v["event5y"].astype(int).to_numpy()
        r_v = sub_v["residual"].to_numpy()
        av, lv, hv, _ = auroc_with_boot_ci(r_v, y_v, 1000, 42)
        boot_v = []
        rng = np.random.default_rng(42)
        for _ in range(2000):
            idx = rng.integers(0, len(y_v), size=len(y_v))
            if len(np.unique(y_v[idx])) < 2:
                continue
            boot_v.append(roc_auc_score(y_v[idx], r_v[idx]))
        alpha = 0.05 / 3
        bf_lo = float(np.percentile(boot_v, 100 * alpha / 2)) if boot_v else float("nan")
        bf_hi = float(np.percentile(boot_v, 100 * (1 - alpha / 2))) if boot_v else float("nan")
        supp[variant_name] = {
            "n_used": int(len(sub_v)),
            "point": av,
            "ci_low_95": lv,
            "ci_high_95": hv,
            "ci_low_bonf3": bf_lo,
            "ci_high_bonf3": bf_hi,
            "verdict_unadjusted": verdict_3x3(av, lv),
            "verdict_bonf3": verdict_3x3(av, bf_lo),
        }

    # spline residual
    try:
        from pygam import LinearGAM, s as gam_s
        sub_s = df.dropna(subset=["sig_dox", "event5y", "P", "H", "B"]).copy()
        if len(sub_s) >= 20:
            Xs = sub_s[["P", "H", "B"]].values
            ys = sub_s["sig_dox"].values
            gam = LinearGAM(gam_s(0) + gam_s(1) + gam_s(2)).fit(Xs, ys)
            sub_s["residual_spline"] = ys - gam.predict(Xs)
            y_s = sub_s["event5y"].astype(int).to_numpy()
            r_s = sub_s["residual_spline"].to_numpy()
            asp, lsp, hsp, _ = auroc_with_boot_ci(r_s, y_s, 1000, 42)
            boot_s = []
            rng = np.random.default_rng(42)
            for _ in range(2000):
                idx = rng.integers(0, len(y_s), size=len(y_s))
                if len(np.unique(y_s[idx])) < 2:
                    continue
                boot_s.append(roc_auc_score(y_s[idx], r_s[idx]))
            alpha = 0.05 / 3
            bf_lo = float(np.percentile(boot_s, 100 * alpha / 2)) if boot_s else float("nan")
            bf_hi = float(np.percentile(boot_s, 100 * (1 - alpha / 2))) if boot_s else float("nan")
            supp["spline_residual"] = {
                "n_used": int(len(sub_s)),
                "point": asp,
                "ci_low_95": lsp,
                "ci_high_95": hsp,
                "ci_low_bonf3": bf_lo,
                "ci_high_bonf3": bf_hi,
                "verdict_unadjusted": verdict_3x3(asp, lsp),
                "verdict_bonf3": verdict_3x3(asp, bf_lo),
            }
    except ImportError:
        supp["spline_residual"] = {"status": "skipped", "note": "pygam not installed"}
    except Exception as e:
        supp["spline_residual"] = {"status": "failed", "error": str(e)}

    log["supplementary_variants"] = supp
    log["finished_utc"] = _dt.datetime.utcnow().isoformat() + "Z"

    out_json = PHASE_4 / "metabric_residual_test.json"
    with open(out_json, "w") as f:
        json.dump(log, f, indent=2, default=str)

    # Summary TSV
    rows = []
    rows.append({"test": "primary (H + H_PAM50basal)", "n": log["primary_model"]["n_used"],
                 "auroc_point": auc, "ci_low_95": lo, "ci_high_95": hi, "verdict": primary_verdict})
    rows.append({"test": "raw sig_dox (no deconfounding)", "n": log["primary_model"]["n_used"],
                 "auroc_point": raw_auc, "ci_low_95": raw_lo, "ci_high_95": raw_hi,
                 "verdict": "informational only"})
    rows.append({"test": "2-gene oracle (MKI67-ESR1)", "n": int(len(sub)),
                 "auroc_point": auc_o, "ci_low_95": lo_o, "ci_high_95": hi_o,
                 "verdict": f"delta vs framework = {delta_obs:+.4f} [95%CI {delta_lo:+.4f}, {delta_hi:+.4f}], DeLong p={p_dl:.4g}"})
    for name, v in supp.items():
        if "point" in v:
            rows.append({"test": f"supp: {name}", "n": v["n_used"],
                         "auroc_point": v["point"], "ci_low_95": v["ci_low_95"], "ci_high_95": v["ci_high_95"],
                         "ci_low_bonf3": v.get("ci_low_bonf3"), "ci_high_bonf3": v.get("ci_high_bonf3"),
                         "verdict": v.get("verdict_bonf3")})
        else:
            rows.append({"test": f"supp: {name}", "n": None, "auroc_point": None, "verdict": v.get("status")})
    summary = pd.DataFrame(rows)
    summary.to_csv(PHASE_4 / "metabric_residual_summary.tsv", sep="\t", index=False, float_format="%.4f")

    # Headline print
    print(f"\n===== METABRIC HEADLINE (n={log['primary_model']['n_used']}) =====")
    print(f"Primary residual AUROC = {auc:.3f}  [95% CI: {lo:.3f}, {hi:.3f}]")
    print(f"Verdict: {primary_verdict}")
    print(f"Raw sig_dox AUROC = {raw_auc:.3f} [{raw_lo:.3f}, {raw_hi:.3f}]")
    print(f"2-gene oracle AUROC = {auc_o:.3f} [{lo_o:.3f}, {hi_o:.3f}]")
    print(f"Framework AUROC = {auc_f:.3f} [{lo_f:.3f}, {hi_f:.3f}]")
    print(f"Delta = {delta_obs:+.4f} [{delta_lo:+.4f}, {delta_hi:+.4f}], DeLong p={p_dl:.4g}")
    print("\nSupplementary variants:")
    for name, v in supp.items():
        if "point" in v:
            print(f"  {name}: AUROC={v['point']:.3f} 95%CI=[{v['ci_low_95']:.3f},{v['ci_high_95']:.3f}] "
                  f"Bonf3=[{v['ci_low_bonf3']:.3f},{v['ci_high_bonf3']:.3f}] verdict={v['verdict_bonf3']}")
        else:
            print(f"  {name}: {v}")


if __name__ == "__main__":
    main()
