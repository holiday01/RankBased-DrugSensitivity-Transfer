#!/usr/bin/env python3
"""Fig 2: Axis-decomposition scatter (GSE32646 + METABRIC).

For each held-out cohort, fit linear model
    sig ~ P + H + B + H_PAM50basal
and scatter framework signature score vs fitted axis-projection score.

Outputs:
    figures/fig2_axis_decomp_scatter.png
    figures/fig2_axis_decomp_scatter.pdf
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score

# ---------- style ----------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "lines.linewidth": 1.5,
})

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/resubmission_v2")
FIG_DIR = __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/manuscript")
FIG_DIR.mkdir(parents=True, exist_ok=True)

COHORTS = [
    {
        "name": "GSE32646",
        "axes": ROOT / "results/phase_2/axes_GSE32646.tsv",
        "score": __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/analysis/corrected_figdata/phase_2/score_DOXORUBICIN_GSE32646.tsv"),
        "verdict_label": "NO-MARGINAL-SIGNAL",
        "verdict_color": "#1f77b4",
        "id_col": "patient_id",
        "outcome": "pcr",
        "label_neg": "Non-pCR",
        "label_pos": "pCR",
        "label_y": "5-yr pCR endpoint",
        "label_title": "GSE32646 (neoadjuvant P-FEC, n=115)",
    },
    {
        "name": "METABRIC",
        "axes": ROOT / "results/phase_4/metabric/axes_METABRIC.tsv",
        "score": __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/analysis/corrected_figdata/phase_4/metabric/score_DOXORUBICIN_METABRIC.tsv"),
        "verdict_label": "orthogonal endpoint,\nlow-power consistency",
        "verdict_color": "#7f7f7f",
        "id_col": "sample_id",
        "outcome": "event5y",
        "label_neg": "No 5-yr event",
        "label_pos": "5-yr OS event",
        "label_y": "5-yr OS event",
        "label_title": "METABRIC (chemo-treated, n=399)",
    },
]


def load_cohort(cfg: dict) -> pd.DataFrame | None:
    if not cfg["axes"].exists() or not cfg["score"].exists():
        print(f"[WARN] {cfg['name']}: input files missing -- skipping")
        return None
    ax = pd.read_csv(cfg["axes"], sep="\t")
    sc = pd.read_csv(cfg["score"], sep="\t")
    # rename ID col
    ax = ax.rename(columns={ax.columns[0]: "id"})
    sc = sc.rename(columns={sc.columns[0]: "id"})
    df = ax.merge(sc, on="id", how="inner")
    # rename outcome col(s) -- drop duplicate outcome from axes if present
    cols = list(df.columns)
    # ensure we have P, H, B, H_PAM50basal, sig_dox, outcome
    need = ["P", "H", "B", "H_PAM50basal", "sig_dox"]
    if not all(c in df.columns for c in need):
        print(f"[WARN] {cfg['name']}: missing required columns {need} -- have {df.columns.tolist()}")
        return None
    # outcome col: pcr or event5y might be in both; collapse
    outc = cfg["outcome"]
    outc_cols = [c for c in df.columns if c == outc or c.startswith(outc + "_")]
    if outc not in df.columns and outc_cols:
        df[outc] = df[outc_cols[0]]
    if outc not in df.columns:
        print(f"[WARN] {cfg['name']}: outcome '{outc}' missing -- skipping")
        return None
    return df


def fit_axes(df: pd.DataFrame) -> tuple[np.ndarray, float, dict]:
    X = df[["P", "H", "B", "H_PAM50basal"]].to_numpy()
    y = df["sig_dox"].to_numpy()
    mdl = LinearRegression().fit(X, y)
    yhat = mdl.predict(X)
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    coefs = {n: float(c) for n, c in zip(["P", "H", "B", "H_PAM50basal"], mdl.coef_)}
    coefs["intercept"] = float(mdl.intercept_)
    return yhat, r2, coefs


def main() -> int:
    panels = []
    for cfg in COHORTS:
        out = load_cohort(cfg)
        if out is None:
            panels.append(None)
            continue
        df = out
        yhat, r2, _ = fit_axes(df)
        panels.append((cfg, df, yhat, r2))

    if all(p is None for p in panels):
        print("[ERR] no cohorts loaded -- nothing to draw")
        return 1

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.2), constrained_layout=True)
    colors = plt.cm.tab10.colors

    for ax, panel in zip(axes, panels):
        if panel is None:
            ax.text(0.5, 0.5, "Data missing", ha="center", va="center")
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        cfg, df, yhat, r2 = panel
        y = df["sig_dox"].to_numpy()
        outcome = df[cfg["outcome"]].to_numpy()

        # color points by outcome
        pos = outcome == 1
        ax.scatter(y[~pos], yhat[~pos], s=14, c=[colors[0]], alpha=0.55,
                   edgecolors="none", label=f"{cfg['label_neg']} (n={int((~pos).sum())})")
        ax.scatter(y[pos], yhat[pos], s=20, c=[colors[3]], alpha=0.75,
                   edgecolors="none", marker="^",
                   label=f"{cfg['label_pos']} (n={int(pos.sum())})")

        lo = min(y.min(), yhat.min())
        hi = max(y.max(), yhat.max())
        pad = (hi - lo) * 0.05
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                ls="--", color="0.4", lw=1.0, label="y = x")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)

        ax.set_xlabel("Framework signature score")
        # y-axis label only on the leftmost panel to avoid centre overlap
        if ax is axes[0]:
            ax.set_ylabel("Axis-projection fit\nP+H+B+Basal")
        ax.set_title(cfg["label_title"])

        # Annotations are computed from the same rows this panel plots, so the figure and
        # its caption cannot drift apart. The previous version read them from a JSON path
        # that did not exist, which is why every panel printed "nan".
        raw_point = roc_auc_score(outcome, y)
        residual_point = roc_auc_score(outcome, y - yhat)
        vshort = cfg["verdict_label"]
        vcolor = cfg["verdict_color"]

        ann = (
            f"$R^2$ = {r2:.3f}\n"
            f"raw AUROC = {raw_point:.3f}\n"
            f"residual AUROC = {residual_point:.3f}\n"
            f"{vshort}"
        )
        ax.text(0.04, 0.96, ann, transform=ax.transAxes,
                va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=vcolor, lw=1.0, alpha=0.9))
        ax.legend(loc="lower right", framealpha=0.9, fontsize=7)
        ax.grid(True, linestyle=":", alpha=0.4)

    fig.suptitle("Framework signature carries no pCR signal beyond the P+H+B+Basal axes",
                 fontsize=12, fontweight="bold")

    png = FIG_DIR / "fig2.png"
    pdf = FIG_DIR / "fig2.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
