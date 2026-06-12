#!/usr/bin/env python3
"""Fig 6: GSEA Hallmark NES heatmap (drug x hallmark).

Loads results/phase_1/hallmark_gsea_NES_matrix.tsv + FDR matrix.
Selects top-15 Hallmark sets by mean |NES| across drugs.
Heatmap rows=drugs, cols=hallmark, color=NES (RdBu_r), annotated with
FDR-significance asterisks (* < 0.05, ** < 0.01, *** < 0.001).

Outputs:
    figures/fig6_gsea_nes_heatmap.png + .pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

ROOT = Path(".")
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

NES_TSV = ROOT / "results/phase_1/hallmark_gsea_NES_matrix.tsv"
FDR_TSV = ROOT / "results/phase_1/hallmark_gsea_FDR_matrix.tsv"


def main() -> int:
    if not NES_TSV.exists():
        print(f"[ERR] missing {NES_TSV}")
        return 1
    nes = pd.read_csv(NES_TSV, sep="\t", index_col=0)
    fdr = (pd.read_csv(FDR_TSV, sep="\t", index_col=0)
           if FDR_TSV.exists() else None)
    print(f"[INFO] NES shape={nes.shape}")

    # Coerce numerics; empty cells -> NaN
    nes = nes.apply(pd.to_numeric, errors="coerce")
    if fdr is not None:
        fdr = fdr.apply(pd.to_numeric, errors="coerce")

    # ---------- select top-15 hallmarks ----------
    mean_abs = nes.abs().mean(axis=0, skipna=True)
    top15 = mean_abs.sort_values(ascending=False).head(15).index.tolist()
    nes_sel = nes[top15].copy()
    if fdr is not None:
        fdr_sel = fdr[top15].copy()
    else:
        fdr_sel = None

    # cluster cols (top15 hallmarks) and rows (drugs) by 1-cos / corr-like
    # use simple agglomerative on abs(NES) correlation
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    def reorder(df, axis):
        X = df.fillna(0.0).to_numpy()
        if axis == 1:
            X = X.T
        if X.shape[0] < 2:
            return df
        d = pdist(X, metric="correlation")
        d = np.nan_to_num(d, nan=1.0)
        Z = linkage(d, method="average")
        order = leaves_list(Z)
        if axis == 0:
            return df.iloc[order, :]
        else:
            return df.iloc[:, order]

    nes_sel = reorder(nes_sel, axis=0)
    nes_sel = reorder(nes_sel, axis=1)
    if fdr_sel is not None:
        fdr_sel = fdr_sel.loc[nes_sel.index, nes_sel.columns]

    # short labels
    def short_h(s: str) -> str:
        return s.replace("HALLMARK_", "").replace("_", " ").title()

    col_labels = [short_h(c) for c in nes_sel.columns]
    row_labels = list(nes_sel.index)

    # ---------- figure ----------
    fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
    vmax = max(np.nanmax(np.abs(nes_sel.to_numpy())), 1e-6)
    im = ax.imshow(nes_sel.to_numpy(), cmap="RdBu_r",
                   vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Hallmark GSEA NES across drug signatures (top-15 |NES|)",
                 fontsize=11)

    # annotate NES + FDR asterisks
    n_rows, n_cols = nes_sel.shape
    for i in range(n_rows):
        for j in range(n_cols):
            v = nes_sel.iat[i, j]
            if pd.isna(v):
                continue
            color = "white" if abs(v) > vmax * 0.55 else "black"
            stars = ""
            if fdr_sel is not None:
                q = fdr_sel.iat[i, j]
                if pd.notna(q):
                    if q < 1e-3:
                        stars = "***"
                    elif q < 1e-2:
                        stars = "**"
                    elif q < 5e-2:
                        stars = "*"
            ax.text(j, i, f"{v:+.1f}{stars}", ha="center", va="center",
                    fontsize=5.5, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("NES", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # FDR significance legend
    ax.text(1.02, -0.18,
            "FDR:  * < 0.05    ** < 0.01    *** < 0.001",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            color="0.3")

    png = FIG_DIR / "fig6_gsea_nes_heatmap.png"
    pdf = FIG_DIR / "fig6_gsea_nes_heatmap.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
