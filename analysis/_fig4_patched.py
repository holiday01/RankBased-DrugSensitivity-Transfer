#!/usr/bin/env python3
"""Fig 4: Hierarchically clustered inter-signature correlation heatmap (GSE25066).

Loads the pre-computed correlation matrix from Phase 1, clusters with
average-linkage on (1-|r|) distance, draws the clustermap with class-coded
row colors.

Signature classes:
    cytotoxic: DOXORUBICIN, ETOPOSIDE, PACLITAXEL, CARBOPLATIN,
               CYCLOPHOSPHAMIDE, FLUOROURACIL, TOP2_consensus
    TKI:       TRAMETINIB, GEFITINIB, LAPATINIB, IMATINIB
    H-positive: FULVESTRANT
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/resubmission_v2")
FIG_DIR = __import__("pathlib").Path("/tmp/figout")
FIG_DIR.mkdir(parents=True, exist_ok=True)

CORR_TSV = __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/analysis/corrected_figdata/phase_1/correlation_matrix_GSE25066.tsv")

CLASS_MAP = {
    "DOXORUBICIN": "cytotoxic",
    "ETOPOSIDE": "cytotoxic",
    "PACLITAXEL": "cytotoxic",
    "CARBOPLATIN": "cytotoxic",
    "CYCLOPHOSPHAMIDE": "cytotoxic",
    "FLUOROURACIL": "cytotoxic",
    "TOP2_consensus": "cytotoxic",
    "TRAMETINIB": "TKI",
    "GEFITINIB": "TKI",
    "LAPATINIB": "TKI",
    "IMATINIB": "TKI",
    "FULVESTRANT": "H-positive",
}

CLASS_COLOR = {
    "cytotoxic": "#d62728",   # red
    "TKI": "#1f77b4",         # blue
    "H-positive": "#2ca02c",  # green
    "null": "#7f7f7f",        # gray
}


def main() -> int:
    if not CORR_TSV.exists():
        print(f"[ERR] missing {CORR_TSV}")
        return 1
    corr = pd.read_csv(CORR_TSV, sep="\t", index_col=0)
    if corr.shape[0] != corr.shape[1]:
        print(f"[ERR] non-square: {corr.shape}")
        return 1
    sigs = corr.columns.tolist()
    print(f"[INFO] loaded {corr.shape[0]}x{corr.shape[1]} matrix; sigs={sigs}")

    # average-linkage on (1-|r|) distance
    dist = 1 - np.abs(corr.values)
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    order_idx = leaves_list(Z)
    ordered_sigs = [sigs[i] for i in order_idx]
    corr_ord = corr.loc[ordered_sigs, ordered_sigs]

    # ---------- figure ----------
    fig = plt.figure(figsize=(9.5, 7.5), constrained_layout=False)
    n = len(ordered_sigs)
    # 3-column layout: [rowclass | heat | cbar]
    # left=0.22 gives ~2.1 in left margin so CYCLOPHOSPHAMIDE (16-char) fits
    # without truncation at fontsize 8.
    gs = fig.add_gridspec(
        nrows=2, ncols=3,
        width_ratios=[0.05, 1.0, 0.06],
        height_ratios=[0.18, 1.0],
        wspace=0.03, hspace=0.04,
        left=0.22, right=0.97, top=0.90, bottom=0.28,
    )
    ax_dendro_top = fig.add_subplot(gs[0, 1])
    ax_rowclass   = fig.add_subplot(gs[1, 0])
    ax_heat       = fig.add_subplot(gs[1, 1])
    ax_cbar       = fig.add_subplot(gs[1, 2])

    # top dendrogram
    from scipy.cluster.hierarchy import dendrogram
    ddata = dendrogram(Z, ax=ax_dendro_top, no_labels=True,
                       color_threshold=0, above_threshold_color="0.3")
    # Align dendro x-axis with heatmap columns.
    # scipy places leaf i at x = 10*i + 5; heatmap column i is at data-x = i.
    # Setting xlim(0, n*10) maps dendro leaf i→(10i+5)/range to the same
    # figure fraction as heatmap column i→(i+0.5)/n.
    ax_dendro_top.set_xlim(0, n * 10)
    ax_dendro_top.set_xticks([])
    ax_dendro_top.set_yticks([])
    for spine in ax_dendro_top.spines.values():
        spine.set_visible(False)

    # row class strip
    rowcolors = [CLASS_COLOR[CLASS_MAP.get(s, "null")] for s in ordered_sigs]
    ax_rowclass.imshow(
        np.array([[0] * len(ordered_sigs)]).T,
        aspect="auto",
        cmap=matplotlib.colors.ListedColormap(["#ffffff"]),
    )
    for i, c in enumerate(rowcolors):
        ax_rowclass.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=c))
    ax_rowclass.set_xlim(0, 1)
    ax_rowclass.set_ylim(len(ordered_sigs) - 0.5, -0.5)
    ax_rowclass.set_xticks([])
    ax_rowclass.set_yticks([])
    for spine in ax_rowclass.spines.values():
        spine.set_visible(False)

    # heatmap
    im = ax_heat.imshow(corr_ord.values, cmap="RdBu_r", vmin=-1, vmax=1,
                        aspect="auto")
    ax_heat.set_xticks(range(len(ordered_sigs)))
    ax_heat.set_xticklabels(ordered_sigs, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks(range(len(ordered_sigs)))
    ax_heat.set_yticklabels(ordered_sigs, fontsize=8, ha="right")
    # pad=35 shifts y-tick labels 35 pt left of the heatmap left edge,
    # clearing the row-class colour strip (~23 pt wide) with room to spare.
    ax_heat.tick_params(axis="y", pad=35)
    # color-code tick labels by signature class (so Imatinib reads as TKI not Cytotoxic)
    for tick_label, sig in zip(ax_heat.get_xticklabels(), ordered_sigs):
        tick_label.set_color(CLASS_COLOR[CLASS_MAP.get(sig, "cytotoxic")])
        tick_label.set_fontweight("bold")
    for tick_label, sig in zip(ax_heat.get_yticklabels(), ordered_sigs):
        tick_label.set_color(CLASS_COLOR[CLASS_MAP.get(sig, "cytotoxic")])
        tick_label.set_fontweight("bold")

    # annotate r values (small)
    for i in range(len(ordered_sigs)):
        for j in range(len(ordered_sigs)):
            v = corr_ord.values[i, j]
            color = "white" if abs(v) > 0.55 else "black"
            ax_heat.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=5.5, color=color)

    # colorbar
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Pearson r", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # title
    fig.suptitle("Inter-signature correlations on GSE25066 (n=508)\n"
                 "Average-linkage clustering on 1−|r|",
                 fontsize=11, fontweight="bold", y=0.99)

    # legend for classes
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=CLASS_COLOR["cytotoxic"], edgecolor="0.2",
              label="Cytotoxic"),
        Patch(facecolor=CLASS_COLOR["TKI"], edgecolor="0.2", label="TKI"),
        Patch(facecolor=CLASS_COLOR["H-positive"], edgecolor="0.2",
              label="H-positive"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               fontsize=9, frameon=False, bbox_to_anchor=(0.60, 0.01))

    png = FIG_DIR / "fig4_inter_signature_correlation.png"
    pdf = FIG_DIR / "fig4_inter_signature_correlation.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
