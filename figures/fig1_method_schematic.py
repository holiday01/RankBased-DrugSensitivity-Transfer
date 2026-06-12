"""
Fig 1 — Method Schematic (programmatic placeholder)

A matplotlib flowchart of the axis-decomposition method pipeline.
For final submission, a polished BioRender/Illustrator version is preferred.

Usage: python fig1_method_schematic.py
Output: figures/fig1_method_schematic.png + .pdf
"""

import os
import hashlib
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def script_sha256():
    with open(__file__, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def make_box(ax, x, y, w, h, text, fc="white", ec="black", lw=1.5, fontsize=8, fontweight="normal"):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.3,rounding_size=0.3",
        edgecolor=ec,
        facecolor=fc,
        linewidth=lw,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=fontweight,
        wrap=True,
    )


def make_arrow(ax, x1, y1, x2, y2, lw=1.5):
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="->",
        mutation_scale=15,
        linewidth=lw,
        color="black",
    )
    ax.add_patch(arrow)


def main(out_dir):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    # Title
    ax.text(
        50,
        96,
        "Pre-specified axis-decomposition method for cohort-free pharmacogenomic signatures",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

    # Stage 1: inputs (left)
    make_box(ax, 2, 65, 18, 15, "CCLE / GDSC\ndrug-AUC training\n(60-gene sig)", fc="#F5F5F5")
    make_box(ax, 2, 40, 18, 15, "Patient expression\n(GSE32646, METABRIC)\nquantile-normalised", fc="#F5F5F5")

    # Stage 2: deployment
    make_box(ax, 30, 52, 18, 16, "singscore\n(Foroutan 2018)\n→ risk score", fc="#E8F4F8")
    make_arrow(ax, 20, 73, 30, 60)
    make_arrow(ax, 20, 47, 30, 60)

    # Stage 3: axis decomposition (THE METHOD) — taller box, higher position
    make_box(
        ax,
        55,
        38,
        25,
        38,
        "AXIS-DECOMPOSITION\n(the method)\n\nCompute 4 axes:\nP (Hallmark prolif)\nH (ESR1+PGR)\nB (ERBB2)\nBasal_prob (PAM50)\n\nFit linear OLS:\nsig ~ P+H+B+Basal\n\nResidual AUROC\n+ bootstrap 95% CI",
        fc="#A23B72",
        ec="black",
        lw=2,
        fontsize=7,
        fontweight="bold",
    )
    make_arrow(ax, 48, 60, 55, 60)

    # Stage 4: decision rule — lower position, no overlap
    make_box(ax, 35, 11, 14, 18, "AXIS-EQUIVALENT\nresid ≤ 0.55\nCI ≤ 0.55", fc="#A8D5BA", lw=1.5, fontsize=7)
    make_box(ax, 52, 11, 14, 18, "GRAY ZONE\n0.55 < pt ≤ 0.60\nOR CI crosses", fc="#FFE699", lw=1.5)
    make_box(ax, 69, 11, 14, 18, "RESIDUAL SIGNAL\nresid > 0.60\nCI > 0.55", fc="#F4A4A4", lw=1.5, fontsize=7)
    # decision arrows from bottom of axis-decomp box (y=38) to top of decision boxes (y=29)
    make_arrow(ax, 60, 38, 42, 29)
    make_arrow(ax, 67, 38, 59, 29)
    make_arrow(ax, 75, 38, 76, 29)

    # Stage 5: GATE companion
    make_box(
        ax,
        85,
        45,
        13,
        30,
        "GATE\nchecklist\n\nWithin-stratum\nAUROC\n(GAP test)\n\nfloor ≈ 0.60",
        fc="#E8F4F8",
        ec="black",
        lw=1.5,
        fontsize=7,
    )
    make_arrow(ax, 80, 60, 85, 60)

    # Footer
    ax.text(
        50,
        4,
        "Pre-specified (exploratory-with-pre-specified-falsifier): protocol.LOCK + protocol sha256; Zenodo deposit pending — see Methods §2.10",
        ha="center",
        va="center",
        fontsize=7,
        style="italic",
        color="#555555",
    )

    out_png = os.path.join(out_dir, "fig1_method_schematic.png")
    out_pdf = os.path.join(out_dir, "fig1_method_schematic.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"sha256: {script_sha256()}")
    print(f"wrote: {out_png}")
    print(f"wrote: {out_pdf}")


if __name__ == "__main__":
    out_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    main(out_dir)
