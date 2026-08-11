#!/usr/bin/env python3
"""fig1 regenerated for the revision.

Copied from revise_bioadv/resubmission_v2/figures/scripts/fig1_v3_schematic.py with the
source box corrected: "trained" -> "derived" (R3-A; the signature is a frozen gene list,
not a fitted classifier), "CCLE/GDSC" -> "CTRP/CTD2" (correction C4; GDSC was never the
source), and "Pearson-r" -> the quantile mean-difference rule actually used in section 2.1.
Output redirected to manuscript/.
"""
"""Figure 1 (v3): incremental-value test schematic — nested-model LRT primary,
residual-AUROC/meta/TOST as a concordant robustness grid, spike-in power."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"]})
FIG = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/manuscript")

# Coordinate system: 0–120 (x) × 0–100 (y); figsize 11.0 × 7.5 in
# → ~0.092 in/unit (x), ~0.075 in/unit (y)
# fontsize 7.5 char width ≈ 0.6×7.5/72 = 0.0625 in ≈ 0.68 units
# → 25-char line ≈ 17 units; fits comfortably in a 30-unit-wide box

def box(ax, x, y, w, h, text, *, fc="#EAEFF5", ec="#445", fs=7.5, lw=1.3, zorder=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.8,rounding_size=1.5",
                                fc=fc, ec=ec, lw=lw, zorder=zorder,
                                clip_on=False))
    ax.text(x + w / 2, y + h / 2, text,
            ha="center", va="center", fontsize=fs, zorder=zorder + 1,
            clip_on=False)


def arrow(ax, x1, y1, x2, y2, lw=1.4, color="#333", ls="-"):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
        lw=lw, color=color, ls=ls, shrinkA=3, shrinkB=3, zorder=1,
        clip_on=False))


fig, ax = plt.subplots(figsize=(11.0, 7.5))
ax.set_xlim(0, 120)
ax.set_ylim(0, 100)
ax.axis("off")

# ── Title ─────────────────────────────────────────────────────────────────────
ax.text(60, 98,
        "Does the cohort-free signature add predictive value beyond the four clinical axes?",
        ha="center", va="top", fontsize=11, fontweight="bold", clip_on=False)

# ── Input boxes (left column, x = 1–27) ─────────────────────────────────────
box(ax, 1, 66, 26, 22,
    "Cell-line–derived signature\n(CTRP/CTD² viability AUC, quantile\nmean-difference, OncoKB, top-30 each way)",
    fc="#FBEEDC")

box(ax, 1, 40, 26, 22,
    "Patient expression matrix\n(GSE32646, METABRIC,\nGSE16446 / 25066 / 22226)",
    fc="#FBEEDC")

# ── Singscore box ─────────────────────────────────────────────────────────────
box(ax, 32, 55, 23, 20,
    "Cohort-free singscore\n(full-transcriptome rank)\n→ one score per patient",
    fc="#E3ECF7")
arrow(ax, 27, 77, 32, 72)
arrow(ax, 27, 51, 32, 62)

# ── Clinical axes box ─────────────────────────────────────────────────────────
box(ax, 32, 30, 23, 20,
    "Compute four\nclinical axes\nP · H · B · Basal",
    fc="#E3ECF7")
arrow(ax, 43.5, 55, 43.5, 50)

# ── PRIMARY INFERENCE box ─────────────────────────────────────────────────────
# PRIMARY INFERENCE box: draw outline only (text placed explicitly below)
# box bottom=22, top=77, separator at y=68; header section y=68–77, body y=22–68
ax.add_patch(plt.matplotlib.patches.FancyBboxPatch(
    (59, 22), 38, 55,
    boxstyle="round,pad=0.8,rounding_size=1.5",
    fc="#D7E9D7", ec="#2E7D32", lw=2.0, zorder=2, clip_on=False))

# Header text — centred in the header band (y=68 to y=77 → mid=72.5)
ax.text(78, 72.5, "PRIMARY INFERENCE",
        ha="center", va="center", fontsize=9.5, fontweight="bold",
        color="#1B5E20", zorder=4, clip_on=False)

# separator line
ax.plot([61, 95], [68.0, 68.0], color="#2E7D32", lw=1.0, zorder=4, clip_on=False)

# Body text — centred in the body band (y=22 to y=68 → mid=45)
ax.text(78, 45,
        "Nested-model LRT\n\n"
        "pCR ~ P+H+B+Basal\n"
        "      vs  +signature\n\n"
        "Added value:\n"
        "  LRT p  ·  OR/SD\n"
        "  ΔAUROC  ·  ΔAIC/BIC",
        ha="center", va="center", fontsize=8.2, zorder=3, clip_on=False)

arrow(ax, 55, 37, 59, 44)                          # axes → primary
arrow(ax, 55, 63, 59, 60, ls=(0, (4, 2)), color="#777")  # singscore → primary (dashed)

# ── Verdict boxes ─────────────────────────────────────────────────────────────
box(ax, 56, 2, 20, 17,
    "NO ADDED VALUE\n(equivalence,\npower-bounded)",
    fc="#A8D5BA", ec="#2E7D32", fs=7.5)
box(ax, 78, 2, 18, 17,
    "RESIDUAL\nSIGNAL\nPRESENT",
    fc="#F4A4A4", ec="#B33", fs=7.5)
arrow(ax, 71, 22, 65, 19)
arrow(ax, 83, 22, 87, 19)

# ── Robustness grid (right column) ────────────────────────────────────────────
box(ax, 101, 52, 18, 30,
    "ROBUSTNESS GRID\n(all concordant)\n\n"
    "• residual-AUROC\n"
    "   decomposition\n"
    "   (sig ~ axes)\n"
    "• fixed-effect\n"
    "   meta-analysis\n"
    "• TOST equivalence",
    fc="#EFE7F5", ec="#6A4C93", fs=6.9)
arrow(ax, 97, 58, 101, 62, color="#6A4C93")

# ── Spike-in power box (right lower) ──────────────────────────────────────────
box(ax, 101, 28, 18, 20,
    "Orthogonal spike-in\n→ power floor\nMDE@80% ≈ 0.60–0.75",
    fc="#FDEBD0", ec="#B9770E", fs=7.0)
arrow(ax, 97, 50, 101, 42, color="#B9770E")

fig.savefig(FIG / "fig1.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "fig1.pdf", bbox_inches="tight")
print("wrote fig1_v3_schematic.png/.pdf")
