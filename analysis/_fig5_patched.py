#!/usr/bin/env python3
"""Fig 5: 9-signature cross-domain verdict composite.

Compiles per-signature residual JSONs (Phase 6/7/8) into a 3-row composite:
  Top: paired bar chart of marginal AUROC + residual AUROC (with CI)
  Middle: 9 x 4 axis-correlation heatmap (axis coefficient values, RdBu_r)
  Right: vertical verdict column (CONFIRMED/GRAY/FALSIFIED color squares)

Outputs:
    figures/fig5_cross_signature_verdict.png + .pdf
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/analysis/corrected_figdata")
_OUTROOT = __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/manuscript")
FIG_DIR = _OUTROOT
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Order: GGI, RS-like, Dox-variant, Mammaprint, TIS, cytolytic, DLDA30, TGFB, HRD
ENTRIES = [
    ("GGI",        ROOT / "phase_6/GGI_gse25066_residual.json",        "GSE25066"),
    ("RS-like",    ROOT / "phase_6/RSlike_gse25066_HRpos_residual.json", "GSE25066 HR+"),
    ("Dox-variant", ROOT / "phase_6/DoxVariant_gse25066_residual.json",  "GSE25066"),
    ("Mammaprint", ROOT / "phase_6/Mammaprint_metabric_residual.json",  "METABRIC"),
    ("TIS",        ROOT / "phase_7/TIS_gse25066_residual.json",         "GSE25066"),
    ("Cytolytic",  ROOT / "phase_7/cytolytic_gse22226_residual.json",   "GSE22226"),
    ("DLDA30",     ROOT / "phase_7/DLDA30_gse16446_residual.json",      "GSE16446"),
    ("TGFB-strom", ROOT / "phase_8/TGFBstrom_gse25066_residual.json",   "GSE25066"),
    ("HR/Fanconi", ROOT / "phase_8/HRD_gse22226_residual.json",         "GSE22226"),
]

AXES_ORDER = ["P", "H", "B", "H_PAM50basal"]
AXES_LABELS = ["P", "H", "B", r"H$_{PAM50basal}$"]


_A2 = None
def _table3_class(name: str, cohort: str):
    """Table 3's pre-specified verdict class, so the figure and the table cannot disagree.

    The deposited `verdict_3x3` predates the C1 correction and fires CONFIRMED on cells with
    no marginal signal at all -- it labels TGFBstrom "axis-equivalent" on a marginal AUROC
    whose interval lies entirely below 0.5. Table 3 is the correct classification.
    """
    global _A2
    if _A2 is None:
        import json as _j
        _A2 = {(c["signature"], c["cohort"]): c["verdict"]
               for c in _j.load(open("/home/holiday01/2026_ISMB_code/revise_bioadv_386/"
                                     "results/A2/A2_grid.json"))["cells"]}
    alias = {"HR/Fanconi": "HRD", "TGFB-strom": "TGFBstrom", "RS-like": "RSlike",
             "Dox-variant": "DoxVariant", "Cytolytic": "cytolytic", "Mammaprint": "MammaPrint"}
    v = _A2.get((alias.get(name, name), cohort))
    return {"NO-MARGINAL-SIGNAL": ("N", "#7f7f7f"),
            "REDUNDANT-WITH-AXES": ("R", "#1f77b4"),
            "INDETERMINATE": ("I", "#bcbd22"),
            "INVERSE-DIRECTION": ("X", "#d62728"),
            "NOT-EVALUABLE": ("--", "#bbbbbb")}.get(v, ("n/a", "#dddddd"))


def verdict_short(v: str) -> tuple[str, str]:
    if v is None:
        return ("?", "#bbbbbb")
    s = v.upper()
    if "CONFIRMED" in s:
        return ("AXIS-EQUIV", "#1f77b4")
    if "GRAY" in s or "GREY" in s:
        return ("GRAY", "#bcbd22")
    if "FALSIFIED" in s:
        return ("RESIDUAL", "#d62728")
    return (v[:10], "#bbbbbb")


def main() -> int:
    rows = []
    for name, jpath, cohort in ENTRIES:
        if not jpath.exists():
            print(f"[WARN] missing {jpath} -- skipping {name}")
            continue
        d = json.loads(jpath.read_text())
        marg = d.get("marginal_auroc", {})
        resd = d.get("residual_auroc", {})
        coefs = d.get("residual_model", {}).get("coefficients", {})
        v = d.get("verdict_3x3", "")
        rows.append({
            "name": name,
            "cohort": cohort,
            "marg_pt": marg.get("point", np.nan),
            "marg_lo": marg.get("ci_low", np.nan),
            "marg_hi": marg.get("ci_high", np.nan),
            "res_pt": resd.get("point", np.nan),
            "res_lo": resd.get("ci_low", np.nan),
            "res_hi": resd.get("ci_high", np.nan),
            "coefs": {a: coefs.get(a, np.nan) for a in AXES_ORDER},
            "verdict": v,
        })

    if not rows:
        print("[ERR] no signatures loaded")
        return 1

    n = len(rows)

    # ---------- figure layout ----------
    fig = plt.figure(figsize=(7.0, 5.0), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2, ncols=3,
        height_ratios=[1.0, 1.05],
        width_ratios=[1.0, 0.10, 0.04],
        hspace=0.45, wspace=0.04,
        left=0.12, right=0.95, top=0.92, bottom=0.10,
    )
    ax_bar = fig.add_subplot(gs[0, 0:2])
    ax_heat = fig.add_subplot(gs[1, 0])
    ax_verdict = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])

    # ===== TOP: paired bar chart of marginal vs residual AUROC =====
    x = np.arange(n)
    width = 0.38
    marg_pt = np.array([r["marg_pt"] for r in rows])
    res_pt = np.array([r["res_pt"] for r in rows])
    marg_err = np.vstack([marg_pt - [r["marg_lo"] for r in rows],
                          [r["marg_hi"] for r in rows] - marg_pt])
    res_err = np.vstack([res_pt - [r["res_lo"] for r in rows],
                         [r["res_hi"] for r in rows] - res_pt])

    ax_bar.bar(x - width / 2, marg_pt, width, yerr=marg_err,
               capsize=2, color="#1f77b4", label="Marginal AUROC",
               edgecolor="white", lw=0.5,
               error_kw=dict(elinewidth=0.8, ecolor="0.25"))
    ax_bar.bar(x + width / 2, res_pt, width, yerr=res_err,
               capsize=2, color="#ff7f0e", label="Residual AUROC (after P+H+B+Basal)",
               edgecolor="white", lw=0.5,
               error_kw=dict(elinewidth=0.8, ecolor="0.25"))
    ax_bar.axhline(0.55, color="0.4", ls="--", lw=1, alpha=0.7)
    ax_bar.axhline(0.60, color="0.4", ls=":", lw=1, alpha=0.7)
    ax_bar.text(n - 0.5, 0.555, "0.55 (axis-equiv cap)", fontsize=6,
                ha="right", color="0.4")
    ax_bar.text(n - 0.5, 0.605, "0.60 (residual-signal floor)", fontsize=6,
                ha="right", color="0.4")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([f"{r['name']}\n({r['cohort']})" for r in rows],
                           rotation=30, ha="right", fontsize=7.5)
    ax_bar.set_ylabel("AUROC (95% CI)")
    ax_bar.set_ylim(0.30, 0.85)
    ax_bar.legend(loc="upper right", fontsize=7, frameon=True)
    ax_bar.grid(True, axis="y", linestyle=":", alpha=0.4)
    ax_bar.set_title(
        "Nine published breast-cancer signatures: marginal vs residual AUROC",
        fontsize=11)

    # ===== MIDDLE: axis-coefficient heatmap (n x 4) =====
    M = np.array([[r["coefs"][a] for a in AXES_ORDER] for r in rows])
    vmax = max(np.nanmax(np.abs(M)), 1e-6)
    im = ax_heat.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        aspect="auto")
    ax_heat.set_xticks(range(len(AXES_ORDER)))
    ax_heat.set_xticklabels(AXES_LABELS, fontsize=8)
    ax_heat.set_yticks(range(n))
    ax_heat.set_yticklabels([r["name"] for r in rows], fontsize=8)
    ax_heat.set_xlabel("Axis coefficient in sig ~ P + H + B + H$_{PAM50basal}$")

    for i in range(n):
        for j in range(len(AXES_ORDER)):
            v = M[i, j]
            color = "white" if abs(v) > vmax * 0.6 else "black"
            ax_heat.text(j, i, f"{v:+.2f}", ha="center", va="center",
                         fontsize=6.5, color=color)

    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("Axis coefficient", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # ===== VERDICT COLUMN =====
    ax_verdict.set_xlim(0, 1)
    ax_verdict.set_ylim(n - 0.5, -0.5)
    # short verdict glyph: C / G / F so it never overflows
    glyph = {"AXIS-EQUIV": "AE", "GRAY": "IND", "RESIDUAL": "RES"}
    for i, r in enumerate(rows):
        short, color = _table3_class(r["name"], r["cohort"])
        ax_verdict.add_patch(plt.Rectangle((0.05, i - 0.35), 0.9, 0.7,
                                           color=color, ec="black", lw=0.4))
        ax_verdict.text(0.5, i, short,
                        ha="center", va="center",
                        fontsize=7.5, color="white", fontweight="bold")
    ax_verdict.set_xticks([])
    ax_verdict.set_yticks([])
    ax_verdict.set_title("Table 3\nclass", fontsize=8)
    for spine in ax_verdict.spines.values():
        spine.set_visible(False)
    # legend strip below the verdict col (outside axes)
    fig.text(0.97, 0.04,
             "N = no marginal signal   R = redundant with axes   I = indeterminate\nX = inverse direction   n/a = outside Family A (Mammaprint is a survival-endpoint cell)",
             ha="right", va="bottom", fontsize=7, color="0.25")

    png = FIG_DIR / "fig5.png"
    pdf = FIG_DIR / "fig5.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
