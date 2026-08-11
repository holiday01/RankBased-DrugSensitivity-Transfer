#!/usr/bin/env python3
"""Fig 3: Forest plot of per-cohort Delta AUROC (2-gene vs framework).

Two-column layout: left = forest plot, right = annotation text panel (no spine).
This avoids annotation text overflowing the plot box.

Outputs:
    figures/fig3_forest_2gene.png
    figures/fig3_forest_2gene.pdf
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
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 9,
    "legend.fontsize": 8,
    "lines.linewidth": 1.5,
})

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv/resubmission_v2")
FIG_DIR = __import__("pathlib").Path("/tmp/figout")
FIG_DIR.mkdir(parents=True, exist_ok=True)

META_JSON = __import__("pathlib").Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/analysis/corrected_figdata/phase_4/critical_gaps/meta_stats.json")


def pretty_p(p: float) -> str:
    if p is None or np.isnan(p):
        return "n/a"
    if p < 1e-9:
        return f"p = {p:.2e}"
    if p < 1e-3:
        return f"p = {p:.2e}"
    return f"p = {p:.3f}"


def main() -> int:
    if not META_JSON.exists():
        print(f"[ERR] missing {META_JSON}")
        return 1
    meta = json.loads(META_JSON.read_text())
    summary = meta["summary"]
    cohort_rows = summary["per_cohort"]
    fe = summary["fixed_effects_meta"]
    stouffer = summary["stouffer_combined"]

    desired = ["GSE16446", "GSE25066", "GSE22226", "GSE32646"]
    rows = []
    by_name = {r["cohort"]: r for r in cohort_rows}
    for n in desired:
        if n in by_name:
            rows.append(by_name[n])
    if not rows:
        print("[ERR] no per-cohort rows")
        return 1

    z975 = 1.959963984540054
    cohorts, deltas, cis_lo, cis_hi, p_vals, n_strs = [], [], [], [], [], []
    for r in rows:
        d = r["delta"]
        se = r["se_delta"]
        cohorts.append(r["cohort"])
        deltas.append(d)
        cis_lo.append(d - z975 * se)
        cis_hi.append(d + z975 * se)
        p_vals.append(r["DeLong_p_two_sided"])
        n_strs.append(f"n={r['n']}, +{r['n_pos']}/-{r['n_neg']}")

    fe_delta = fe["delta"]
    fe_lo, fe_hi = fe["ci_95"]

    # ── Layout: left = forest plot, right = annotation text ──────────────────
    fig = plt.figure(figsize=(9.0, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.1], wspace=0.0,
                          left=0.10, right=0.98, top=0.88, bottom=0.13)
    ax = fig.add_subplot(gs[0])          # forest plot
    ax_txt = fig.add_subplot(gs[1])      # annotation panel (no axis)
    ax_txt.axis("off")

    labels = cohorts + ["Fixed-effects meta"]
    n_total = len(labels)
    y_pos = np.arange(n_total)[::-1]    # 4, 3, 2, 1, 0  (top → bottom)

    # ── Cohort markers ────────────────────────────────────────────────────────
    for i, (d, lo, hi) in enumerate(zip(deltas, cis_lo, cis_hi)):
        ax.errorbar(d, y_pos[i], xerr=[[d - lo], [hi - d]],
                    fmt="s", color="#1f77b4", ecolor="#1f77b4",
                    elinewidth=1.5, capsize=3, markersize=6,
                    markeredgecolor="white", markeredgewidth=0.5)

    # ── Summary diamond ───────────────────────────────────────────────────────
    diam_y = y_pos[-1]
    ax.fill([fe_lo, fe_delta, fe_hi, fe_delta],
            [diam_y, diam_y + 0.32, diam_y, diam_y - 0.32],
            color="#d62728", alpha=0.85, ec="black", lw=1.0, zorder=5)

    # ── Reference line at Δ = 0 ───────────────────────────────────────────────
    ax.axvline(0, color="0.3", ls="--", lw=1.0, zorder=1)

    # ── Axes styling ──────────────────────────────────────────────────────────
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$\Delta$ AUROC  (2-gene $-$ framework)", fontsize=10)

    xmax = max(cis_hi + [fe_hi])
    xmin = min(cis_lo + [fe_lo, 0])
    pad = (xmax - xmin) * 0.08
    ax.set_xlim(xmin - pad, xmax + pad * 2)
    ax.set_ylim(-0.75, n_total - 0.25)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)

    # ── Annotation text in right panel ────────────────────────────────────────
    # y positions in the annotation panel mirror the forest plot
    # ax_txt shares the same y-axis range via twinning of the data space
    # We use axes-fraction coordinates by converting y_pos to fractions.
    ylim_lo, ylim_hi = -0.75, n_total - 0.25
    y_range = ylim_hi - ylim_lo

    col_headers = ["  Δ AUROC [95% CI]", "n (+pCR / −pCR)  p-value"]
    # header at top (just above y_pos[0]=4)
    y_hdr = (y_pos[0] + 0.55 - ylim_lo) / y_range
    ax_txt.text(0.02, y_hdr, col_headers[0], transform=ax_txt.transAxes,
                va="bottom", ha="left", fontsize=8, fontweight="bold", color="#333")
    ax_txt.text(0.02, y_hdr - 0.06, col_headers[1], transform=ax_txt.transAxes,
                va="bottom", ha="left", fontsize=7.5, fontweight="bold", color="#333")

    for i, (d, lo, hi, p, ns) in enumerate(
            zip(deltas, cis_lo, cis_hi, p_vals, n_strs)):
        yf = (y_pos[i] - ylim_lo) / y_range
        line1 = f"{d:+.3f} [{lo:+.3f}, {hi:+.3f}]"
        line2 = f"{ns},  {pretty_p(p)}"
        ax_txt.text(0.02, yf + 0.018, line1, transform=ax_txt.transAxes,
                    va="center", ha="left", fontsize=8.5)
        ax_txt.text(0.02, yf - 0.030, line2, transform=ax_txt.transAxes,
                    va="center", ha="left", fontsize=7.5, color="#444")

    # Summary diamond row
    yf_fe = (diam_y - ylim_lo) / y_range
    fe_line1 = f"{fe_delta:+.3f} [{fe_lo:+.3f}, {fe_hi:+.3f}]"
    fe_line2 = f"Stouffer {pretty_p(stouffer['p_one_sided_2g_gt_fw'])} (1-sided)"
    ax_txt.text(0.02, yf_fe + 0.018, fe_line1, transform=ax_txt.transAxes,
                va="center", ha="left", fontsize=8.5, fontweight="bold")
    ax_txt.text(0.02, yf_fe - 0.030, fe_line2, transform=ax_txt.transAxes,
                va="center", ha="left", fontsize=7.5, fontweight="bold", color="#d62728")

    # ── Figure title (spanning both columns) ──────────────────────────────────
    fig.suptitle(r"Two-gene proliferation$-$ER floor matches the framework across cohorts",
                 fontsize=11, fontweight="bold", y=0.97)

    png = FIG_DIR / "fig3_forest_2gene.png"
    pdf = FIG_DIR / "fig3_forest_2gene.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
