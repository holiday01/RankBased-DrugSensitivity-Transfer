#!/usr/bin/env python3
"""Fig 3: Forest plot of per-cohort Delta AUROC (2-gene vs framework).

Reads results/phase_4/critical_gaps/meta_stats.json + summary tsv,
draws a forest plot with per-cohort points + 95% CI and a fixed-effects
meta-Delta summary diamond. Annotates Stouffer combined p.

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
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "lines.linewidth": 1.5,
})

ROOT = Path(".")
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

META_JSON = ROOT / "results/phase_4/critical_gaps/meta_stats.json"


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

    # cohort order (top to bottom on y-axis)
    desired = ["GSE16446", "GSE25066", "GSE22226", "GSE32646"]
    rows = []
    by_name = {r["cohort"]: r for r in cohort_rows}
    for n in desired:
        if n in by_name:
            rows.append(by_name[n])
    if not rows:
        print("[ERR] no per-cohort rows")
        return 1

    # ---------- compute per-cohort 95% CI from se_delta ----------
    z975 = 1.959963984540054
    cohorts = []
    deltas = []
    cis_lo = []
    cis_hi = []
    p_vals = []
    n_strs = []
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

    # ---------- forest plot ----------
    fig, ax = plt.subplots(figsize=(3.5, 3.5), constrained_layout=True)

    labels = cohorts + ["Fixed-effects meta"]
    n_total = len(labels)
    y_pos = np.arange(n_total)[::-1]  # top->bottom

    # cohort markers
    for i, (d, lo, hi) in enumerate(zip(deltas, cis_lo, cis_hi)):
        ypos = y_pos[i]
        ax.errorbar(d, ypos, xerr=[[d - lo], [hi - d]],
                    fmt="s", color="#1f77b4", ecolor="#1f77b4",
                    elinewidth=1.5, capsize=3, markersize=6,
                    markeredgecolor="white", markeredgewidth=0.5)

    # summary diamond at bottom
    diamond_y = y_pos[-1]
    diamond_x = [fe_lo, fe_delta, fe_hi, fe_delta]
    diamond_yy = [diamond_y, diamond_y + 0.3, diamond_y, diamond_y - 0.3]
    ax.fill(diamond_x, diamond_yy, color="#d62728", alpha=0.85, ec="black",
            lw=1.0, zorder=5)

    # vertical line at delta=0
    ax.axvline(0, color="0.3", ls="--", lw=1.0, zorder=1)

    # y-axis labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel(r"$\Delta$ AUROC  (2-gene $-$ framework)")
    ax.set_title("Two-gene IHC4-proxy baseline matches the framework across cohorts",
                 fontsize=11)

    # annotate per-cohort delta + p
    xmax = max(cis_hi + [fe_hi])
    xmin = min(cis_lo + [fe_lo, 0])
    pad = (xmax - xmin) * 0.05
    ax.set_xlim(xmin - pad, xmax + pad + (xmax - xmin) * 0.35)

    for i, (d, lo, hi, p, ns) in enumerate(
            zip(deltas, cis_lo, cis_hi, p_vals, n_strs)):
        ypos = y_pos[i]
        txt = f" {d:+.3f} [{lo:+.3f}, {hi:+.3f}]\n {ns}, {pretty_p(p)}"
        ax.text(xmax + pad, ypos, txt, va="center", ha="left", fontsize=7)

    txt = (f" {fe_delta:+.3f} [{fe_lo:+.3f}, {fe_hi:+.3f}]\n"
           f" Stouffer {pretty_p(stouffer['p_one_sided_2g_gt_fw'])} (1-sided)")
    ax.text(xmax + pad, diamond_y, txt, va="center", ha="left",
            fontsize=7, fontweight="bold")

    ax.grid(True, axis="x", linestyle=":", alpha=0.4)
    ax.set_axisbelow(True)

    png = FIG_DIR / "fig3_forest_2gene.png"
    pdf = FIG_DIR / "fig3_forest_2gene.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote {png}\n[OK] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
