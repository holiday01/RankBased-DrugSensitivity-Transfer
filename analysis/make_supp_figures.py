#!/usr/bin/env python3
"""Regenerate Supplementary Figures S1 and S2 from the corrected deposits.

The versions shipped with the previous submission were computed before the C1 scorer
correction, and their image paths in SUPPLEMENTARY.md resolve to a directory that does not
exist, so neither figure rendered. Both are rebuilt here into revise_bioadv_386/figures/,
which is where '../figures/...' resolves to from manuscript/.

  S1  per-cell residual AUROC forest plot with the pooled fixed-effect estimate
      source: results/A0/A0i_cluster_bootstrap.json  (patient-level cluster bootstrap,
              B = 2000, seed 42; the intervals the addendum requires)
  S2  nested-LRT detection curve against planted residual AUROC, per cohort
      source: results/A1/lrt_power_and_delong_CORRECTED.json -> lrt_power_curve
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386")
FIGDIR = ROOT / "figures"


def forest() -> None:
    d = json.loads((ROOT / "results/A0/A0i_cluster_bootstrap.json").read_text())
    cells = d["cells"]
    labels = [f'{c["name"]} — {c["cohort"]}'
              + (f' ({c["stratum"]})' if c.get("stratum") not in (None, "all") else '')
              + f'\n n = {c["n"]}, {c["events"]} events'
              for c in cells]
    pts = [c["point"] for c in cells]
    los = [c["ci_cluster_percentile"][0] for c in cells]
    his = [c["ci_cluster_percentile"][1] for c in cells]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    y = list(range(len(cells), 0, -1))
    for yi, p, lo, hi in zip(y, pts, los, his):
        ax.plot([lo, hi], [yi, yi], color="#444444", lw=1.4, solid_capstyle="butt")
        ax.plot([p], [yi], marker="o", ms=5.5, color="#1f4e79", zorder=3)

    pooled = d["observed"]["pooled"]
    plo, phi = d["ci_boot_percentile"]
    ax.axhline(0.35, color="#bbbbbb", lw=0.8)
    ax.plot([plo, phi], [-0.35, -0.35], color="#b03030", lw=2.2, solid_capstyle="butt")
    ax.plot([pooled], [-0.35], marker="D", ms=7, color="#b03030", zorder=3)

    ax.axvline(0.5, color="#888888", ls="--", lw=1.0)
    ax.axvline(0.55, color="#cc8800", ls=":", lw=1.0)
    ax.set_yticks(y + [-0.35])
    ax.set_yticklabels(labels + [f"Pooled (k = {len(cells)})\n {pooled:.3f} [{plo:.3f}, {phi:.3f}]"],
                       fontsize=7.2)
    ax.set_xlabel("Residual AUROC after axis decomposition (P + H + B + PAM50-basal)", fontsize=8.5)
    ax.set_xlim(0.30, 0.80)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_title("Figure S1. Per-cell residual AUROC, cluster-bootstrap intervals",
                 fontsize=9.5, loc="left")
    ax.text(0.552, len(cells) + 0.4, "0.55 margin", fontsize=7, color="#cc8800")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"supp_S11_forest.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  S1 forest: {len(cells)} cells, pooled {pooled:.4f} [{plo:.4f}, {phi:.4f}]")


def power() -> None:
    d = json.loads((ROOT / "results/A1/lrt_power_and_delong_CORRECTED.json").read_text())
    curves = d["lrt_power_curve"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    colours = {"GSE25066": "#1f4e79", "GSE22226": "#2e7d32",
               "GSE32646": "#b03030", "GSE16446": "#8a6d3b"}
    for coh, c in sorted(curves.items(), key=lambda kv: -kv[1]["n"]):
        ax.plot(c["grid"], c["detection"], marker="o", ms=4, lw=1.6,
                color=colours.get(coh, "#555555"),
                label=f'{coh} (n = {c["n"]}, MDE$_{{80}}$ = {c["mde80_lrt"]:.2f})')
    ax.axhline(0.80, color="#888888", ls="--", lw=1.0)
    ax.text(0.502, 0.815, "80% power", fontsize=7.5, color="#666666")
    ax.set_xlabel("Planted residual AUROC (axis-orthogonal spike-in)", fontsize=9)
    ax.set_ylabel("Nested-LRT detection probability", fontsize=9)
    ax.set_ylim(0, 1.03)
    ax.set_title("Figure S2. Nested-LRT power against a planted axis-orthogonal signal",
                 fontsize=9.5, loc="left")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.tick_params(labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"supp_S11_power.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  S2 power: {len(curves)} cohorts")


if __name__ == "__main__":
    FIGDIR.mkdir(parents=True, exist_ok=True)
    forest()
    power()
    print(f"wrote into {FIGDIR}")
