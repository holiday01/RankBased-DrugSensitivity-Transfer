#!/usr/bin/env python3
"""A14 stage a -- build GSE41998 and evaluate the OUTCOME-BLIND gate.

Addendum section 7, as amended by section 13.5. GSE41998 is the optional second independent
neoadjuvant pCR cohort. Its gate is outcome-blind by construction and one-way:

    >= 200 pCR calls
    >= 20 events per class
    >= 0.60 up-set AND down-set gene resolution
    a documented patient-overlap check against all four existing cohorts (13.5)
    verdict and UTC timestamp written to decisions_addendum.md BEFORE any residual AUROC or
    LRT is computed
    hard date 2026-08-26; if it passes, abandoned without trace in the manuscript

**This script deliberately contains no code that computes a residual AUROC, a nested LRT, or
any other outcome-dependent quantity.** The ordering is enforced structurally rather than by
discipline: stage b (A14b) holds that code and refuses to run until the gate verdict exists.
Event counts are computed here because they are gate criteria, not analysis.

Provenance note: NCBI/GEO returns 403 to this host, so the series was retrieved from the EBI
BioStudies mirror (E-GEOD-41998) instead. Same submission, same 279 samples.

Platform is A-AFFY-37 = Affymetrix HG-U133A (GPL96), the same platform as GSE25066, so the
already-deposited GPL96 annotation is reused and no further external request is made.

**The defective HGNC alias table of section 13.1 is NOT applied at matrix-build time.** That is
what collapsed HGF, SMO, CDH1 and MET into unrelated genes in the existing cohorts. GSE41998 is
built with approved symbols left intact, so it does not inherit the corruption; the difference
is reported.

Writes external_data/GSE41998/{expression,axes,phenotype}.tsv and results/A14/A14a_gate.json.
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import A0e_dedup as D  # noqa: E402

G = D.ROOT / "external_data/GSE41998"
ANNOT = D.ROOT / "external_data/GSE25066/platform_annot.annot.gz"   # GPL96, same platform
GMT = D.ROOT / "h.all.v2025.1.Hs.symbols.gmt"
CENTROIDS = D.ROOT / "external_data/pam50_centroids.tsv"
OUT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A14")
HALLMARK = ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
            "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1"]
MIN_CALLS, MIN_EVENTS, MIN_RES = 200, 20, 0.60


def build_expression() -> pd.DataFrame:
    files = sorted((G / "samples").glob("*_sample_table.txt"))
    print(f"assembling {len(files)} per-sample tables")
    cols = {}
    for f in files:
        gsm = f.name.split("_")[0]
        t = pd.read_csv(f, sep="\t", index_col=0)
        cols[gsm] = pd.to_numeric(t.iloc[:, 0], errors="coerce")
    probes = pd.DataFrame(cols)
    print(f"  probe matrix {probes.shape}")

    ann = pd.read_csv(ANNOT, sep="\t", skiprows=27, dtype=str, low_memory=False)
    ann = ann[["ID", "Gene symbol"]].dropna()
    # a probe mapping to several genes ("A///B") is ambiguous and is dropped, not split
    ann = ann[~ann["Gene symbol"].str.contains("///", na=False)]
    m = dict(zip(ann["ID"], ann["Gene symbol"]))
    probes = probes.loc[probes.index.intersection(list(m))]
    sym = pd.Series([m[i] for i in probes.index], index=probes.index)
    # collapse duplicate probes to the highest-IQR probe per symbol (the deposited convention)
    iqr = probes.quantile(0.75, axis=1) - probes.quantile(0.25, axis=1)
    keep = iqr.groupby(sym).idxmax()
    expr = probes.loc[keep.values]
    expr.index = keep.index
    expr = expr.T
    print(f"  gene matrix {expr.shape} (samples x genes)")
    return expr


def load_gmt() -> dict:
    out = {}
    for line in GMT.read_text().splitlines():
        p = line.split("\t")
        if p and p[0] in HALLMARK:
            out[p[0]] = [g for g in p[2:] if g]
    return out


def unsigned_singscore(rank_pct: pd.DataFrame, genes: list[str]) -> pd.Series:
    g = [x for x in genes if x in rank_pct.columns]
    n = len(g)
    mn = (1 + n) / (2 * rank_pct.shape[1])
    return (rank_pct[g].mean(axis=1) - mn) / ((1 - mn) - mn)


def z(s):
    return (s - s.mean()) / s.std(ddof=0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expr = build_expression()
    expr.to_csv(G / "GSE41998_expression.tsv", sep="\t")

    # ---- phenotype from the SDRF ----
    s = pd.read_csv(G / "E-GEOD-41998.sdrf.txt", sep="\t", dtype=str)
    s["patient_id"] = s["Source Name"].str.split().str[0]
    ph = pd.DataFrame({
        "patient_id": s["patient_id"],
        "pcr_raw": s["Characteristics [pcr]"].astype(str).str.strip(),
        "pcr_rcb01_raw": s["Characteristics [pcrrcb1]"].astype(str).str.strip(),
        "er_status_ihc": s["Characteristics [er]"].astype(str).str.strip(),
        "pr_status_ihc": s["Characteristics [pr]"].astype(str).str.strip(),
        "her2_status": s["Characteristics [her2stat]"].astype(str).str.strip(),
        "age_years": pd.to_numeric(s["Characteristics [age]"], errors="coerce"),
        "treatment_arm": s["Characteristics [treatment arm]"].astype(str).str.strip(),
    }).set_index("patient_id")
    ph["pCR"] = ph["pcr_raw"].map({"Yes": 1, "No": 0})
    ph.to_csv(G / "GSE41998_phenotype.tsv", sep="\t")

    common = expr.index.intersection(ph.index)
    ph, expr = ph.loc[common], expr.loc[common]
    lab = ph["pCR"].dropna()
    n_calls, n_pos, n_neg = len(lab), int((lab == 1).sum()), int((lab == 0).sum())
    print(f"pCR calls {n_calls}  events {n_pos}  non-events {n_neg}")

    # ---- axes ----
    rank_pct = expr.rank(axis=1, pct=True)
    sets = load_gmt()
    P = pd.concat([unsigned_singscore(rank_pct, g) for g in sets.values()], axis=1).mean(axis=1)
    def col(g):
        return expr[g] if g in expr.columns else pd.Series(np.nan, index=expr.index)
    H = pd.concat([z(col("ESR1")), z(col("PGR"))], axis=1).mean(axis=1)
    axes = pd.DataFrame({"P": z(P), "H": H, "H_ESR1only": z(col("ESR1")),
                         "B": z(col("ERBB2"))}, index=expr.index)

    pam = OUT / "GSE41998_pam50.tsv"
    subprocess.run([sys.executable, str(D.ROOT / "scripts/compute_pam50.py"),
                    "--cohort", "GSE41998",
                    "--expression", str(G / "GSE41998_expression.tsv"),
                    "--centroids", str(CENTROIDS), "--out-tsv", str(pam),
                    "--out-log", str(OUT / "GSE41998_pam50_log.json")], check=True)
    p50 = pd.read_csv(pam, sep="\t", index_col=0)
    bcol = [c for c in p50.columns if "basal" in c.lower()][0]
    axes["H_PAM50basal"] = z(p50[bcol].reindex(axes.index))
    axes["pCR"] = ph["pCR"]
    axes.dropna(subset=["P", "H", "B", "H_PAM50basal", "pCR"]).to_csv(
        G / "GSE41998_axes.tsv", sep="\t")
    n_axes = int(axes.dropna(subset=["P", "H", "B", "H_PAM50basal", "pCR"]).shape[0])
    print(f"complete four-axis + pCR rows: {n_axes}")

    # ---- gene resolution for the frozen framework signature ----
    fw = pd.read_csv(D.ROOT / "frozen_signatures/solid_only/quantile_30_oncokb_all/"
                              "DOXORUBICIN.tsv", sep="\t")
    up = list(fw.loc[fw.direction == "sensitivity", "gene_symbol"])
    dn = list(fw.loc[fw.direction == "resistance", "gene_symbol"])
    cols = {c.upper() for c in expr.columns}
    fu = sum(g.upper() in cols for g in up) / len(up)
    fd = sum(g.upper() in cols for g in dn) / len(dn)
    print(f"gene resolution  up {fu:.3f}  down {fd:.3f}")
    collided = [g for g in ["HGF", "SMO", "CDH1", "MET"] if g.upper() in cols]
    print(f"symbols preserved that the alias table destroys elsewhere: {collided}")

    overlap = json.loads((G / "OVERLAP_CHECK.json").read_text())
    tot_overlap = sum(v["gsm_overlap"] + v["identifier_overlap"]
                      for k, v in overlap.items() if isinstance(v, dict) and "gsm_overlap" in v)

    crit = {
        "pcr_calls": dict(value=n_calls, threshold=MIN_CALLS, pass_=n_calls >= MIN_CALLS),
        "events_positive": dict(value=n_pos, threshold=MIN_EVENTS, pass_=n_pos >= MIN_EVENTS),
        "events_negative": dict(value=n_neg, threshold=MIN_EVENTS, pass_=n_neg >= MIN_EVENTS),
        "up_resolution": dict(value=round(fu, 4), threshold=MIN_RES, pass_=fu >= MIN_RES),
        "down_resolution": dict(value=round(fd, 4), threshold=MIN_RES, pass_=fd >= MIN_RES),
        "patient_overlap": dict(value=tot_overlap, threshold=0, pass_=tot_overlap == 0),
    }
    verdict = "PASS" if all(c["pass_"] for c in crit.values()) else "FAIL"
    res = dict(cohort="GSE41998", source="EBI BioStudies E-GEOD-41998 (NCBI returns 403 here)",
               platform="A-AFFY-37 / GPL96 (HG-U133A), annotation reused from GSE25066",
               n_samples=int(expr.shape[0]), n_genes=int(expr.shape[1]),
               n_complete_axes_pcr=n_axes,
               alias_table_applied=False,
               alias_note=("Section 13.1: the HGNC alias table is defective and its application "
                           "at matrix-build time is what collapsed HGF/SMO/CDH1/MET into "
                           "unrelated genes in the existing cohorts. It is not applied here."),
               symbols_preserved=collided,
               criteria=crit, verdict=verdict,
               utc=datetime.now(timezone.utc).isoformat(),
               gate_note=("Outcome-blind. No residual AUROC, nested LRT or other "
                          "outcome-dependent quantity is computed in this stage; the code for "
                          "those lives in A14b and is gated on this verdict."))
    (OUT / "A14a_gate.json").write_text(json.dumps(res, indent=2))
    print("\n--- GATE ---")
    for k, c in crit.items():
        print(f"  {k:18s} {str(c['value']):>8s}  vs {c['threshold']}  "
              f"{'PASS' if c['pass_'] else 'FAIL'}")
    print(f"  VERDICT: {verdict}   ({res['utc']})")
    print(f"\nwrote {OUT / 'A14a_gate.json'}")
    print("Next: write the verdict to decisions_addendum.md, THEN run A14b.")


if __name__ == "__main__":
    main()
