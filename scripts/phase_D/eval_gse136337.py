"""Track D replication — GSE136337 (UAMS Total Therapy long-term cohort, gene-
level GPL27143, n~426). Bortezomib-containing frontline TT protocols + survival.
CAVEAT: UAMS TT -> likely PARTIAL patient overlap with GSE2658 (not fully
independent); frontline/TT/survival setting (same that failed). Reported anyway.

OS from mixed-format dates (ISO YYYY-MM-DD or Excel serial). Cox HR per SD.
"""
from __future__ import annotations
import sys, json, re
from datetime import datetime, timezone, date
from pathlib import Path
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

REPO = Path(".")
sys.path.insert(0, str(REPO / "." / "scripts"))
from parse_gse_to_matrix import read_series_matrix
from score_external_cohort import compute_singscore, load_signature

COHORT = REPO / "." / "external_data" / "GSE136337"
PHASE_D = REPO / "." / "." / "results" / "phase_D"
OUT = PHASE_D / "gse136337_eval"; OUT.mkdir(parents=True, exist_ok=True)
EXCEL_EPOCH = date(1899, 12, 30)


def parse_date(s):
    s = str(s).strip().strip('"')
    if s in (".", "", "na", "NA", "nd"):
        return None
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})$", s)
    if m:
        return date(int(m[1]), int(m[2]), int(m[3]))
    if re.match(r"^\d+(\.0)?$", s):
        n = int(float(s))
        if 10000 < n < 60000:
            return date.fromordinal(EXCEL_EPOCH.toordinal() + n)
    return None


def char_dict(meta):
    ids = meta["sample_ids"]; lines = meta["sample_chars"]["characteristics"]
    cols = {}
    for cells in lines:  # each is already a per-sample value list
        if not cells:
            continue
        key = cells[0].split(":")[0].strip() if ":" in cells[0] else None
        if key and len(cells) == len(ids):
            cols[key] = [c.split(":", 1)[1].strip() if ":" in c else "" for c in cells]
    df = pd.DataFrame(cols); df.index = ids
    return df


def cox(score, time, event):
    df = pd.DataFrame({"t": time, "e": event, "z": (score - np.nanmean(score)) / np.nanstd(score)}).dropna()
    df = df[df["t"] > 0]
    if df["e"].sum() < 5 or len(df) < 20:
        return None
    c = CoxPHFitter(); c.fit(df, "t", "e")
    return {"HR_per_SD": float(np.exp(c.params_["z"])),
            "HR_ci": [float(np.exp(c.confidence_intervals_.iloc[0, 0])), float(np.exp(c.confidence_intervals_.iloc[0, 1]))],
            "p": float(c.summary.loc["z", "p"]), "c_index": float(c.concordance_index_),
            "n": int(len(df)), "events": int(df["e"].sum())}


def main():
    res = {"started_utc": datetime.now(timezone.utc).isoformat(), "cohort": "GSE136337",
           "caveat": "UAMS Total Therapy; likely partial overlap with GSE2658; frontline/TT/survival"}
    expr_g, meta = read_series_matrix(COHORT / "series_matrix.txt.gz")  # gene x sample
    if "sample_ids" not in meta:
        meta["sample_ids"] = list(expr_g.columns)
    expr = expr_g.T  # sample x gene
    expr.columns = [str(c).upper() for c in expr.columns]
    res["n_samples"], res["n_genes"] = int(expr.shape[0]), int(expr.shape[1])

    cd = char_dict(meta)
    # OS: baseline=dateenrolled (fallback datesamplernbx); death=datedeath; censor=datelastcontact
    base_col = "dateenrolled" if "dateenrolled" in cd else ("datesamplernbx" if "datesamplernbx" in cd else None)
    os_time, os_evt = {}, {}
    for sid in cd.index:
        base = parse_date(cd.loc[sid, base_col]) if base_col else None
        dd = parse_date(cd.loc[sid, "datedeath"]) if "datedeath" in cd else None
        lc = parse_date(cd.loc[sid, "datelastcontact"]) if "datelastcontact" in cd else None
        end = dd or lc
        if base and end:
            os_time[sid] = (end.toordinal() - base.toordinal()) / 30.44
            os_evt[sid] = 1 if dd else 0
    res["n_with_OS"] = len(os_time)
    if "ftprotocol" in cd:
        res["protocol_counts"] = cd["ftprotocol"].value_counts().head(12).to_dict()

    common = [s for s in expr.index if s in os_time]
    expr = expr.loc[common]
    t = np.array([os_time[s] for s in common]); e = np.array([os_evt[s] for s in common])
    rank_pct = expr.rank(axis=1, method="average", pct=True)

    res["signatures"], res["baselines"] = {}, {}
    for label, path in [("pan_cancer", PHASE_D / "BORTEZOMIB.tsv"),
                        ("heme_matched", PHASE_D / "BORTEZOMIB_heme.tsv")]:
        sens, resg, _, info = load_signature(path)
        try:
            sc, sd = compute_singscore(expr, sens, resg)
            res["signatures"][label] = cox(sc.loc[common].to_numpy(), t, e)
            res["signatures"][label]["genes_present"] = [sd["n_up_genes_present_in_cohort"], sd["n_down_genes_present_in_cohort"]]
        except Exception as ex:
            res["signatures"][label] = {"error": str(ex)}
    for g in ["MKI67", "TNFRSF17"]:
        if g in rank_pct.columns:
            res["baselines"][g] = cox(rank_pct[g].to_numpy(), t, e)

    res["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT / "gse136337_eval.json").write_text(json.dumps(res, indent=2))
    print(f"n={res['n_samples']} genes={res['n_genes']} with_OS={res['n_with_OS']}")
    print("protocols:", res.get("protocol_counts"))
    for label, c in res["signatures"].items():
        if c and "HR_per_SD" in c:
            print(f"  {label:12s} HR/SD {c['HR_per_SD']:.3f} [{c['HR_ci'][0]:.3f},{c['HR_ci'][1]:.3f}] p={c['p']:.3g} C={c['c_index']:.3f} (n={c['n']},ev={c['events']})")
    for g, c in res["baselines"].items():
        if c and "HR_per_SD" in c:
            print(f"  baseline {g:9s} HR/SD {c['HR_per_SD']:.3f} p={c['p']:.3g} C={c['c_index']:.3f}")
    print("NOTE: HR/SD>1 = higher score -> worse OS (proliferation-like); sensitivity predictor would be HR<1.")
    print("Saved ->", OUT / "gse136337_eval.json")


if __name__ == "__main__":
    main()
