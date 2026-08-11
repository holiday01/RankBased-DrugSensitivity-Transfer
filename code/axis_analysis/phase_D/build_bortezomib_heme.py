"""Track D — lineage-matched BORTEZOMIB signature trained ONLY on hematologic
cell lines (the lineage of the myeloma validation cohort). Removes the
pan-cancer 'heme-vs-solid lineage detector' confound that dominates the
pan-cancer signature (top sensitivity genes were all heme markers: LCP1,
IKZF1, PTPRC...). Same quantile_30_oncokb_all gene-selection algorithm,
but expr pool restricted to heme lines; per-group minimum relaxed to 30
(documented deviation; ~40/group available from 133 heme lines at q=0.30).

NOTE: CTRP bortezomib AUC has 133 heme lines but 0 labeled myeloma; heme =
leukemia/lymphoid/plasma. This is the closest available lineage match.
"""
from __future__ import annotations
import sys, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
import pandas as pd

REPO = Path("/home/holiday01/2026_ISMB_code")
sys.path.insert(0, str(REPO / "revise_bioadv" / "scripts"))
import freeze_signatures as fs

OUT_DIR = REPO / "revise_bioadv" / "resubmission_v2" / "results" / "phase_D"
OUT_DIR.mkdir(parents=True, exist_ok=True)
STRATEGY = dict(fs.STRATEGY_BY_NAME["quantile_30_oncokb_all"])
DRUG = "BORTEZOMIB"

# documented deviation: relax per-group min for lineage-restricted training
fs.MIN_SAMPLES_PER_GROUP = 30
fs.MIN_CELL_LINES = 80


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    log = {"started_utc": datetime.now(timezone.utc).isoformat(), "variant": "heme_only",
           "strategy": STRATEGY, "min_samples_per_group": fs.MIN_SAMPLES_PER_GROUP}

    ctd = pd.read_csv(REPO / "CCLE25Q3" / "Drug_sensitivity_AUC_(CTD^2)_subsetted.csv")
    id_col = ctd.columns[0]
    if not ctd[id_col].astype(str).str.startswith("ACH-").mean() > 0.5:
        for c in ctd.columns:
            if ctd[c].astype(str).str.startswith("ACH-").mean() > 0.5:
                id_col = c; break
    bort_col = [c for c in ctd.columns if c.upper().startswith("BORTEZOMIB")][0]
    auc = ctd.set_index(id_col)[bort_col].astype(float).dropna()
    auc = auc[~auc.index.duplicated(keep="first")]

    # heme model ids
    model_df = pd.read_csv(fs.CCLE_MODEL_PATH)
    heme = fs.identify_heme_model_ids(model_df)
    auc_heme = auc[[m for m in auc.index if m in heme]]
    log["n_heme_lines_with_bort_auc"] = int(len(auc_heme))

    ccle_exp = pd.read_csv(fs.CCLE_EXPR_PATH)
    oncokb = pd.read_csv(fs.ONCOKB_PATH, sep="\t")
    oncogenes_all, _ = fs.build_oncokb_sets(oncokb)

    drug_frame = auc_heme.reset_index()
    drug_frame.columns = ["ModelID", "BORTEZOMIB"]
    expr_base = fs.build_expr_base(drug_frame, ccle_exp)
    expr_oncokb = expr_base[[g for g in expr_base.columns if g in oncogenes_all]]
    log["heme_pool_cells"] = int(expr_base.shape[0])
    log["oncokb_all_gene_pool_size"] = int(expr_oncokb.shape[1])

    sig = fs.build_signature(DRUG, STRATEGY, None, None, expr_oncokb, auc_override=auc_heme)
    if sig is None or "_skip_reason" in sig:
        log["result"] = "SKIP"; log["skip"] = sig
        (OUT_DIR / "build_bortezomib_heme_log.json").write_text(json.dumps(log, indent=1))
        raise SystemExit(f"heme signature build skipped: {sig}")

    rows = []
    for rank, (g, d) in enumerate(sig["up_sensitive"], 1):
        rows.append({"gene_symbol": g, "direction": "sensitivity", "rank_in_direction": rank, "delta": d})
    for rank, (g, d) in enumerate(sig["up_resistant"], 1):
        rows.append({"gene_symbol": g, "direction": "resistance", "rank_in_direction": rank, "delta": d})
    out_tsv = OUT_DIR / "BORTEZOMIB_heme.tsv"
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False, float_format="%.9g")

    log.update(result="OK", n_sensitive_cells=sig["n_sensitive_cells"],
               n_resistant_cells=sig["n_resistant_cells"], n_cell_lines=sig["n_cell_lines"],
               signature_sha256=sha(out_tsv),
               top10_sensitivity=[g for g, _ in sig["up_sensitive"][:10]],
               top10_resistance=[g for g, _ in sig["up_resistant"][:10]],
               finished_utc=datetime.now(timezone.utc).isoformat())
    (OUT_DIR / "build_bortezomib_heme_log.json").write_text(json.dumps(log, indent=1))
    print("DONE heme variant ->", out_tsv, "| n_cell_lines:", sig["n_cell_lines"])
    print("top sensitivity:", log["top10_sensitivity"])
    print("top resistance :", log["top10_resistance"])


if __name__ == "__main__":
    main()
