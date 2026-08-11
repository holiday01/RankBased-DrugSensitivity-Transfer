"""Track D — build bidirectional CCLE-derived BORTEZOMIB (proteasome-inhibitor)
signature using the EXACT frozen-signature methodology (freeze_signatures.py),
via auc_override, for the new-drug-class validation (myeloma / GSE9782).

Primary config (matches original primary): pan_cancer / quantile_30_oncokb_all
 - pan_cancer pool INCLUDES heme/myeloma cell lines (correct lineage for the
   bortezomib->myeloma transfer, unlike the solid-tumour Dox case).
 - q=0.30 quantile split, OncoKB-all gene universe, top-30 per direction,
   bidirectional (sensitivity + resistance).

AUC source: CCLE25Q3/Drug_sensitivity_AUC_(CTD^2)_subsetted.csv (BORTEZOMIB col),
keyed by DepMap ModelID. Same CTRP/CTD^2 provenance as the existing frozen sigs.
"""
from __future__ import annotations
import sys, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

REPO = Path("/home/holiday01/2026_ISMB_code")
sys.path.insert(0, str(REPO / "revise_bioadv" / "scripts"))
import freeze_signatures as fs  # reuse exact build functions

OUT_DIR = REPO / "revise_bioadv" / "resubmission_v2" / "results" / "phase_D"
OUT_DIR.mkdir(parents=True, exist_ok=True)

STRATEGY = fs.STRATEGY_BY_NAME["quantile_30_oncokb_all"]
DRUG_NORM = "BORTEZOMIB"


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main():
    log = {"started_utc": datetime.now(timezone.utc).isoformat(),
           "strategy": STRATEGY, "drug": DRUG_NORM,
           "freeze_script_sha256": sha(REPO / "revise_bioadv/scripts/freeze_signatures.py")}

    # --- 1. bortezomib AUC keyed by ModelID (CTD^2) ---
    ctd_path = REPO / "CCLE25Q3" / "Drug_sensitivity_AUC_(CTD^2)_subsetted.csv"
    ctd = pd.read_csv(ctd_path)
    # detect id column
    id_col = None
    for c in ctd.columns:
        cl = c.lower()
        if cl in ("modelid", "depmap_id", "depmapid") or cl.startswith("ach-"):
            id_col = c; break
    if id_col is None:
        # first column usually the id / 'Unnamed: 0'
        first = ctd.columns[0]
        if ctd[first].astype(str).str.startswith("ACH-").mean() > 0.5:
            id_col = first
    if id_col is None:
        raise SystemExit(f"could not find ModelID column; cols={list(ctd.columns)[:6]}")
    bort_col = [c for c in ctd.columns if c.upper().startswith("BORTEZOMIB")]
    if not bort_col:
        raise SystemExit("BORTEZOMIB column not found in CTD^2")
    bort_col = bort_col[0]
    auc = ctd.set_index(id_col)[bort_col].astype(float).dropna()
    auc = auc[~auc.index.duplicated(keep="first")]
    log["ctd_path"] = str(ctd_path); log["ctd_sha256"] = sha(ctd_path)
    log["bort_col"] = bort_col; log["id_col"] = id_col
    log["n_cells_with_bort_auc"] = int(len(auc))

    # --- 2. CCLE expression + OncoKB, build pan-cancer oncokb_all pool ---
    ccle_exp = pd.read_csv(fs.CCLE_EXPR_PATH)
    ccle_exp = ccle_exp.rename(columns={ccle_exp.columns[0]: "ModelID"}) \
        if "ModelID" not in ccle_exp.columns else ccle_exp
    oncokb = pd.read_csv(fs.ONCOKB_PATH, sep="\t")
    oncogenes_all, oncogenes_only = fs.build_oncokb_sets(oncokb)

    # build_expr_base wants a drug frame whose col0 = model ids (to subset expr)
    drug_frame = auc.reset_index()
    drug_frame.columns = ["ModelID", "BORTEZOMIB"]
    expr_base = fs.build_expr_base(drug_frame, ccle_exp)  # pan-cancer (all lineages w/ bort AUC)
    expr_oncokb = expr_base[[g for g in expr_base.columns if g in oncogenes_all]]
    log["pan_cancer_pool_cells"] = int(expr_base.shape[0])
    log["oncokb_all_gene_pool_size"] = int(expr_oncokb.shape[1])

    # heme composition diagnostic (myeloma lineage alignment)
    try:
        model_df = pd.read_csv(fs.CCLE_MODEL_PATH)
        heme = fs.identify_heme_model_ids(model_df)
        log["n_heme_lines_in_pool"] = int(sum(1 for m in expr_base.index if m in heme))
        # myeloma-specific count
        dis = model_df.set_index("ModelID")["OncotreePrimaryDisease"].astype(str).str.lower()
        mm_ids = set(dis[dis.str.contains("myeloma", na=False)].index)
        log["n_myeloma_lines_in_pool"] = int(sum(1 for m in expr_base.index if m in mm_ids))
    except Exception as e:
        log["heme_diag_error"] = str(e)

    # --- 3. build signature via auc_override (identical algorithm) ---
    sig = fs.build_signature(DRUG_NORM, STRATEGY, None, None, expr_oncokb, auc_override=auc)
    if sig is None or "_skip_reason" in sig:
        log["result"] = "SKIP"; log["skip"] = sig
        (OUT_DIR / "build_bortezomib_log.json").write_text(json.dumps(log, indent=1))
        raise SystemExit(f"signature build skipped: {sig}")

    rows = []
    for rank, (g, d) in enumerate(sig["up_sensitive"], 1):
        rows.append({"gene_symbol": g, "direction": "sensitivity", "rank_in_direction": rank, "delta": d})
    for rank, (g, d) in enumerate(sig["up_resistant"], 1):
        rows.append({"gene_symbol": g, "direction": "resistance", "rank_in_direction": rank, "delta": d})
    sig_df = pd.DataFrame(rows)
    out_tsv = OUT_DIR / "BORTEZOMIB.tsv"
    sig_df.to_csv(out_tsv, sep="\t", index=False, float_format="%.9g")

    log["result"] = "OK"
    log["n_sensitive_cells"] = sig["n_sensitive_cells"]
    log["n_resistant_cells"] = sig["n_resistant_cells"]
    log["n_cell_lines"] = sig["n_cell_lines"]
    log["dropoff_ratio_sensitivity"] = sig["dropoff_ratio_sensitivity"]
    log["dropoff_ratio_resistance"] = sig["dropoff_ratio_resistance"]
    log["signature_sha256"] = sha(out_tsv)
    log["top10_sensitivity"] = [g for g, _ in sig["up_sensitive"][:10]]
    log["top10_resistance"] = [g for g, _ in sig["up_resistant"][:10]]
    log["finished_utc"] = datetime.now(timezone.utc).isoformat()
    (OUT_DIR / "build_bortezomib_log.json").write_text(json.dumps(log, indent=1))
    print("DONE: bortezomib signature ->", out_tsv)
    print("n_cell_lines:", sig["n_cell_lines"],
          "| heme:", log.get("n_heme_lines_in_pool"),
          "| myeloma:", log.get("n_myeloma_lines_in_pool"))
    print("top sensitivity:", log["top10_sensitivity"])
    print("top resistance :", log["top10_resistance"])


if __name__ == "__main__":
    main()
