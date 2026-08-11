#!/usr/bin/env python3
"""How many genes per direction? An OUTCOME-BLIND criterion.

The manuscript's original justification for N = 30 was that it ranked 6th of 36 in the
ablation grid and sat on a flat AUROC optimum. On the corrected scorer it ranks 30th of 36
(A7b), so that justification is gone. Choosing N by AUROC on the evaluation cohort is not an
option: it is selection on the outcome, which is the exact practice this paper cautions
against, and with n = 488 and ~99 events the whole 36-cell grid (0.457-0.590) sits inside the
noise on a single AUROC estimate.

This script therefore asks a different question, one that never touches patient outcome:
**at what N does gene selection become stable?** Cell lines are resampled with replacement,
the signature is re-derived on each resample using the frozen selection rule, and the
bootstrap set is compared with the full-data set by Jaccard index. Nothing here reads pCR,
any patient cohort, or any clinical endpoint.

One confound is handled explicitly. Jaccard rises with N mechanically, because two larger
sets drawn from the same pool overlap more by chance. For two independent uniform draws of
size N from a pool of M genes, E[|A n B|] = N^2/M and E[|A u B|] = 2N - N^2/M, so
E[J] = N / (2M - N). That baseline is reported next to every observed value, and the excess
over it is what a stability argument may actually lean on.

Outputs results/A7c/gene_count_stability.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/holiday01/2026_ISMB_code/revise_bioadv/scripts")
import freeze_signatures as F  # noqa: E402

DRUG = "DOXORUBICIN"
STRATEGY = "quantile_30_oncokb_all"   # the frozen signature's rule
N_GRID = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100]
N_BOOT = 100                           # matches the deposited convention
SEED = 42
OUT = Path("/home/holiday01/2026_ISMB_code/revise_bioadv_386/results/A7c")


def expected_jaccard(n: int, m: int) -> float:
    """E[Jaccard] for two independent uniform size-n draws from a pool of m."""
    return n / (2 * m - n)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("loading CCLE ...", flush=True)
    ccle_drug = pd.read_csv(F.CCLE_DRUG_PATH)
    ccle_exp = pd.read_csv(F.CCLE_EXPR_PATH)
    oncokb = pd.read_csv(F.ONCOKB_PATH, sep="\t")
    model_df = pd.read_csv(F.CCLE_MODEL_PATH, low_memory=False)

    onco_all, _ = F.build_oncokb_sets(oncokb)
    expr_pan = F.build_expr_base(ccle_drug, ccle_exp)
    heme = F.identify_heme_model_ids(model_df)
    expr_solid = expr_pan.loc[~expr_pan.index.isin(heme)]          # frozen signature is solid-only
    pool = expr_solid[[g for g in expr_solid.columns if g in onco_all]]
    drug_map, _ = F.build_ccle_drug_map_with_collisions(ccle_drug)
    strat = F.STRATEGY_BY_NAME[STRATEGY]
    m = int(pool.shape[1])
    print(f"  solid-only cell lines: {pool.shape[0]}, OncoKB-all gene pool: {m}", flush=True)

    rows = []
    original_n = F.QUANTILE_TOP_N
    try:
        for n in N_GRID:
            F.QUANTILE_TOP_N = n
            sig = F.build_signature(DRUG, strat, ccle_drug, drug_map, pool)
            if sig is None or "_skip_reason" in sig:
                print(f"  N={n}: skipped ({sig})", flush=True)
                continue
            orig_s = {g for g, _ in sig["up_sensitive"]}
            orig_r = {g for g, _ in sig["up_resistant"]}
            st = F.stability_jaccard(DRUG, strat, ccle_drug, drug_map, pool,
                                     orig_s, orig_r, n_bootstrap=N_BOOT, seed=SEED)
            ej = expected_jaccard(n, m)
            row = {
                "N": n,
                "gene_pool_size": m,
                "jaccard_sensitivity_median": st["jaccard_sens_median"],
                "jaccard_sensitivity_iqr": st["jaccard_sens_iqr"],
                "jaccard_resistance_median": st["jaccard_res_median"],
                "jaccard_resistance_iqr": st["jaccard_res_iqr"],
                "expected_jaccard_random": ej,
                "excess_sensitivity": st["jaccard_sens_median"] - ej,
                "excess_resistance": st["jaccard_res_median"] - ej,
                "n_bootstrap_successful": st["n_bootstrap_successful"],
            }
            rows.append(row)
            print(f"  N={n:>3}: sens {row['jaccard_sensitivity_median']:.3f} "
                  f"res {row['jaccard_resistance_median']:.3f} "
                  f"(random baseline {ej:.4f})", flush=True)
    finally:
        F.QUANTILE_TOP_N = original_n

    out = {
        "question": "At what number of genes per direction does selection become stable? "
                    "Outcome-blind: no patient cohort, pCR call or clinical endpoint is read.",
        "drug": DRUG,
        "strategy": STRATEGY,
        "cell_line_pool": "CCLE solid-only",
        "n_bootstrap": N_BOOT,
        "seed": SEED,
        "expected_jaccard_formula": "N / (2M - N) for two independent uniform size-N draws "
                                    "from a pool of M genes",
        "rows": rows,
    }
    (OUT / "gene_count_stability.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {OUT / 'gene_count_stability.json'}")


if __name__ == "__main__":
    main()
