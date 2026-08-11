"""Regenerate results/phase_9_revalidation/meta_pools.json from meta_cells.tsv.

Fixed-effect inverse-variance pools of the per-cell residual AUROCs:
  * "ALL k=11"      : every cell (reported in Supplementary S11 for completeness only).
  * "pCR-only k=9"  : drops BOTH survival-endpoint METABRIC cells (the two Mammaprint
                      cells: METABRIC cohort or name starting with 'Mammaprint'); this is
                      the headline pool used in the abstract / Results / Figure 5 / S11.

The earlier revalidation.py "exclude METABRIC" pass filtered on cohort=='METABRIC' only,
which dropped one of the two survival cells (k=10); this script drops both (k=9), matching
the headline pooled residual 0.501 [0.474, 0.529]. Pooled point/CI are stored rounded to
4 dp (as deposited); TOST p-values are full precision.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parents[1] / "results" / "phase_9_revalidation"
cells = pd.read_csv(R / "meta_cells.tsv", sep="\t")
is_surv = (cells["cohort"] == "METABRIC") | (cells["name"].str.startswith("Mammaprint"))


def pool(sub):
    se = (sub["hi"] - sub["lo"]) / (2 * 1.96)
    w = 1.0 / se ** 2
    m = float((w * sub["point"]).sum() / w.sum())
    se_p = float(np.sqrt(1.0 / w.sum()))
    lo, hi = m - 1.96 * se_p, m + 1.96 * se_p
    Q = float((w * (sub["point"] - m) ** 2).sum())
    dfree = len(sub) - 1
    I2 = max(0.0, (Q - dfree) / Q) * 100 if Q > 0 else 0.0
    return {
        "k": int(len(sub)),
        "pooled": round(m, 4),
        "ci": [round(lo, 4), round(hi, 4)],
        "I2": round(I2, 4),
        "tost055_p": float(stats.norm.cdf((m - 0.55) / se_p)),
        "tost060_p": float(stats.norm.cdf((m - 0.60) / se_p)),
    }


def main():
    out = {"ALL k=11": pool(cells), "pCR-only k=9": pool(cells[~is_surv])}
    with open(R / "meta_pools.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
