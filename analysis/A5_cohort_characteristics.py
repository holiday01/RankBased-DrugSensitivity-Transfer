#!/usr/bin/env python3
"""A5 -- cohort characteristics table (Reviewer 3 point 3).

R3-3 asks how many ER+, HER2+ and TNBC patients each cohort contains, because without
those counts a reader cannot judge whether the results are driven by subtype imbalance.
Pure tabulation from deposited phenotype files; no model is fitted and no outcome is
predicted.
"""
import json, os
from pathlib import Path
import pandas as pd

_DR = Path(os.environ.get("DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
_OR = Path(os.environ.get("OUT_ROOT", Path(__file__).resolve().parents[1] / "results"))
OUT = _OR / "A5"; AXES = ["P","H","B","H_PAM50basal"]

SRC = {
 "GSE16446": (_DR/"external_data/GSE16446/GSE16446_pcr.tsv", _DR/"resubmission_v2/results/phase_0/axes_GSE16446.tsv"),
 "GSE22226": (_DR/"external_data/GSE22226/GSE22226_pcr.tsv", _DR/"resubmission_v2/results/phase_0/axes_GSE22226.tsv"),
 "GSE25066": (_DR/"external_data/GSE25066/GSE25066_pcr.tsv", _DR/"resubmission_v2/results/phase_0/axes_GSE25066.tsv"),
}

def pos(x): return x.astype(str).str.lower().str.startswith(("p","1","y","t"))

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows=[]
    for coh,(pp,ap) in SRC.items():
        pcr=pd.read_csv(pp,sep="\t").set_index("patient_id")
        ax=pd.read_csv(ap,sep="\t",index_col=0).dropna(subset=AXES+["pCR"])
        common=pcr.index.intersection(ax.index)
        d=pcr.loc[common]; y=ax.loc[common,"pCR"].astype(int)
        er,her2=pos(d["er_status"]),pos(d["her2_status"])
        tnbc=(~er)&(~her2)
        rows.append(dict(cohort=coh, n_profiled=len(pcr), n_analysed=len(common),
            events=int(y.sum()), pcr_rate=round(float(y.mean()),3),
            ER_pos=int(er.sum()), ER_neg=int((~er).sum()),
            HER2_pos=int(her2.sum()), HER2_neg=int((~her2).sum()),
            TNBC=int(tnbc.sum()), HRpos_HER2neg=int((er&~her2).sum())))
    df=pd.DataFrame(rows); df.to_csv(OUT/"A5_cohort_characteristics.tsv",sep="\t",index=False)
    json.dump(rows,open(OUT/"A5_cohort_characteristics.json","w"),indent=2)
    print(df.to_string(index=False))
    print("\nnear-degenerate axes (a contrast that is almost absent):")
    for r in rows:
        if r["HER2_pos"]<=10: print(f"  {r['cohort']}: only {r['HER2_pos']} HER2-positive of {r['n_analysed']} -> B axis near-constant")
        if r["ER_pos"]<=5:    print(f"  {r['cohort']}: only {r['ER_pos']} ER-positive of {r['n_analysed']} -> H axis near-constant")
    print(f"\nwrote {OUT}/A5_cohort_characteristics.tsv")

if __name__=="__main__": main()
