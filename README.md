# RankBased Drug-Sensitivity Transfer — code & data

Rank-based, cohort-free pharmacogenomic-signature method (breast cancer).

## Layout
- `scripts/` — signature construction, cohort-free singscore deployment, cohort
  download/parse, PAM50, axis-decomposition test, GATE check, 2-gene MKI67-ESR1
  baseline, multi-drug, GSEA, ablation, metrics. Sub-folders `phase_2` (GSE32646),
  `phase_4` (METABRIC), `phase_6/7/8` (nine-signature replication), `phase_D`
  (bortezomib / multiple myeloma).
- `figures/` — figure generators.
- `frozen_signatures/` — **data**: frozen bidirectional cell-line drug signatures
  (4 strategies) + `MANIFEST.json` (sha256). `scripts/phase_D/*.tsv` = bortezomib signatures.

Patient cohorts (GEO) and CCLE/GDSC/CTRP/OncoKB are public; download + parse is scripted.
Python 3.11 (numpy, pandas, scipy, scikit-learn, lifelines, matplotlib). Seeds pinned.
