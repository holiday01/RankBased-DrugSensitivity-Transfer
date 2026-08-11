# RankBased-DrugSensitivity-Transfer

Two related studies share this repository.

**The framework.** A rank-based method for transferring cell-line drug-response
signatures to patient tumours, without cohort-specific retraining. The notebooks at
the repository root (`0. data_preprocess.ipynb`, `1. sig_enrich.ipynb`,
`2. validation.ipynb`) and `CODE_MANIFEST.md` are that work.

**The benchmark.** A follow-up asking what those signatures actually add. Every
directory below is that study: it benchmarks cell-line drug-response signatures and
nine published breast-cancer transcriptomic signatures against four routine clinical
axes (proliferation, ER/PR, HER2, PAM50-basal) for neoadjuvant
pathological-complete-response prediction, across five cohorts.

Archived with a DOI at **10.5281/zenodo.20726402** (concept DOI; resolves to the
current version).

## Layout

```
analysis/          the A0-A14 analyses, plus corrected_figdata/
code/
  framework/       cohort download and parse, PAM50, signature freeze and its
                   verifier, the six pre-registered primary tests
  axis_analysis/   the scripts defining the published signatures the analyses
                   re-score, two cohort downloaders, and phase_D
results/           analysis outputs, and the derived tables the analyses read
figures/           supplementary Figures S1 and S2
frozen_signatures/ MANIFEST.json (sha256) and the signature gene lists used here
external_data/     no data — which cohort files each script needs, and where to get them
```

## Running the analyses

```bash
export DATA_ROOT="$PWD" OUT_ROOT="$PWD/results"
```

No cohort data is in this repository; every cohort is public and belongs to the group
that published it. `external_data/README.md` lists what each script needs and which
script downloads it.

The derived tables the analyses read — per-cohort clinical axes, pCR and phenotype
tables, meta cells — are under `results/`, so `analysis/` runs without re-deriving
them. Two cohort expression matrices (GSE32646, METABRIC) are read directly and must
be fetched first.

Environment: Python 3.11 (numpy, pandas, scipy, scikit-learn, statsmodels, lifelines,
matplotlib; genefu/R for PAM50). Seeds pinned: bootstrap 42, permutation 1337/1338.

## Which analysis answers which review point

Each script's docstring states this; in summary:

| analysis | question |
|---|---|
| `A2_grid.py` | the full signature × cohort added-value grid |
| `A3_baselines.py` | are four axes enough, or are stronger baselines needed |
| `A4_within_subtype.py` | does the signature predict inside a subtype, or only separate subtypes |
| `A5_cohort_characteristics.py` | ER+/HER2+/TNBC counts per cohort |
| `A6b_gse41998_endpoint.py` | does the verdict change under a harmonised endpoint |
| `A7b`, `A7c` | why 30 genes per direction; is the choice stable |
| `A10b_required_n_corrected.py` | what n would be needed to exclude a small effect |
| `A11` (in `A_remaining.py`) | multi-drug regimens versus a single-drug signature |
| `A13_native_concordance.py` | is a null the signature, or our re-implementation |
| `A14a/b/c` | a second independent cohort, GSE41998 |

## Integrity

`python3 code/framework/verify_frozen.py` re-hashes each signature against
`frozen_signatures/MANIFEST.json`. The freeze produced 32 drug signatures per pool;
this repository carries the 14 the analyses use, so the verifier reports the rest as
absent.

## Citation

Cite the paper for the science and the concept DOI for the code:
`https://doi.org/10.5281/zenodo.20726402`.

## License

See `LICENSE`.
