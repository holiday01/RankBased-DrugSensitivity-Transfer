# Cell-line drug-response and published breast-cancer transcriptomic signatures, benchmarked against four routine clinical axes

Code, result summaries and pre-registration for a benchmark study asking whether
cell-line-derived drug-response signatures and published breast-cancer transcriptomic
signatures carry neoadjuvant pathological-complete-response information beyond four
routine clinical axes (proliferation, ER/PR, HER2, PAM50-basal) in breast cancer.

**Cite:** 10.5281/zenodo.20726402 — the concept DOI, which resolves to the current version.
**Author:** Yen-Jung Chiu (ORCID 0000-0002-2087-2266), Chang Gung University /
Chang Gung Memorial Hospital, Taiwan. **License:** CC-BY-4.0.

## Layout

```
analysis/          the A0-A14 analyses, plus corrected_figdata/
code/
  framework/       cohort download and parse, PAM50, signature freeze and its
                   verifier, and the six pre-registered primary tests
  axis_analysis/   the scripts defining the published signatures the analyses
                   re-score, the two remaining cohort downloaders, and phase_D
results/           result summaries, JSON/TSV, one directory per analysis
figures/           supplementary Figures S1 and S2
frozen_signatures/ MANIFEST.json (sha256) and the signature gene lists used here
external_data/     no data — which cohort files each script needs, and where to get them
MANIFEST.txt       sha256 of every file in this deposit
```

## Running the analyses

Inputs and outputs resolve through two environment variables:

```bash
export DATA_ROOT="$PWD" OUT_ROOT="$PWD/results"
```

Cohort data is not distributed here. `external_data/README.md` lists what each script
needs and which script downloads it.

1. `code/framework/` — freeze signatures, lock the pre-registration, download and parse
   the cohorts, PAM50, then the primary tests T1-T6.
2. `analysis/` — A0 through A14 and the decision curve, each writing to the matching
   `results/` directory. The derived tables they read — the per-cohort clinical axes,
   the pCR and phenotype tables, the meta cells — are deposited under `results/`, so
   this step runs without re-deriving them.

`results/phase_9_revalidation/` carries `revalidation.py` and `compute_meta_pools.py`,
which produce the nested-LRT, meta-pool and TOST results the supplement reports.

Environment: Python 3.11 (numpy, pandas, scipy, scikit-learn, statsmodels, lifelines,
matplotlib; genefu/R for PAM50). Seeds pinned: bootstrap 42, permutation 1337/1338.

## Integrity checks

`python3 code/framework/verify_frozen.py` re-hashes each signature against
`frozen_signatures/MANIFEST.json` and reports a mismatch if any differs. The freeze
produced 32 drug signatures per pool; this deposit carries the 14 the analyses use,
so the verifier also reports the rest as absent. `MANIFEST.txt` carries the sha256 of
every file in this deposit.

## Paths quoted in the paper

Most resolve literally from this directory. Three do not:

- `scripts/within_stratum_immune.py`, `scripts/within_stratum_cytolytic.py` and
  `scripts/lrt_power_and_delong.py` are in `code/axis_analysis/`.
- the download and eval scripts described as being in `results/phase_D/` are in
  `code/axis_analysis/phase_D/`; the results are in `results/phase_D/`.
- a few scripts give 10.5281/zenodo.20726403 in their docstrings, the version DOI of v1.0.0. The
  concept DOI is 10.5281/zenodo.20726402, above.

## Not included

- **Cohort data.** Every cohort belongs to the group that published it and is public;
  see `external_data/README.md`. Per-patient *scores* — one value per sample, keyed by
  public accession — are in `results/A0/score_*.tsv` and `results/phase_D/scores_*.tsv`.
- **Intermediate tables** from the earlier phases, except where the paper names them as
  a source. `code/` regenerates them from the cohorts.
- **The pre-registration**, its lock files and the deviation log, which accompany the
  article rather than this deposit.
- **The article**, which is available from the journal.

## Disclosures

Held-out results are equivalence / no-detectable-signal results bounded by cohort
power, not positive validations. The primary estimator (nested-model LRT) was adopted
in response to methodological review; the originally pre-specified residual-AUROC test
is retained as a concordant robustness column. Neither the residual-AUROC decomposition
nor the LRT was pre-registered — the hash-locked parent plan specifies a
PAM50-basal-adjusted logistic OR/SD Wald test and names no other estimator. PROBAST
risk of bias is self-rated HIGH on the analysis domain.
