# Pre-Registration: External Validation of CCLE-Derived Drug-Sensitivity Signatures

**Project**: BioAdv resubmission of BIOADV-2026-206
**Date**: 2026-05-28
**Status**: DRAFT (sha256 fingerprints + final commit pending Stage 0 completion)
**Public location** (after final lock): OSF / Zenodo URL TBD

---

## 0. Scope of clinical claim (locked)

This framework is registered as a **research-grade biomarker** (AUROC tier
0.55–0.80) for **predicting pathologic complete response (pCR) under
anthracycline-class neoadjuvant regimens**. It is NOT a clinical decision-support
tool and is NOT a pan-drug intrinsic sensitivity profiler. Specifically:

- AUROC ≥ 0.80 with prospective calibration would be required for trial-candidate
  status (we explicitly do not target this).
- AUROC ≥ 0.85 with prospective validation and decision-curve net benefit would
  be required for clinical decision support (we explicitly do not target this).
- On lock of this pre-registration, the following manuscript sentences/clauses
  are deleted or rewritten, bound by quoted text (not line numbers, which are
  fragile across edits):
  - DELETE: "could support treatment selection decisions" (Discussion)
  - DELETE: "compares favorably with AUROC values reported for published domain
    adaptation methods" along with the "exceeded the AUROC values reported for
    CODE-AE, TRANSACT, and TUGDA" claim (Results / Discussion)
  - REWRITE: "AUROC of 0.780 [doxorubicin]" → must be restated as a TCGA-only
    observation with explicit GATE caveat; not a clinical performance claim
  - REWRITE: Abstract AUROC 0.752 headline → external-cohort AUROC framing
    after Stage 5 results are unblinded

The narrowed claim per hypothesis:
- H1/H2/H3 are tests of "predictor of pCR under anthracycline-containing regimen",
  not "intrinsic drug sensitivity".

## Purpose

This document pre-registers the EXACT primary analyses for external validation of
the bidirectional drug-class signature framework, BEFORE downloading any external
cohort data. Any analysis not specified here is exploratory and cannot serve as a
headline result.

The pre-registration is necessary because:
1. The prior submission's headline AUROC 0.752 was retracted as a base-rate
   artifact (GATE diagnostic, 2026-05-19).
2. A reviewer's central concern was that with 30 drugs × 3 strategies × multiple
   metrics × multiple ablations, "somewhere will produce a good number." This
   document forecloses that route by fixing the primary analyses in advance.

This document references frozen signature artifacts produced by
`scripts/freeze_signatures.py` and verified by `scripts/verify_frozen.py`. The
exact sha256 of each artifact is listed in `frozen_signatures/MANIFEST.json` and
fingerprinted into Section 9 below before final lock.

---

## 1. Primary hypothesis

**H1 (anthracycline-class transferability, GSE16446):**
The frozen Doxorubicin signature (treated as anthracycline-class proxy) — derived
from solid-tumor-only CCLE cell lines using strategy `quantile_30_oncokb_all` —
applied via full-transcriptome rank-based singscore to GSE16446 (single-arm
epirubicin, ER−/TN BC, pCR endpoint) yields AUROC > 0.55 for predicting pCR with
95% non-parametric bootstrap CI lower bound > 0.50.

**Success criterion (two-tier, locked)**:
- **PRIMARY (clinical-mechanism evidence)**: PAM50-basal-adjusted logistic
  regression `glm(pcr ~ signature_score + pam50_basal_score, family=binomial)`
  yields OR per SD of signature_score ≥ 1.3 with Wald p < 0.0083 (Bonferroni
  α = 0.05/6 primary tests, see §6).
- **SUPPORTING (discrimination evidence)**: raw AUROC > 0.55 with 95% bootstrap
  CI lower bound > 0.50.

Note: at expected effect size ≈ 0.55 and n ≈ 120, the SUPPORTING criterion is
near-impossible without inflated luck. PRIMARY is the headline; SUPPORTING is
reported alongside.

**Anthracycline-class proxy disclosure** (locked): Epirubicin (GSE16446 agent)
and Doxorubicin (CCLE signature source) share TOP2-poisoning mechanism but
differ pharmacokinetically (cardiotoxicity profile, optimal dose). Signature
transfer is **class-level**, not drug-identical. The same disclosure is added
to manuscript Methods.

**Anthracycline-proxy paired test (H3b, EXPLORATORY — not in α=0.0083 family)**:
GSE22226 patients receive AC-T (anthracycline + taxane). Within these patients,
the Doxorubicin signature should outperform the Paclitaxel signature by ≥ 0.05
AUROC for pCR (paired bootstrap 95% CI).
- Statistical caveat (per Audit #4): with n≈150 and pCR ~25%, paired-bootstrap
  ΔAUROC CI half-width is ~0.08–0.10 — clearing the 0.05 threshold is unlikely
  without strong mechanism-specific signal. H3b is explicitly framed as
  EXPLORATORY mechanism evidence; failure does NOT retract H3 primary, and
  success is reported alongside qualitative caveats rather than as a
  falsification of taxane mechanism in this cohort.

**Rationale for choosing solid-tumor pool** (per Round 2 audit 2026-05-28): The
pan-cancer Doxorubicin signature is dominated by hematopoietic-lineage markers
(IKZF1, LCP1, CXCR4, WAS, PTPRC), a known CCLE pan-cancer artifact. The
solid-tumor pool removes 149 heme cell lines (Myeloid + Lymphoid + Plasma Cell
lineages) at signature derivation time, surfacing the bona-fide TOP2-poison
sensitivity marker SLFN11 and resistance marker ABCC3 (anthracycline efflux).

**Rationale for choosing quantile_30_oncokb_all as primary strategy**: It was the
highest-performing single strategy in the original manuscript (AUROC 0.574),
uses the most biologically interpretable gene universe (OncoKB cancer genes),
and demonstrates the highest bootstrap-Jaccard stability among solid-only
strategies (median 0.463; DOXORUBICIN specifically: sens=0.579, res=0.714).

**Rationale for choosing single-drug Doxorubicin over TOP2-poison consensus**
(per Audit #4 of Round 2): The consensus signature (Dox + Etop rank-average)
loses anthracycline-specific markers (NQO1 quinone metabolism, AXL, DKK1, BIRC3)
that epirubicin (also an anthracycline) shares with Doxorubicin. The single-drug
signature is the better mechanism match for epirubicin; consensus is reported
as SECONDARY/sensitivity analysis only.

**Primary-eligibility gate (Round 2 audit recommendation):** A signature is
primary-eligible iff bootstrap Jaccard ≥ 0.40 in BOTH sensitivity and resistance
directions on its source pool. DOXORUBICIN solid_only/quantile_30_oncokb_all
clears this gate (sens=0.579, res=0.714).
Drugs failing this gate (CYCLOPHOSPHAMIDE 0.25, FULVESTRANT 0.25, DACARBAZINE
0.28) are exploratory only and not used in primary inference.

## 2. Secondary hypotheses

**H2 (cross-cohort robustness, GSE25066 TFAC):** Same signature scored on
GSE25066, stratified by clinical subtype:
- Primary stratum: TN (triple-negative); AUROC for pCR via primary scoring rule
  = anthracycline-only signature score
- Secondary stratum: HR+ (luminal); same scoring rule
- Combo regimen scoring rule pre-specified in Section 4 below

**H3 (prospective cohort confirmation, GSE22226 I-SPY1):** Same signature
applied to GSE22226 (AC-T neoadjuvant, ~150 patients, pCR + RCB endpoints):
- Primary endpoint: pCR AUROC (matched to H1/H2)
- Secondary endpoint: RCB index (continuous) — Spearman ρ with signature score

**H4 (negative-control sanity):** Wrong-mechanism signatures (Dasatinib,
Gefitinib, Erlotinib — Imatinib excluded due to Jaccard < 0.40) applied to
GSE16446 (anthracycline cohort) should yield AUROC indistinguishable from 0.5
(95% CI contains 0.5). If any of these reaches AUROC > 0.6, cohort-level
confound is suspected and H1 result is qualified.

## 2b. TCGA-side CV protocol (committed)

For any TCGA-based analysis retained in the manuscript (e.g., supplementary
reproduction of the original 0.553/0.613/0.752 with corrected Spearman sign):

- **GroupKFold by `patient_id`**, k=5, shuffle=True, random_state=42
- No patient appears in train + test
- All signature derivation done on CCLE only (no TCGA labels touched at fit time)
- Methods text adds: "Signatures were derived from CCLE drug-AUC contrasts
  using only cell-line data; TCGA samples were used solely for held-out scoring
  via within-sample rank singscore. No TCGA outcome label entered signature
  derivation."

## 3.1 PAM50 subtype computation (pinned, used in §5 + §6.2)

For GSE16446 (no a priori PAM50 metadata), compute the basal class membership
as follows:

- **Centroid source**: PAM50 centroid matrix from the `genefu` R package
  (Bioconductor; Parker et al. 2009 centroids). Specific release: `genefu`
  version recorded at scoring time in `score_external_cohort.py` output JSON
  sidecar.
- **Algorithm**: Spearman correlation of each patient's cohort-internal-rank
  expression with each of the 5 PAM50 centroids; assign patient to the class
  with the highest correlation. For permutation null #1 binarization (basal vs
  non-basal): basal = highest-correlation class is `Basal`; non-basal otherwise.
- **Tie-breaking**: deterministic by class-name alphabetical order; ties are
  expected to be rare (<1% per Parker 2009).
- **Gene-symbol alignment**: PAM50 centroid genes must map to ≥40 of 50 cohort
  HGNC symbols after alias resolution; if <40, GSE16446 stratification falls
  back to ER−/HER2− indicator (uniformly the cohort) which collapses the
  stratum permutation to single-class — in that fallback case Null #1 is
  reported as `not_applicable` and only Null #2 (gene-set) counts toward the
  primary decision rule for H1.

(For GSE25066 and GSE22226, intrinsic subtype is taken from cohort metadata,
no PAM50 recomputation required.)

## 3. Scoring methodology

For all hypotheses:
- **Full-transcriptome rank-based singscore** (per manuscript Methods Eq. 1).
- Rank computed over the entire external cohort expression matrix
  (post probe-to-symbol collapse), NOT over the signature gene subset.
- Up-genes (sensitivity) and down-genes (resistance) drawn from frozen artifact;
  no re-derivation, no re-fitting.
- Probe-to-symbol mapping: max-IQR probe per HGNC symbol; HGNC alias resolution
  applied before signature lookup. Alias table sha256 pinned in Section 9.
- No FSQN normalization at primary scoring (D1: rank-only primary; FSQN as
  parallel sensitivity comparator).

## 4. GSE25066 combination regimen scoring rule (committed before unblinding)

GSE25066 patients received TFAC = Paclitaxel + 5FU + Doxorubicin + Cyclophosphamide.
Docetaxel was considered but EXCLUDED for low stability (Jaccard sens=0.395 <
0.40 in solid_only pool).

- **Primary score for pCR prediction**: frozen Doxorubicin signature score
  (anthracycline-class proxy; PRIMARY-eligible).
- **Secondary score**: frozen Paclitaxel signature score (taxane proxy;
  PRIMARY-eligible).
- **Exploratory composite**: max(Doxorubicin_score, Paclitaxel_score) per patient.
- **NOT permitted as primary**: arithmetic average of 5 component signatures;
  per-drug AUROC cherry-picking; Docetaxel-based scoring; Cyclophosphamide/5FU
  signatures (both Jaccard < 0.40).

## 5. Stratification rules

- **GSE16446 (PAM50/basal covariate, per Round 2 audit):** The solid-only
  signature's resistance axis is enriched for EMT genes (CAV1, AXL, MET, EPCAM,
  CDH1), which correlate with basal/claudin-low intrinsic subtype. To prevent
  spurious correlation with subtype, the primary analysis adjusts for
  PAM50-derived basal score (or, if PAM50 metadata unavailable, the cohort's
  ER−/HER2− status indicator) as a covariate in a logistic regression rather
  than reporting raw AUROC. Both adjusted and unadjusted AUROC are reported.
- **GSE25066:** stratified primary analysis by intrinsic subtype (TN, HR+,
  HER2+). Pooled AUROC reported only as a secondary metric with explicit
  base-rate caveat.
- **GSE22226 (I-SPY1):** stratified by HR/HER2 status if metadata permits;
  otherwise pooled with caveat.

## 6. Statistical inference

### 6.1 Estimation

- **Bootstrap 95% CI**: non-parametric, **1000 resamples**, fixed bootstrap seed
  **= 42**. (Sensitivity comparison with N=5000 reported as Supplementary if any
  primary test's CI lower bound falls in [0.49, 0.52], otherwise N=1000 suffices.)
- **Calibration metrics (R1-M5, R2-2)**: Brier score, slope+intercept from
  isotonic regression of `predicted_pCR_probability ~ true_pCR_label`, and a
  10-bin reliability curve. Reported PER stratum where stratification applies.
- **Class-imbalance-appropriate discrimination (R1-M5, R2-2)**: PR-AUC
  (average-precision) and F1 at Youden-J threshold reported alongside ROC-AUC
  for every primary test. PR-AUC reported with bootstrap 95% CI (same n=1000
  resamples, same seed=42). Per Reviewer 2's R2-2 request, F1 + precision +
  recall at the Youden-J threshold are reported in the primary results table,
  not relegated to supplementary.

### 6.2 Permutation nulls

- **Permutation null #1 (within-stratum label permutation)**: Shuffle pCR labels
  WITHIN each pre-specified stratum, recompute AUROC, 1000× with **seed = 1337**.
  Empirical p = (perm_AUROCs ≥ observed) / 1001.
  - Stratum variable per cohort:
    - **GSE16446**: PAM50 basal vs non-basal class (computed from cohort
      expression matrix using `genefu`/PAM50 centroid correlation; algorithm
      pinned in §3.1).
    - **GSE25066**: clinical intrinsic subtype (TN, HR+, HER2+) from cohort
      metadata.
    - **GSE22226**: HR/HER2 status from cohort metadata.
  - The implementation refuses to fall back to global label-shuffle (script
    `score_external_cohort.py` raises ValueError if strata is None).
- **Permutation null #2 (gene-set permutation)**: Replace signature with random
  same-size gene set from the same gene universe (OncoKB-all / OncoKB-onco-only /
  protein-coding as appropriate per strategy; sha256 of universe lists in §9),
  recompute AUROC, 1000× with **seed = 1338**.
- **Decision rule**: BOTH nulls must pass at Bonferroni-adjusted α = 0.0083 (see
  §6.3) for the primary hypothesis to be confirmed.

### 6.3 Multiple testing — locked test family

The primary-claim test family contains 6 tests:

| # | Test | α/6 = 0.0083 |
|---|---|---|
| T1 | H1 (GSE16446, PAM50-adjusted OR ≥ 1.3, Wald p) | ✓ |
| T2 | H2 (GSE25066, TN stratum, OR ≥ 1.3, p) | ✓ |
| T3 | H2 (GSE25066, HR+ stratum, OR ≥ 1.3, p) | ✓ |
| T4 | H3 (GSE22226, pCR endpoint, OR ≥ 1.3, p) | ✓ |
| T5 | H3 (GSE22226, RCB index, Spearman ρ p) | ✓ |
| T6 | H4 (GSE16446, max wrong-MoA AUROC < 0.6) | ✓ |

H4 is interpreted as a *constraint* (max ≥ 0.6 invalidates H1), not a hypothesis;
listed in the family for completeness. If any stratum has n_responders < 5 OR
n_non_responders < 5, that test is dropped from the family pre-unblinding (the
α threshold tightens accordingly).

Secondary / exploratory tests (FSQN sensitivity, consensus signature,
quantile_20_protein, spearman_top20_oncokb, H3b paired anthracycline-vs-taxane
test, ablation grid, calibration metrics, DCA) are NOT in the α=0.0083 family
and are reported with point estimates + 95% CIs only.

### 6.4 Decision-Curve Analysis (DCA, clinical utility — H1 QUALIFIER)

Decision-Curve Analysis (Vickers net-benefit) is computed PER primary
hypothesis. Net-benefit curve over threshold range 5–40%
(predicted-pCR-probability) is reported.

**H1 outcome levels (per Audit #4 of Round 4):**
- **CONFIRMED**: PRIMARY criterion (adjusted OR ≥ 1.3, p < 0.0083) PASSES AND
  DCA shows net benefit > "treat-all" at the cohort's empirical pCR base rate.
- **SUPPORTED-ONLY**: PRIMARY criterion PASSES BUT DCA fails. Headline is
  framed as "discrimination evidence; clinical-utility net benefit not yet
  demonstrated."
- **NOT SUPPORTED**: PRIMARY criterion fails (DCA irrelevant; H1 retracts).

This 3-state outcome directly answers R2-2's "0.75 not clinically useful" by
making clinical net benefit an explicit qualifier on the headline framing.

## 7. Comparison to existing methods (same-cohort head-to-head)

For H1 (GSE16446) primary cohort:
- **pRRophetic (epirubicin if available; doxorubicin otherwise)**:
  AUROC + bootstrap CI on the SAME GSE16446 patients.
- **oncoPredict (GDSC2 doxorubicin model)**: same.
- **CODE-AE pretrained on GSE16446** — committed as a REQUIRED comparator
  (audit-driven upgrade from "1-day budget" framing). If the pretrained
  checkpoint or inference pipeline is not reproducible within 3 calendar days,
  the response letter explicitly states: "we attempted CODE-AE reproduction
  on GSE16446; the public release lacks XYZ component required for inference
  in our setting; we therefore acknowledge this as a non-comparable method and
  remove the cross-study claim from contribution." This explicit-failure
  framing replaces silent omission.
- **TRANSACT / TUGDA**: not run; explicit rationale in response letter — these
  methods require unlabeled target-cohort transcriptomes during training, which
  contradicts our cohort-free deployment claim; running them would test a
  different methodology, not our framework. This is documented as out-of-scope,
  not retreat.

No method comparison is required for H2/H3 — those are cross-cohort robustness
checks of our framework only.

## 8. Trivial baselines (R1-M4)

On all three external cohorts, the framework AUROC will be compared to four
trivial baselines:
- **Random scoring** (uniform random per patient; expected AUROC ≈ 0.5).
- **Single-gene MKI67 score**: generic proliferation marker.
- **Hallmark E2F_TARGETS singscore**: curated proliferation pathway score.
- **Hallmark G2M_CHECKPOINT singscore**: same.

**Success criteria (revised per Audit #2 of Round 3)** — both required:
- (a) **Paired DeLong test**: framework AUROC > each baseline individually, with
  paired DeLong p < 0.0083 (Bonferroni-adjusted) in ≥ 2 of 3 cohorts.
- (b) **Δ-AUROC**: framework AUROC − best baseline AUROC ≥ +0.03 (relaxed from
  +0.05 because expected effect size is ~0.55 and CI half-width is ~0.10 — the
  +0.05 + non-overlapping-CIs combination was self-defeating).

"Best baseline" is the *single pre-specified* proliferation baseline with the
highest empirical AUROC in the cohort (no max-over-4 selection bias — the
specific baseline winning per cohort is reported alongside).

## 8b. Ablation robustness grid (R1-O3, locked)

After primary tests complete (and only as a robustness check, not for primary
inference), report the full grid in a Supplementary Table:

- Quantile threshold q ∈ {0.20, 0.25, 0.30}
- Top-N per direction ∈ {10, 20, 30, 50}
- Gene universe ∈ {OncoKB-all, OncoKB-onco-only, protein-coding}
- Up:down direction asymmetry ∈ {30:30, 30:10, 10:30}
- Significance gate ∈ {none, FDR<0.10 on per-gene t-test, FDR<0.05}
- Selected baseline: the pre-registered primary config
  (q=0.30, top-N=30, OncoKB-all, 30:30, no FDR gate)

**Robustness criterion**: the pre-registered primary config must rank in the
**top quartile of AUROC across all grid configs** on each cohort; if not, the
pre-registered choice was non-robust and the Discussion explicitly notes this
(no headline claim is changed).

## 9. Frozen signature artifacts (sha256 fingerprints, locked 2026-05-28)

All signatures live under
`revise_bioadv/frozen_signatures/solid_only/quantile_30_oncokb_all/` unless
noted. SHA256 abbreviated to first 16 hex chars for readability; full hashes
in `MANIFEST.json` and verifiable by `scripts/verify_frozen.py`.

### Primary (H1, H2, H3)
- **PRIMARY anthracycline-class proxy** — `DOXORUBICIN.tsv`
  - sha256: `dc53c93a73503e7e...`
  - Jaccard: sens=0.579, res=0.714 (PRIMARY-ELIGIBLE)

### Secondary / sensitivity analyses
- TOP2-poison consensus (Dox+Etop rank-average) — `TOP2_POISON_CONSENSUS.tsv`
  - sha256: `bd2bbfc6dca15d70...`
  - Jaccard: sens=0.500, res=0.593 (eligible; SECONDARY by design — anthracycline-specific dilution per Audit #4)
- Taxane proxy (H2 GSE25066 secondary stratum) — `PACLITAXEL.tsv`
  - sha256: `ed7fa3734ef8831f...`
  - Jaccard: sens=0.500, res=0.538 (eligible)
- TOP2 reference (sensitivity analysis) — `ETOPOSIDE.tsv`
  - sha256: `0edf7673b1c76604...`
  - Jaccard: sens=0.538, res=0.622 (eligible)

### Negative controls (H4, GSE16446 wrong-MoA tests)
- `DASATINIB.tsv` sha256: `57938ff6fed90990...` (sens=0.714, res=0.579)
- `GEFITINIB.tsv` sha256: `03c2b654ef7803d5...` (sens=0.714, res=0.714)
- `ERLOTINIB.tsv` sha256: `016e654b5e890afe...` (sens=0.765, res=0.690)
- IMATINIB EXCLUDED from negative controls (sens=0.304 fails Jaccard ≥ 0.40)

### Scripts
- `revise_bioadv/scripts/freeze_signatures.py`
  - sha256: `cc6cde39b4976fd1...`
- `revise_bioadv/scripts/verify_frozen.py`
  - sha256: `c4c11b7a2224de91...`

### Manifest
- `revise_bioadv/frozen_signatures/MANIFEST.json`
  - (sha256 captured in `preregistration.LOCK` companion file at lock time)

### Auxiliary scripts (also fingerprinted in MANIFEST.auxiliary_artifacts)
- `revise_bioadv/scripts/freeze_signatures.py`
- `revise_bioadv/scripts/verify_frozen.py`
- `revise_bioadv/scripts/score_external_cohort.py`
- `revise_bioadv/scripts/parse_gse_to_matrix.py`
- `revise_bioadv/scripts/download_gse.py`

### Gene universes (referenced by Permutation null #2)
- OncoKB-all (1202 genes after CCLE intersection) — sha256 derived from
  `OncoKB/cancerGeneList.tsv` (in MANIFEST.source_files.oncokb)
- OncoKB-onco-only (417 genes) — same source
- Protein-coding (19215 genes) — sha256 derived from
  `CCLE25Q3/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
  (in MANIFEST.source_files.ccle_expression)

Any artifact whose full sha256 differs from `MANIFEST.json` invalidates the
pre-registration. Run `python scripts/verify_frozen.py --archive-only` to
confirm bundled artifacts; `--strict` additionally verifies external CCLE/OncoKB
sources match the recorded hashes at freeze time.

## 10. What this pre-registration does NOT cover (exploratory)

The following are explicitly exploratory and cannot be claimed as primary
evidence regardless of outcome:
- FSQN-aligned scoring (sensitivity analysis only, parallel-reported to rank-only)
- TOP2-poison consensus signature (sensitivity analysis — supplements H1; not
  primary because it dilutes anthracycline-specific features per Audit #4)
- Pan-cancer pool signatures (acknowledged in Discussion as the prior version
  retracted due to heme confound)
- All ablation studies (B5; including up:down direction asymmetry and
  FDR-thresholded selection)
- Per-drug AUROC ranking
- Strata with n_responders < 5 OR n_non_responders < 5
- Quantile_20_protein, spearman_top20_oncokb strategies for primary external
  validation (used only in robustness reporting)
- Low-stability drugs (Jaccard sens_median OR res_median < 0.40 on source pool):
  CYCLOPHOSPHAMIDE, FULVESTRANT, DACARBAZINE on solid_only/quantile_30_oncokb_all
  pool — see MANIFEST.json `stability` field for the canonical list.

## 10b. Explicit primary-eligibility table (locked 2026-05-28)

Threshold: signature is PRIMARY-eligible iff sens_median ≥ 0.40 AND res_median
≥ 0.40 (per Audit #3 of Round 2). Pool = solid_only/quantile_30_oncokb_all
unless noted.

| Drug                  | Sens-Jaccard | Res-Jaccard | Eligible? | Role |
|---|---|---|---|---|
| DOXORUBICIN           | 0.579        | 0.714       | ✓         | **PRIMARY (H1, H2 anthracycline, H3 anthracycline)** |
| PACLITAXEL            | 0.500        | 0.538       | ✓         | H2 GSE25066 taxane secondary; H3 GSE22226 taxane |
| ETOPOSIDE             | 0.538        | 0.622       | ✓         | Sensitivity analysis (TOP2 reference) |
| TOP2_POISON_CONSENSUS | 0.500        | 0.593       | ✓         | SECONDARY (Audit #4: not primary — anthracycline dilution) |
| DASATINIB             | 0.714        | 0.579       | ✓         | Negative control (wrong-MoA TKI) |
| GEFITINIB             | 0.714        | 0.714       | ✓         | Negative control (wrong-MoA EGFR TKI) |
| ERLOTINIB             | 0.765        | 0.690       | ✓         | Negative control (wrong-MoA EGFR TKI) |
| DOCETAXEL             | 0.395        | 0.463       | **✗**     | Excluded (sens<0.40); not used for primary |
| IMATINIB              | 0.304        | 0.395       | **✗**     | Excluded — **both** sens AND res fail the ≥ 0.40 gate. Not used as negative control. |
| CYCLOPHOSPHAMIDE      | (low)        | (low)       | **✗**     | Exploratory only |
| FULVESTRANT           | (low)        | (low)       | **✗**     | Exploratory only |
| DACARBAZINE           | (low)        | (low)       | **✗**     | Exploratory only |

Adjustment from prior plan: GSE25066 taxane proxy uses Paclitaxel only (Docetaxel
excluded for instability). GSE16446 negative-control panel reduced from 4 to 3 TKIs
(Imatinib excluded). Both adjustments documented here, not post-hoc.

## 11. Pre-registration commitment

Once Section 9 is filled with sha256 values from a successful run of
`freeze_signatures.py` and confirmed by `verify_frozen.py`, and Section 11.1
below is filled with the self-sha256 + timestamp, this document is locked.
Any deviation in the primary analyses requires:
1. Explicit written rationale logged in `decisions.md`.
2. Re-classification of the affected hypothesis as exploratory.
3. Disclosure in the manuscript and response letter.

### 11.1 Self-lock

Lock timestamp (UTC): captured in companion file `preregistration.LOCK`
Lock sha256 of this document: captured in companion file `preregistration.LOCK`

(Run `scripts/lock_preregistration.py` after final edits; this generates the
LOCK file with this document's sha256, the MANIFEST.json sha256, and the UTC
timestamp. Upload BOTH files together to Zenodo / OSF as a single deposit so
the third-party-issued DOI + deposit timestamp constitute the immutable
external witness.)
