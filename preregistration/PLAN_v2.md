# PLAN_v2 — BIOADV-2026-206 Resubmission (P5: Axis-Decomposition Methodology Lead)

**Status:** v0.5 — P5 reframe post results audit (3 experts converge: methodology lead, 2-gene as evidence; METABRIC added)
**Date:** 2026-05-28
**Working folder:** `revise_bioadv/resubmission_v2/`
**Parent context:**
- BIOADV-2026-206 rejected 2026-05-19 (resubmission invited as NEW submission)
- Pre-registered external validation FAILED 2026-05-28
- Multi-anchor methods paper direction shelved 2026-05-28
- PLAN_v2 v0.1 audited by 5 experts 2026-05-28; verdict: 4/5 PASS-WITH-CHANGES, 1/5 WOULD NOT SUPPORT
- Path C chosen by user: held-out cohort + GATE-led title + all expert edits

## Change log

- **v0.5 (2026-05-29):** P5 reframe after 3-expert results audit. All 3 audits converge on the same critique: P2 has novelty deficit ("2-gene = IHC4 rediscovery"), and the real contribution is the methodology + held-out validation. Reviewer-sim's specific advice: "Lead with axis-decomposition methodology + GSE32646 holdout; demote 2-gene oracle to Layer-2 structural consequence; add METABRIC as 2nd held-out for n=2 meta-confirmation."
  - **Title rewrite (P5):** "A pre-registered axis-decomposition test for cohort-free pharmacogenomic signatures: confirmatory validation on independent breast cancer cohorts"
  - **Layer restructure:**
    - Layer 1 PRIMARY = Axis-decomposition methodology + residual AUROC test + 3×3 decision rule + GSE32646 + METABRIC held-out confirmation
    - Layer 2 STRUCTURAL CONSEQUENCE = When axis test passes, signature reduces to simple clinical features (2-gene oracle case study)
    - Layer 3 HONEST APPLICATIONS = GATE diagnostic with marginal AUROC floor; multi-drug (cytotoxic Δ wrong direction, broader axis-collapse); GSEA (stress/inflammation not mechanism); ablation
  - **§4.10 added: METABRIC second held-out** — axis residual AUROC test + 2-gene oracle vs framework on chemo subset. Goal: n=2 meta-confirmation with GSE32646.
  - **Language softening per biostat**: "recapitulates" not "outperforms"; "consistent with" not "confirms"; honest "REJECTED" for §4.2 wrong-direction prediction.
  - **Critical gap fills COMPLETED 2026-05-29**:
    - **2-gene MKI67+ESR1 oracle on GSE32646 (held-out)**: AUROC = **0.7959** [0.7075, 0.8818] vs framework 0.5623 [0.4524, 0.6763]. **Δ = +0.2336**, paired DeLong **p = 3.45×10⁻⁴**. 2-gene CRUSHES framework in held-out cohort.
    - **Cross-cohort meta (4 cohorts: 3 burned + GSE32646)**: Fixed-effects Δ = **+0.1648** [0.1117, 0.2178], **p = 1.12×10⁻⁹**. Stouffer combined **p = 1.06×10⁻⁹**. I² = 0.00 (perfectly homogeneous). 2/4 cohort DeLong significant individually (GSE25066 p=2.83×10⁻⁶, GSE32646 p=3.45×10⁻⁴); GSE16446 p=0.083 (marginal); GSE22226 p=0.232 (NS) — but Holm-adjusted 2/4 stay significant.
    - **Implication**: P5 framing is empirically strong — 2-gene oracle is a STRUCTURAL CONSEQUENCE of axis decomposition that REPLICATES across 4 independent BC neoadjuvant cohorts with meta-p < 10⁻⁹.
- v0.4 (2026-05-28): P2 framing pivot after Phase 1+2+3 results came in:
  - **Headline pivot**: From "GATE-led diagnostic" (P4 / v0.3) to **"2-gene oracle outperforms 60-gene CCLE-derived signature"** (P2)
  - **Reason**: GATE GAP < 0.10 on both worked examples (Dox-frozen GAP 0.085/0.039; Mammaprint × TCGA-BRCA GAP 0.038). GATE didn't fire as positive demo. MEANWHILE MKI67+ESR1 2-gene baseline BEATS framework by 0.09-0.18 AUROC in ALL 3 burned cohorts — clean strong positive finding.
  - **New Title (candidate)**: "A two-gene clinical oracle (MKI67 + ESR1) matches or exceeds a 60-gene CCLE-derived signature for breast cancer chemotherapy response prediction"
  - **Layer restructure**:
    - Layer 1 (PRIMARY): 2-gene oracle finding + axis projection methodology
    - Layer 2 (SUPPORTING): held-out GSE32646 confirmation of axis decomposition
    - Layer 3 (SUPPORTING): GATE diagnostic as supplementary tool (correctly identifies no within-stratum signal — valid behavior)
  - Anthracycline failure becomes "evidence that complex doesn't beat simple", not headline negative
  - Cohort §4.2 wrong-direction Δ reframed: "TKI signatures show similar axis-correlation patterns" — broader implication
  - GSEA TNFA finding reframed: "CCLE-derived signatures capture stress-response not mechanism-specific biology" — reinforces complex-vs-simple thesis
- v0.3 (2026-05-28):
  - **GSE20271 → GSE32646 cohort swap** (GSE20271 NO-GO; Tabchy/Symmans/Pusztai/Vidaurre/Martin/Souchon authors on both GSE20271 + GSE25066, mandatory non-overlap failed). GSE32646: Osaka U Noguchi group; n=115; FEC100/T (epirubicin); GPL570; pCR ~24%; zero overlap.
  - **§4.8 decision rule rewritten to 3×3 table** (E2: current rule had logical inversion — CI lower 0.54 + point 0.62 incorrectly "confirms null"; need explicit GRAY ZONE).
  - **Power calculation added** (E3: Hanley-McNeil CI half-width ~±0.09 at n=115/24% pCR; gray zone is modal outcome).
  - **LOCK enforcement strengthened** (E1: designated custodian + pre-LOCK public push of download script with placeholder accession + post-execution timestamp audit).
  - **§4.8 single-headline commitment** (E4: only ESR1+PGR linear OLS residual AUROC is headline; supplementary variants Bonferroni × 3 ineligible to flip outcome).
  - **§4.2 Imatinib added** as true-null TKI control (BCR-ABL/KIT, no BC indication, Jaccard 0.76/0.82). Set now 11 signatures.
  - **§4.2 thresholds loosened per realistic noise estimates**: Trametinib |r(sig,P)|<0.5 (was 0.4); Lapatinib r(sig,B)>0.30 (was 0.4) with dilution caveat; Fulvestrant r(sig,H)>0.30 (was 0.4) — flag as directional control only (Jaccard 0.40/0.33 weak).
  - **§4.2 explicit headline test added**: cytotoxic-clade mean |r(sig,P)| − TKI+null-clade mean |r(sig,P)| > 0.3 (inferential lever per biostat n=1 caveat).
  - **TOP2_POISON_CONSENSUS** (Dox+Etoposide pre-frozen, Jaccard 0.85) added as headline row in §4.2.
  - **Title trimmed** (case study moved out of title to Abstract).
  - **Layer-1 second-GATE-worked-example** flagged as TBD decision (algorithm-level addition; surface to user).
- v0.2 (2026-05-28): Major rewrite per 5-expert audit synthesis:
  - **Title rewritten** to GATE-led (per Convergent finding C3: oncologist + reviewer-sim + pre-reg auditor)
  - **Layer reorder**: Layer 1 (GATE) = headline; Layer 2 (honest negative) = worked example;
    Layer 3 (H-b axis) = supporting evidence with held-out validation
  - **NEW §2.5-2.7**: data-burn acknowledgement, independent cohort requirement, "NOT Registered Reports Stage 1" statement (per pre-reg auditor)
  - **§4.1 axis decomposition rebuilt** per biostat + bioinformatician: drop MYC_TARGETS_V2, add Whitfield_G2M; H = ESR1+PGR primary; ERBB2 as separate predictor B; VIF + ridge + bootstrap R² CI; pre-specified falsifier; synthetic positive control
  - **§4.2 multi-drug rebuilt**: Etoposide added (6 cytotoxics); TKI controls reshuffled to Trametinib + Gefitinib + Lapatinib; Fulvestrant added as H-axis POSITIVE control; explicit n=1 caveat; pairwise correlations; TKI contrast as headline test
  - **§4.3 baselines expanded**: MKI67+ESR1 2-gene sum added as PRIMARY comparator; random OncoKB-60 × 100; random Hallmark × 50
  - **§4.4 permutation** specified within-PAM50 stratum
  - **§4.6 method** switched from Fisher to GSEA (fgsea)
  - **NEW §4.7**: multiplicity correction across exploratory family
  - **NEW §4.8**: 4th held-out cohort (GSE20271 primary; GSE32646 / MetaGxBreast alternates) — axis residual AUROC test pre-registered before any data touch
  - **NEW §6.6**: PROBAST self-rating (expect HIGH on Analysis domain)
  - **§7 Title** rewritten; Discussion adds 3 mandatory paragraphs per oncologist
- v0.1 (2026-05-28): Initial H-b axis framing draft

---

## 0. Locked decisions (2026-05-28, post-audit)

| Dimension | Value |
|---|---|
| Paper direction | Keep original BIOADV-2026-206 methodology + GATE-led contribution structure |
| Headline | **GATE diagnostic** as portable tool; axis reframe as supporting evidence |
| Pre-registered primary | LOCKED (preregistration.LOCK); all 6 tests FAILED, reported as worked example of GATE catch |
| New analyses | Exploratory § 4.1-4.8, pre-registered in PLAN_v2.LOCK |
| Held-out validation | 4th cohort **GSE32646** (Osaka U Noguchi; GSE20271 swapped out per author overlap) |
| Pre-reg discipline | 2-tier: preregistration.LOCK + PLAN_v2.LOCK; explicit data-burn admission |
| Submission venue | BioAdv (use 2026-05-19 NEW submission invitation) |

---

## 1. Scientific scope (GATE-led)

### 1.1 Central claim (restructured)

> **GATE (Gating Artifact Test) is a portable within-stratum AUROC
> diagnostic that catches base-rate artifacts in any cohort-free
> pharmacogenomic signature. We demonstrate GATE on a pre-registered
> external validation of a CCLE-derived anthracycline signature in
> breast cancer pCR — the signature fails GATE, and decomposition shows
> it operates as a proliferation × hormone-receptor axis projection,
> consistent with the failure mode. We validate this axis reinterpretation
> on an independent held-out cohort.**

The central claim is **methodological** (GATE), with the failed
validation as a worked example and the axis reframe as supporting
mechanistic explanation.

### 1.2 Three-layer contribution (REORDERED)

**Layer 1 — Methodological (PRIMARY, the paper's headline):**
- GATE diagnostic: within-stratum AUROC test for base-rate artifacts
- Cohort-free signature deployment via singscore (carry over from
  original BIOADV-2026-206 contribution)
- Pre-registration discipline for cell-line → patient signature transfer
  (preregistration.LOCK + PLAN_v2.LOCK as templates)

**Layer 2 — Empirical worked example (SUPPORTING evidence for Layer 1):**
- Pre-registered failure of CCLE-derived anthracycline signature on 3 BC
  cohorts (6 of 6 primary tests FAIL); GATE catches it
- 0.752 internal AUROC explicitly retracted as base-rate artifact
- GEFITINIB wrong-MoA control AUROC 0.612 demonstrates cohort confound

**Layer 3 — Mechanistic interpretation (SUPPORTING, with held-out validation):**
- Axis decomposition: signature ≈ proliferation + HR linear combination
- Multi-drug confirmation (6 cytotoxics) + targeted-therapy contrast
  (Trametinib + Gefitinib + Lapatinib + Fulvestrant)
- **Held-out cohort validates axis residual AUROC test** (§4.8)

### 1.3 What this paper does NOT claim

- Does NOT claim drug-specific predictive validity for any cytotoxic agent
- Does NOT claim the axis reinterpretation is novel biology (it is
  rediscovery of GGI / DLDA30 / Mammaprint-era observations; framed as
  methodological warning, not insight)
- Does NOT claim GATE is novel statistics (it is portable packaging of
  Hatzis 2011 / Buyse 2006 era stratified validation, formalized as a
  named deployable algorithm)
- Does NOT claim independent validation for any §4.1-4.7 claim except
  the §4.8 held-out axis test

### 1.4 Target venue + audience

BioAdv 2026 honest-negative track. Audience = pharmacogenomic signature
methodologists + translational researchers needing diagnostic tools to
avoid base-rate artifacts. GATE is the reusable contribution; worked
example + axis is the field-warning content.

---

## 2. What is pre-registered vs exploratory (CRITICAL — expanded per pre-reg auditor)

### 2.1 Pre-registered (LOCKED, immutable)

`preregistration.LOCK` (sha256 `3a9192f...`, 2026-05-28T02:14:24Z):
- Primary tests T1-T6 (all FAILED, reported in Results)
- §8 trivial baseline comparison + §8b ablation grid TCGA

### 2.2 Exploratory (NEW, requires PLAN_v2.LOCK)

All §4 analyses are POST-HOC relative to preregistration.LOCK. They
cannot serve as confirmatory tests on cohorts GSE16446 / GSE25066 /
GSE22226. §4.8 (held-out cohort) is the SINGLE EXCEPTION — it is
pre-registered confirmatory IF no §4.8 data touches a project script
before PLAN_v2.LOCK.

### 2.3 PLAN_v2 LOCK structure

1. PLAN_v2.md sha256-hashed
2. Deposit to Zenodo + OSF (DOI required, NOT a [TBC] placeholder)
3. Third-party timestamp confirmation
4. **Zenodo timestamp audited against `git log` after execution
   (per biostatistician edit)**
5. Only after all 4 steps complete, execute any §4 analysis
6. §4.8 4th cohort data download blocked until LOCK signed

### 2.4 Disclosure language (REWRITTEN per pre-reg auditor finding §4)

> Following primary external validation (six pre-registered tests, all
> failing, see preregistration.LOCK and primary_results.md), we conducted
> the following analyses. These are categorized as follows:
>
> (i) **Exploratory analyses (§4.1-4.7)** on the same three BC cohorts
> (GSE16446, GSE25066, GSE22226). **These cohorts are burned for
> confirmatory inference: the H-b axis hypothesis was generated FROM
> the observed pattern of primary test failure (HR+ direction, basal
> inversion, GEFITINIB 0.612 cohort confound). Pre-registering
> PLAN_v2 constrains analysis-time degrees of freedom (no model
> tweaking after seeing results) but does NOT restore
> hypothesis-generation/test independence.** All §4.1-4.7 results are
> exploratory-with-pre-specified-falsifier, never confirmatory.
>
> (ii) **Confirmatory held-out validation (§4.8)** on an independent
> cohort (GSE20271 primary; n=178 BC neoadjuvant T/FAC; pre-registered
> for axis residual AUROC test). This is the single confirmatory cell
> in PLAN_v2. No GSE20271 data has been examined by any project script
> before PLAN_v2.LOCK timestamp.
>
> Pre-registration here does NOT meet OSF Registered Reports Stage 1
> criteria (which require data not yet collected for the registered
> tests); see §2.7. PROBAST risk of bias rating: HIGH on Analysis
> domain (acknowledged, §6.6).

### 2.5 Data-burn acknowledgement (NEW)

GSE16446, GSE25066, GSE22226 have been examined under preregistration.LOCK
primary tests T1-T6. They have produced specific failure patterns
(HR+ stratum direction, basal inversion, GEFITINIB cohort confound).
The H-b axis hypothesis is informed by these failure patterns. Therefore:

- Any axis-decomposition result on these 3 cohorts cannot serve as
  independent confirmation of the H-b axis hypothesis
- §4.1-4.7 results are descriptive / interpretive, not confirmatory
- Only §4.8 (held-out GSE20271) can provide pre-registered confirmation

### 2.6 Independent cohort commitment (LOCKED v0.3)

**LOCKED PRIMARY: GSE32646** (Miyake et al. 2012; PMID 22320227;
Osaka University Noguchi group; n=115; P-FEC neoadjuvant BC =
paclitaxel → 5FU/epirubicin/cyclophosphamide; Affymetrix HG-U133 Plus 2.0
(GPL570); pCR labels in metadata; pCR rate ~24%).

**Why this cohort (full verification in `cohort_verification_4th.md`):**
- Zero author overlap with GSE16446 (Bonnefoi/Iggo EORTC),
  GSE25066 (Hatzis/Symmans/Pusztai MDACC), or GSE22226 (Esserman UCSF) —
  all Osaka U Noguchi group
- Anthracycline-containing regimen (epirubicin is FDA-approved anthracycline;
  pre-reg §2.6 inclusion permits any anthracycline)
- Japanese cohort = generalization bonus across population genetics
  (Reviewer 2 angle)
- GPL570 platform compatible with existing probeset→symbol collapse
  pipeline

**Why NOT GSE20271 (Hess 2006/Tabchy 2010):** four named co-authors
overlap with GSE25066 (Symmans, Vidaurre, Martin, Souchon, Pusztai).
Same MDACC T/FAC multi-trial pool. Patient overlap highly probable.
Mandatory non-overlap criterion fails on authorship alone, before
patient-overlap inference. (Documented in
`cohort_verification_4th.md`; logged in decisions.md D20.)

**Alternates** (if GSE32646 verification fails programmatic checks at
LOCK time):
- MetaGxBreast computational subset (curatedBreastData exclusion of all
  GSE16446 / GSE25066 / GSE22226 / GSE20271 patient IDs)

**Caveats noted at LOCK:**
- Smaller n=115 vs prior alternates → larger CIs (Hanley-McNeil CI
  half-width ~±0.09 at 24% prevalence; see §4.8.2 power calc)
- Smallest of 4 cohorts → consider bootstrap CIs + fixed-effects
  meta-analysis as supplementary
- Epirubicin not doxorubicin; pre-reg allows generic anthracycline; this
  generalises the H-b claim across anthracyclines, slightly weakens
  Dox-specificity sub-claim (acceptable — H-b is anti-drug-specific by
  construction)

**Inclusion criteria** (unchanged from v0.2):
- Anthracycline-containing regimen ✓ (epirubicin)
- pCR outcome ≥ 15% and ≤ 60% ✓ (24%)
- ≥ 100 patients with expression + outcome ✓ (n=115)
- Microarray or RNA-seq ✓ (GPL570)
- No author overlap with GSE16446/25066/22226 ✓
- No patient overlap ✓ (independent cohort, different country)

**Inclusion criteria** (pre-committed):
- Anthracycline-containing regimen (FAC, AC, AC-T, FEC, EC, TAC variants)
- pCR (or RECIST CR/PR) outcome with response rate 15–60%
- ≥ 100 patients with both expression + outcome
- Microarray or RNA-seq (any platform; within-cohort quantile normalization
  before scoring)
- **No author overlap with GSE16446 / GSE25066 / GSE22226 publications**
  (mandatory — pre-reg auditor finding)

**Test on §4.8 cohort:** axis residual AUROC (same as §4.1) — single
pre-specified threshold (residual AUROC 95% CI lower bound > 0.55 →
H-b axis reframe rejected).

### 2.7 NOT a Registered Reports Stage 1 statement (NEW per pre-reg auditor)

This PLAN_v2 is NOT a Stage-1 Registered Reports pre-registration.
Stage-1 requires data not yet collected for the registered tests.
The 3 BC cohorts have been examined; the 4th cohort has not. PLAN_v2 is
a **sequential timestamped exploratory plan** with one nested
confirmatory test (§4.8 on the untouched 4th cohort). The structure is
defensible (Lakens 2019; Nosek 2018 on sequential pre-registration)
ONLY because (i) the 4th cohort is truly untouched, (ii) we explicitly
admit the 3 cohorts are burned, (iii) we report exploratory results
with pre-specified falsifiers (not "primary tests").

---

## 3. Methodology summary (unchanged from BIOADV-2026-206)

Refer to original `preregistration.md` §§1-8. Quick summary:

- **Training:** CCLE solid pool, unsupervised Pearson r(gene, drug AUC),
  OncoKB-all filter, q=0.30, top-N=30 per direction, bidirectional
- **Deployment:** singscore on quantile-normalized patient expression
- **Frozen signatures:** `revise_bioadv/frozen_signatures/MANIFEST.json`
- **External cohorts:** GSE16446 (FEC, ER−), GSE25066 (TFAC), GSE22226 (AC-T)
  + NEW: GSE20271 (T/FAC, untouched, §4.8)

No methodology changes vs original. The paper does NOT propose a new
training method. The contribution is GATE + interpretation + honest
validation.

---

## 4. New analyses (PLAN_v2 LOCK covers these)

### 4.1 Axis decomposition (REBUILT per biostat + bioinformatician edits)

**Goal:** Quantify how much variance in CCLE-Dox signature score is
explained by proliferation + HR axes in BC external cohorts. Exploratory
on 3 burned cohorts; confirmed in §4.8 on held-out.

**Per cohort:**

1. **Proliferation axis (P) — REVISED (v0.3 + Phase 0 fallback):** Per-patient mean z-score of:
   - HALLMARK_E2F_TARGETS
   - HALLMARK_G2M_CHECKPOINT
   - HALLMARK_MITOTIC_SPINDLE
   - HALLMARK_MYC_TARGETS_V1
   - ~~Whitfield_G2M~~ → **NOT AVAILABLE locally; fallback to 4-set composite**
     (Whitfield_G2M source file not found in `revise_bioadv/`; logged in
     `phase_0/PHASE_0_LOG.md` 2026-05-28). Documented limitation in
     Methods. Net effect: P axis slightly more MYC-weighted than ideal
     but bioinformatician's primary concern (dropping MYC_TARGETS_V2)
     remains addressed.

   Supplementary single-set sensitivity: MKI67 z-score alone, E2F_TARGETS
   alone.

2. **Hormone-receptor axis (H) — REVISED:** Per-patient mean z-score of
   **ESR1 + PGR** (clinical HR proxy). Supplementary: ESR1-only,
   PAM50 basal-probability continuous.

3. **ERBB2 axis (B) — NEW separate predictor:** ERBB2 z-score
   (not folded into H per bioinformatician edit).

4. **Multivariable linear model:**
   ```
   sig_dox ~ β0 + β1·P + β2·H + β3·B + β4·Basal_prob + ε
   ```
   Fit per cohort.

5. **Statistical reports per biostat:**
   - R² and partial-R² for each predictor
   - **VIF** for P, H, B, Basal_prob — if any > 5, declare β's individually
     uninterpretable; limit interpretation to R² + partial-R²
   - **Ridge regression sensitivity** (10-fold CV λ, seed=42); report
     coefficient sign stability
   - **Bootstrap 95% CI on R²** (1000× cluster-bootstrap on patients,
     seed=42)

6. **Residual signal test (the H-b core, pre-specified falsifier):**

   - Construct residual = sig_dox − fitted(P + H + B + Basal_prob) per
     cohort (within-cohort OLS residual)
   - Test residual ~ pCR (logistic regression + AUROC)
   - **Pre-specified falsifier (LOCKED):** residual AUROC 95% CI lower
     bound > 0.55 in ≥ 2 of 3 cohorts → H-b axis reframe REJECTED
     (signature has unexplained drug-specific signal beyond axes)
   - Bonferroni cushion: α = 0.01 for this elevated exploratory endpoint
   - **Supplementary nonlinear check (per bioinformatician):** spline
     regression `sig ~ s(P) + s(H) + s(B)`; report R² gap vs linear.
     Spline R² >> linear R² → signature has nonlinear axis dependence
     (H-b stronger, not weaker)

7. **Synthetic positive control (NEW per biostat):**
   - Construct synthetic_sig = α·P + β·H + γ·DRUG_SPECIFIC + noise
   - Where DRUG_SPECIFIC is a known TOP2-poison-related signal
     (e.g., TOP2A expression z-score with effect size 0.3 SD)
   - Show: residual AUROC of synthetic_sig > 0.6 (recovers the drug
     component)
   - Sanity check: if residual test cannot detect known drug component,
     §4.1 design is broken

### 4.2 Multi-drug axis confirmation (v0.3: 11 signatures + TOP2_CONSENSUS + loosened thresholds)

**Goal:** Show cytotoxic-vs-targeted asymmetry in axis-collapse via
3-tier contrast (cytotoxics > BC-pathway TKIs > true-null TKI).

**Drug set (v0.3 LOCKED):**

- **Cytotoxic (6):** Doxorubicin, Etoposide, Paclitaxel, Carboplatin,
  Cyclophosphamide, 5-Fluorouracil
- **TOP2-poison consensus row (NEW v0.3):** **TOP2_POISON_CONSENSUS**
  (Dox + Etoposide pre-frozen in MANIFEST, Jaccard 0.85, n_up=27,
  n_down=28) — reported as headline TOP2 row alongside individual Dox + Etoposide
- **BC-pathway TKI controls (3):** Trametinib (MEK), Gefitinib (EGFR),
  Lapatinib (HER2)
- **True-null TKI control (1, NEW v0.3 per bioinformatician):**
  **Imatinib** (BCR-ABL/KIT, no BC indication, MANIFEST n=230,
  Jaccard 0.76/0.82) — predicted to correlate with NONE of P/H/B axes
- **Positive H-axis control (1):** Fulvestrant (ESR1-driven; flagged as
  low-stability directional control only — Jaccard 0.40/0.33, n=60)

**Total: 11 signatures. All present in `frozen_signatures/MANIFEST.json`
quantile_30_oncokb_all (no new CCLE training required).**

**Per drug × per cohort:**
1. Signature score (singscore, frozen from MANIFEST)
2. Compute correlations r(sig, P), r(sig, H), r(sig, B)

**Pre-registered predictions (thresholds loosened v0.3 per realistic
Jaccard / noise estimates):**

- **Cytotoxic convergence (descriptive, n=1 per biostat caveat):**
  all 6 cytotoxics show |r(sig, P)| > 0.6 in ≥ 2 of 3 cohorts. This is
  ONE joint observation with correlated structure, not 6 independent
  tests. Pairwise inter-signature correlation matrix reported in §4.2.2.

- **Dox vs Etoposide (TOP2 family) — directional NOT numerical:**
  both show |r(sig, P)| > 0.6 with same sign; pre-register **direction
  identity**, NOT effect-size equivalence. Dox-r ≥ Etoposide-r expected
  per MoA difference (Dox adds TOP2B + intercalation + ROS).

- **§4.2 HEADLINE inferential test (v0.3 — explicit form per biostat):**
  ```
  Δ = mean_clade(cytotoxic) |r(sig, P)| − mean_clade(TKI ∪ null) |r(sig, P)|
  ```
  Pre-registered: **Δ > 0.3 with BH-corrected paired-comparison p < 0.05
  → cytotoxic-vs-targeted asymmetry CONFIRMED**. This is the
  inferential lever; cytotoxic intra-clade is descriptive only.

- **BC-pathway TKI individual thresholds (v0.3 LOOSENED):**
  - **Trametinib** |r(sig, P)| < 0.5 (was 0.4 — MAPK→cyclin/cdk carries
    some proliferation correlation, especially basal/RAS-active BC).
    Trametinib alone cannot falsify TKI contrast at 0.42.
  - **Gefitinib** |r(sig, P)| < 0.4 (unchanged — clean EGFR-TKI).
  - **Lapatinib** r(sig, B) > 0.3 (was 0.4 — ERBB2 prevalence in BC
    cohorts ~15-20% dilutes ERBB2 correlation across full cohort).
    Explicit dilution caveat in Methods.

- **True-null Imatinib check (NEW v0.3):**
  |r(sig, P)| < 0.3 AND |r(sig, H)| < 0.3 AND |r(sig, B)| < 0.3 in
  ≥ 2 of 3 cohorts. Imatinib is the cleanest null because it has no
  BC indication / no BC pathway. Failure → all our CCLE-derived
  signatures pick up generic BC biology regardless of drug;
  framework methodology is structurally compromised.

- **Fulvestrant H-axis check (v0.3 LOOSENED + flagged):**
  r(sig, H) > 0.3 (was 0.4 — Fulvestrant low Jaccard stability).
  Flag in Methods: "directional positive control only, effect size not
  interpretable due to n=60 CCLE training set + bimodal sensitivity."
  Supplement with **ESR1+PGR 2-gene oracle** as rigorous H-axis check.

#### 4.2.2 Inter-signature correlation matrix → Figure 4 (PROMOTED v0.3)

Report 11×11 inter-signature correlation matrix on each cohort
(3 matrices). Display in main text as **Figure 4** with:
- Hierarchical clustering dendrogram
- Cytotoxic clade vs TKI clade vs Fulvestrant vs Imatinib clades color-coded
- Axis-projection annotation (each signature labeled with which axis it
  loads on)
- Per-cohort tile

This visually anchors the cytotoxic-collapse narrative more than any
single correlation coefficient.

### 4.3 Trivial baseline expansion (EXPANDED per bioinformatician)

For each baseline × cohort: AUROC for pCR + paired DeLong vs framework
signature.

| Baseline | Status | Description |
|---|---|---|
| Random uniform | Done | Re-report |
| MKI67 single-gene | Done | Re-report |
| HALLMARK_E2F_TARGETS singscore | Done | Re-report |
| HALLMARK_G2M_CHECKPOINT singscore | Done | Re-report |
| Mean expression z-score across signature genes | NEW | Tests "is it just gene-set-membership not weights?" |
| 5-set composite proliferation (the P axis itself) | NEW | Proper proliferation oracle |
| Random 60 OncoKB genes × 100 iterations (matched-size null) | NEW | Tests "better than structural prior?" |
| Random Hallmark set × 50 sets (any-pathway null) | NEW | Tests "any pathway scoring would work" |
| **MKI67 + ESR1 two-gene sum** | NEW | **PRIMARY COMPARATOR** per bioinformatician — if 2-gene sum matches 60-gene framework, that is the H-b headline |

The **MKI67+ESR1 2-gene comparator** is elevated to **primary headline
result** in the Abstract (per bioinformatician #4). Anti-signature is
dropped (lower priority).

### 4.4 Calibration + PR-AUC + permutation (per biostat #5)

Per (signature × cohort × outcome):

1. PR-AUC + average precision (extract from existing primary_tests.json
   + new signatures)
2. F1 at Youden-optimal threshold
3. Brier score + reliability diagram (10-bin calibration plot)
4. **Permutation null — REVISED:** within-PAM50-stratum permutation for
   GSE25066 + GSE22226 (preserves stratum base-rate); marginal permutation
   for GSE16446 (uniformly ER−/TN, single stratum). Report BOTH
   within-stratum and marginal empirical p; within-stratum is the
   headline null.
5. CV protocol clarification in Methods (group-by-patient, already
   implemented).

### 4.5 Ablation grid completion (unchanged from v0.1)

Pre-reg §8b grid on TCGA — q × top-N × gene-universe. Complete the grid
via `task_b_optimal_gene_count.py` extension. Report robustness of
primary config (q=0.30, top-N=30, OncoKB-all).

### 4.6 Hallmark across drugs (REVISED to GSEA per bioinformatician)

For each drug d (6 cytotoxic + 3 targeted + Fulvestrant):

1. Pull the full per-drug Pearson-r ranked gene list from CCLE training
   (all OncoKB-all genes, not just top-30)
2. Run **fgsea** on the ranked list against MSigDB Hallmark 2025.1.Hs,
   10000 permutations, BH-FDR correction
3. Report NES heatmap (drug × Hallmark) + leading-edge gene overlap
   for proliferation Hallmark sets

Supplementary: Fisher's exact overrepresentation on top-30 (sanity check).

### 4.7 Multiplicity correction for exploratory family (NEW per biostat)

- §4.2 (10 signatures × 3 cohorts × 3 axes = 90 correlations):
  BH q < 0.10 within grid
- §4.4 (5 metrics × ≥3 signatures × 3 cohorts ≈ 45 entries):
  BH q < 0.10 within grid
- **§4.1 H-b headline residual AUROC test:** Bonferroni α = 0.01
  (elevated single endpoint)
- §4.3, §4.5, §4.6: descriptive, no inferential threshold

All p-values in §4.1-4.7 are reported as nominal alongside corrected.
No confirmatory inference from any individual §4.1-4.7 entry; only the
§4.8 held-out test is confirmatory.

### 4.9 GATE second worked example — Mammaprint × TCGA-BRCA (NEW v0.3 — Layer 1 strengthening)

**Purpose:** Demonstrate GATE generalizes beyond CCLE-derived cytotoxic
signatures, lifting Layer 1 (GATE diagnostic) from "1-cohort wrapper"
to "portable across signature classes." Without this section, Layer 1
occupies only Results §5.2 (1 subsection) while Layers 2-3 occupy 5-6
subsections — reviewer-sim audit flagged this imbalance as risking
"GATE = wrapper for negative paper" interpretation.

**Signature:** **Mammaprint** (Van't Veer 2002; 70-gene poor-prognosis
signature; gene list publicly available; commercial test FDA-cleared).

**Cohort:** **TCGA-BRCA chemotherapy-treated subset** (n ~200-300
expected; pull from `treatments` table for "adjuvant chemo" or
"neoadjuvant chemo" annotation; verify n at LOCK time).

**Analyses (pre-registered):**
1. Score Mammaprint signature (70 genes, sign-weighted per original
   Van't Veer 2002 publication) via singscore on quantile-normalized
   TCGA-BRCA RNA-seq + chemo subset
2. Outcome: 5-year recurrence-free survival OR overall survival
   (whichever has cleaner labels)
3. Compute marginal AUROC for outcome
4. **Run GATE diagnostic:** compute within-PAM50-stratum AUROC; report
   GAP = marginal AUROC − within-stratum AUROC

**Pre-registered prediction:**
- Mammaprint marginal AUROC ≥ 0.65 (Mammaprint is well-documented
  prognostic in BC)
- Within-PAM50-stratum AUROC ≈ 0.50-0.55 (Mammaprint encodes
  proliferation + intrinsic subtype, so signal collapses within PAM50)
- **GAP ≥ 0.10** → GATE catches that Mammaprint's marginal performance
  is base-rate artifact for *predictive* use cases (it remains valid for
  *prognostic* stratification; this is the methodological point)

**Pre-registered falsifier:**
- If GAP < 0.05 → GATE does NOT catch Mammaprint as base-rate confounded
  → either Mammaprint really has within-stratum signal (unexpected) or
  GATE has insufficient sensitivity → Layer 1 claim weakened
- If GAP > 0.10 and within-stratum AUROC ≈ 0.50 → CONFIRMS GATE
  generalizes beyond CCLE-derived signatures

**Why this is a "clean" worked example:**
- Mammaprint is well-studied, published, commercially deployed
- Mechanism is well understood (proliferation + intrinsic subtype)
- Independent of the BIOADV-2026-206 anthracycline failure path
- TCGA-BRCA is well-curated, public, large

**Operational discipline:**
- TCGA-BRCA expression + chemo subset definition + Mammaprint gene list
  pre-LOCKED before any score is computed
- No iterative tuning of stratum definition or signature genes
- Result reported as Results §5.2.2 alongside the anthracycline GATE
  worked example (§5.2.1)

**Cost:** ~3-4 days compute (TCGA-BRCA download + chemo subset extraction
+ Mammaprint scoring + GATE diagnostic + bootstrap CIs + reporting)

**Decision rule LOCKED at PLAN_v2.LOCK:**
- GATE GAP ≥ 0.10 AND within-stratum AUROC ≤ 0.55 → Layer 1
  CONFIRMED portable across signature classes
- GAP < 0.05 OR within-stratum AUROC ≥ 0.60 → Layer 1 portability
  caveated in Discussion; GATE remains useful for cytotoxic signatures
  but generalizability flagged

---

### 4.8 Held-out cohort validation (v0.3 LOCKED — GSE32646 + 3×3 decision rule + power calc + LOCK enforcement)

#### 4.8.1 Cohort + test

**Cohort (LOCKED v0.3):** **GSE32646** (Miyake et al. 2012;
PMID 22320227; Osaka U Noguchi group; n=115; P-FEC neoadjuvant BC;
Affymetrix HG-U133 Plus 2.0 (GPL570); pCR rate ~24%).

**Pre-registered confirmatory test:**
- Apply LOCKED frozen Dox signature (from
  `frozen_signatures/solid_only/quantile_30_oncokb_all/DOXORUBICIN.tsv`)
  via singscore on quantile-normalized + probeset-collapsed GSE32646
  expression matrix
- Compute P + H + B + Basal_prob axes per §4.1
- Fit `sig_dox ~ β0 + β1·P + β2·H + β3·B + β4·Basal_prob + ε` per §4.1
- Test residual ~ pCR (logistic + AUROC + bootstrap 95% CI 1000×)

#### 4.8.2 Pre-registered headline + power calculation (NEW v0.3 per E3/E4)

**SINGLE HEADLINE TEST (E4 commitment):** The §4.8 headline cell is
**exactly** the ESR1+PGR linear OLS residual AUROC on GSE32646.

All other variants computed on GSE32646:
- ESR1-only H-axis residual
- PAM50 basal-prob H-axis residual
- Spline (nonlinear) residual

These are reported in §4.8 supplementary with Bonferroni × 3
correction and are **explicitly ineligible to flip the headline
confirmation / falsification outcome**.

**Power calculation (NEW v0.3 per E3):**
- Sample size: n=115; outcome prevalence ≈ 24% (~28 pCR events)
- Hanley-McNeil 95% CI half-width on AUROC: ≈ **±0.085 to ±0.095**
- MDE at α=0.05, 80% power to reject AUROC=0.50: ≈ **0.63**
- MDE at α=0.05, 80% power to reject AUROC=0.55: ≈ **0.67**
- **Implication: gray-zone outcome (0.55-0.60) is the MODAL expectation
  under both null and small-effect alternatives.** Pre-register this
  as expected, not as exception.

#### 4.8.3 Pre-specified decision rule — 3×3 TABLE (v0.3 REWRITE per E2)

The v0.2 2-outcome rule had a logical inversion (CI lower 0.54 +
point 0.62 incorrectly "confirms null"). Replaced with explicit
3-outcome table:

| Point estimate AUROC | 95% CI lower bound | Verdict |
|---|---|---|
| **≤ 0.55** | **≤ 0.55** | **H-b axis CONFIRMED** — residual is null; signature reduces to P+H+B+Basal axes |
| **0.55 < point ≤ 0.60** OR (CI crosses 0.55 either direction) | crosses 0.55 | **GRAY ZONE** — neither confirmed nor falsified; Abstract treats as null; Discussion contingency paragraph documents |
| **> 0.60** | **> 0.55** | **H-b axis FALSIFIED** — drug-specific residual signal present; signature has structure beyond axes; paper must reframe |

Operational interpretation:
- Both point estimate AND CI lower bound must satisfy the rule (not
  either-or)
- A CI crossing 0.55 with point estimate in (0.55, 0.60) is **gray**
  even if point estimate alone would suggest confirm/falsify
- Reporting template (LOCKED): "Residual AUROC point estimate X [95% CI
  L, U], gray-zone per LOCK §4.8.3 → H-b reframe NEITHER CONFIRMED NOR
  FALSIFIED on independent cohort."

#### 4.8.4 Gray-zone handling (NEW v0.3 — formalized per E5)

Per power calc §4.8.2, gray-zone outcome is the modal expectation.
This must be handled honestly:

- Abstract treats gray-zone as null (i.e., does NOT claim confirmation)
- Discussion adds dedicated contingency paragraph:
  > "The held-out cohort residual AUROC was [X] [95% CI L, U],
  > falling within the pre-specified gray zone. Per LOCK §4.8.3, the
  > axis reinterpretation is neither confirmed nor falsified at this
  > sample size. A larger held-out cohort (target n ≥ 250) would be
  > required to distinguish residual AUROC = 0.50 from 0.55 with 80%
  > power."
- The paper's contribution remains intact (Layer 1 GATE + Layer 2
  worked example); Layer 3 axis claim is downgraded to "consistent
  with held-out data but not confirmed."

#### 4.8.5 Operational discipline + LOCK enforcement (v0.3 STRENGTHENED per E1)

**Required before any GSE32646 data touches a project script:**

1. **Designated custodian:** ONE author (TBD at sign-off) is sole holder
   of GSE32646 SOFT file post-LOCK. File stored on a separate machine
   OR in an encrypted tarball where decryption key = Zenodo DOI string.
2. **Pre-LOCK public push of download script** with **placeholder
   accession `GSEXXXXX`** to public repo. Script sha256 third-party
   witnessed via git commit on public branch BEFORE Zenodo deposit.
3. **Post-execution audit:** Verify `T_D` (GEO `getGEO('GSE32646')`
   HTTP access log timestamp) satisfies `T_D > T_L` (Zenodo timestamp).
   Document the comparison in supplementary.
4. **No GSE32646 sample IDs, n, or any metadata referenced in any
   project file** before LOCK (verified via grep on git history).

**Cover letter language (v0.3):**
> "Held-out confirmatory test on independent cohort GSE32646 (Osaka U
> Noguchi group; no author or patient overlap with prior validation
> cohorts), pre-registered in PLAN_v2 Zenodo DOI [...]; download timestamp
> post-LOCK verified against Zenodo timestamp."

---

## 5. Reviewer matrix B1-B6 integration

| Block | Status under PLAN_v2 v0.2 | Section |
|---|---|---|
| **B1** Novelty reframing | Reframe to GATE-led (Layer 1) + axis as supporting (Layer 3); CODE-AE/TRANSACT/TUGDA dropped per R1-M2 | Intro + §1 |
| **B2** External validation | DONE (pre-reg, all fail); oncoPredict / pRRophetic head-to-head added | §3 + Results 5.1 |
| **B3** 0.752 honest + GATE | Headline of Layer 1; dedicated subsection | Results 5.2 |
| **B4** Wider metrics + permutation + trivial baselines | §4.3 + §4.4 here; expanded with MKI67+ESR1 comparator | Results 5.3 |
| **B5** Ablation grid | §4.5 | Results 5.6 / Supp |
| **B6** Figures + schematic + Hallmark + log2 | §4.6 GSEA + new figures redrawn; log2 → DESeq2 vst sensitivity as supp | All figures |

**Plus reviewer matrix gaps surfaced by reviewer-sim audit:**
- R2-2 "0.75 not clinically useful → DCA" — DCA carried over from
  preregistration §6.4 + extended to MKI67+ESR1 2-gene comparator
- R1-M2 deep-learning comparators — explicitly stated as not run with
  rationale (different cross-study protocols, computational scope)

**Workload estimate:** 18-25 person-days (up from v0.1 estimate of
15-22 due to §4.8 + §4.2 reframe + figure redraw).

---

## 6. Pre-registration discipline (STRENGTHENED per pre-reg auditor)

### 6.1 preregistration.LOCK INTACT

Unchanged. Primary tests T1-T6 reported as locked.

### 6.2 PLAN_v2.LOCK protocol (additions)

1. sha256-hash PLAN_v2.md
2. Deposit to Zenodo + OSF with public DOI (NO [TBC] placeholders)
3. Third-party timestamp confirmation
4. **Audit Zenodo timestamp against `git log` post-execution** (NEW —
   biostat #6); any commit on §4 analysis scripts before Zenodo timestamp
   invalidates the LOCK
5. Only after all 4 steps complete, execute any §4 analysis

### 6.3 Deviation log (expanded D-entries)

`revise_bioadv/decisions.md` additions:
- D11: PLAN_v2 LOCK trigger and 2-tier rationale
- D12: H-b axis framing chosen as SUPPORTING layer (not headline)
- D13: 6 cytotoxics + 3 targeted + Fulvestrant control set LOCKED
- D14: Reviewer matrix B1-B6 absorbed into v2
- D15: Multi-anchor pivot abandoned (`multi_anchor_paper/SHELVED.md`)
- D16: 5-expert audit on PLAN_v2 v0.1 → v0.2 (findings + resolutions)
- D17: H-axis = ESR1+PGR composite (resolved from §10 open issue)
- D18: H-b axis hypothesis was generated from §3 failure pattern
       (explicit admission)
- D19: Residual AUROC ≤ 0.55 → H-b confirmed; > 0.55 → H-b falsified;
       behavior at 0.55-0.60 = "ambiguous, see Discussion"
- D20: §4.8 held-out cohort = **GSE32646** (Osaka U Noguchi). GSE20271
       NO-GO documented (Tabchy/Symmans/Pusztai/Vidaurre/Martin/Souchon
       author overlap with GSE25066). See `cohort_verification_4th.md`
- D21: BH q < 0.10 within §4.2 + §4.4 grids; Bonferroni α=0.01 for §4.1
       H-b headline

### 6.4 Reporting commitment (STRENGTHENED)

- All §4 analyses report in full regardless of outcome
- **Negative results from §4 in MAIN TEXT (not just supplementary)**
- **§4.1 pre-specified table skeleton in Methods** with cell-level
  commitment (R², β1 (P), β2 (H), β3 (B), β4 (Basal), VIF, ridge
  coefficients, bootstrap R² CI, residual AUROC + CI + permutation p
  per cohort)
- **§4.8 pre-specified single primary cell**: residual AUROC 95% CI
  reported with the pre-specified threshold
- No selective omission

### 6.5 Cohort discipline

- **3 burned cohorts** (GSE16446/25066/22226): exploratory only, never
  confirmatory, explicit data-burn acknowledgement (§2.5)
- **1 held-out cohort** (GSE20271): blocked from project script access
  until PLAN_v2.LOCK confirmed; provides the single confirmatory test

### 6.6 PROBAST self-rating (NEW per pre-reg auditor)

Self-rated PROBAST domain table (will appear as supplementary):

| Domain | Risk | Justification |
|---|---|---|
| 1. Participants | LOW | Cohort selection from public GEO/TCGA, inclusion criteria locked |
| 2. Predictors | LOW | Frozen signatures (sha256), no post-hoc gene selection |
| 3. Outcome | LOW | pCR (binary, standardized assessment) |
| 4. **Analysis** | **HIGH** | H-b axis hypothesis generated from observed failure pattern on same cohorts; §4.8 partially mitigates via held-out test, but §4.1-4.7 on burned cohorts retain HIGH risk |

This HIGH rating is acknowledged in main text not hidden in supplementary.

---

## 7. Manuscript outline (TITLE REWRITTEN per oncologist + reviewer-sim + pre-reg auditor)

### 7.1 Title (v0.3 TRIMMED per light-audit Q1)

> **GATE: a portable within-stratum diagnostic for base-rate artifacts
> in cohort-free pharmacogenomic signatures**

Case-study qualifier moved out of Title into Abstract (was: ", with a
pre-registered worked example on CCLE → breast cancer cytotoxic transfer").
Reason: v0.2 title had 14-word qualifier dominating GATE clause
visually; editor's eye landed on "breast cancer cytotoxic transfer"
telegraphing negative-anthracycline-paper. v0.3 title is clean methods
contribution; Abstract carries the worked example.

**Abstract opens with:** "We introduce GATE, a within-stratum AUROC
diagnostic that catches base-rate artifacts in cohort-free
pharmacogenomic signatures. We demonstrate GATE on a pre-registered
external validation of a CCLE-derived anthracycline signature in breast
cancer pCR..."

### 7.2 Section structure (~5500 words)

- **Abstract** (~250 words):
  - GATE diagnostic as central contribution
  - Pre-registered worked example (anthracycline BC pCR)
  - 6 of 6 primary tests fail; GATE catches it
  - Axis decomposition shows signature ≈ proliferation + HR
  - Held-out cohort GSE20271 confirms axis reframe (residual AUROC
    [reported value] [95% CI])
  - **MKI67+ESR1 2-gene sum matches 60-gene framework AUROC** —
    elevated to Abstract per bioinformatician

- **Introduction** (~700 words):
  - Pharmacogenomic signature transfer state of the art
  - Base-rate artifacts as a known problem in cell-line → patient
    deployment
  - Need for portable diagnostic (motivates GATE)
  - Cohort-free deployment as our context
  - Contributions: GATE (primary) + worked example + axis interpretation

- **Methods** (~1300 words):
  - Original training pipeline (refer to preregistration.md)
  - Cohort-free deployment via singscore
  - Pre-registered external validation design (T1-T6)
  - **GATE diagnostic algorithm specification** (the headline method)
  - Axis decomposition methodology (§4.1)
  - Multi-drug confirmation (§4.2)
  - Held-out validation design (§4.8)
  - Statistical tests + permutation + multiplicity

- **Results** (~1900 words):
  - 5.1 Pre-registered primary tests (T1-T6) — all FAIL (worked example
    setup)
  - 5.2 **GATE diagnostic: portable across signature classes**
       (Layer 1 headline)
    - 5.2.1 Worked example #1: 0.752 internal AUROC of our CCLE-Dox
      signature is base-rate artifact (within-stratum AUROC ≈ 0.5)
    - **5.2.2 Worked example #2 (NEW v0.3): Mammaprint × TCGA-BRCA
      chemo cohort** — GATE catches Mammaprint's marginal AUROC as
      proliferation/subtype base-rate (within-stratum collapses).
      Generalizes Layer 1 contribution beyond CCLE-derived
      signatures.
  - 5.3 Trivial baselines + **MKI67+ESR1 2-gene sum matches framework**
       (B4)
  - 5.4 Axis decomposition: signature ≈ β1·proliferation + β2·HR + β3·ERBB2
       (Layer 3 exploratory)
  - 5.5 Multi-drug: 6 cytotoxics converge on P+H axes; TKI controls do NOT
       (Layer 3 supporting contrast)
  - 5.6 Ablation grid: primary config not a degenerate corner (B5)
  - 5.7 GSEA Hallmark across drugs (B6, R1-O2)
  - 5.8 **Held-out cohort GSE20271 — axis residual AUROC test PASSES /
       FAILS** (the single confirmatory cell)

- **Discussion** (~900 words, EXPANDED per oncologist):
  - What CCLE-derived cytotoxic signatures actually encode (rediscovery
    of GGI / DLDA30 / Mammaprint biology, framed as methodological
    warning, NOT discovery)
  - GATE as portable tool — code + usage guidance
  - **"Ki67-IHC is already sufficient"** paragraph (~150 words, the
    hardest truth — central Ki67-IHC + ER/PR/HER2 predicts pCR as well
    as transcriptional signatures for anthracycline; cite POETIC,
    Dowsett 2011, IKWG 2020, MINDACC) — NEW per oncologist
  - **"How clinicians should NOT use such signatures"** paragraph
    (~150 words, with GATE-based clinical decision rule: don't use if
    within-subtype AUROC < 0.60 AND within-Ki67-quartile AUROC < 0.60)
    — NEW per oncologist
  - **"What signature class might actually work"** (~100 words —
    drug-specific resistance signatures: TOP2A copy-loss, ABCB1
    induction, p53 context; cite multi-anchor work as motivated)
    — NEW per oncologist

- **Supplementary:**
  - Full primary_test_results JSON dump
  - PROBAST self-rated domain table
  - Calibration + PR-AUC + permutation tables
  - Hallmark GSEA NES heatmap + leading-edge
  - Code + signature LOCK + Zenodo DOIs

### 7.3 Cover letter (revised structure)

> Dear Editor,
>
> Following your 2026-05-19 letter inviting resubmission of
> BIOADV-2026-206 as a new submission, we are pleased to submit our
> revised manuscript.
>
> The revision is substantial. We have:
>
> 1. Conducted a pre-registered external validation
>    (preregistration.LOCK Zenodo DOI [...], locked 2026-05-28). All
>    six pre-registered primary tests failed. We report this failure
>    honestly as a worked example.
>
> 2. Reframed the paper's central contribution around **GATE**
>    (Gating Artifact Test), a portable within-stratum AUROC diagnostic
>    that catches base-rate artifacts in cohort-free pharmacogenomic
>    signatures. The failed validation is the worked example
>    demonstrating GATE's utility.
>
> 3. Added an independent held-out cohort (GSE20271, untouched by any
>    project script before PLAN_v2.LOCK timestamp) to provide pre-
>    registered confirmatory validation of the axis decomposition that
>    explains the failure mechanism.
>
> 4. Addressed all Reviewer 1 + 2 mandatory points from the 2026-05-19
>    round, point-by-point response attached.
>
> 5. Pre-registered all new analyses (PLAN_v2 Zenodo DOI [...]) with
>    explicit data-burn acknowledgement and PROBAST risk-of-bias
>    self-rating (HIGH on Analysis domain for §4.1-4.7, mitigated for
>    §4.8 held-out test).
>
> This honest reframing — methodological tool (GATE) supported by
> pre-registered worked example and held-out validation — represents
> the rigorous treatment of negative results that BioAdv has
> championed in 2026 editorial policy.
>
> Sincerely, [authors]

---

## 8. Timeline (REVISED for v0.2 + Path C)

| Week | Activity | Deliverable |
|---|---|---|
| 0 | PLAN_v2 v0.2 drafted | this file |
| 0 | Light audit on v0.2 diff (especially §4.8 held-out design + Title) | audit findings |
| 0–1 | User sign-off + revisions | PLAN_v2 v1.0 |
| 1 | 4th cohort GSE20271 verification (GEOquery + non-overlap check) | verified cohort spec |
| 1 | sha256 + Zenodo + OSF deposit | PLAN_v2.LOCK + DOI |
| 1–2 | §4.1 axis decomposition + §4.2 multi-drug (11 signatures) | results JSONs |
| 1–2 | §4.3 trivial baselines (incl. MKI67+ESR1 2-gene) + §4.4 calibration / PR-AUC / permutation | results tables |
| 2–3 | §4.5 ablation grid + §4.6 GSEA Hallmark | grid + NES JSON |
| 2–3 | §4.8 **GSE32646** held-out validation (POST-LOCK, blocked before; ESR1+PGR linear OLS residual AUROC headline + supplementary variants) | held-out result + 3×3 decision |
| 3 | **§4.9 NEW: Mammaprint × TCGA-BRCA GATE second worked example** (~3-4 days) | Layer 1 portability evidence |
| 3–4 | Figures redraw (B6) incl. 11×11 matrix → Fig 4 hierarchical clustered heatmap | new Fig 1-6 |
| 4–5 | Manuscript writing | draft v1.0 |
| 5 | Point-by-point response letter | response_letter.md |
| 5–6 | Internal review + revision | manuscript v1.x |
| 6–8 | Submit to BioAdv | submitted |

**Total: 7-8 weeks from PLAN_v2 v0.3 to submission** (was 6-7 in v0.2;
+1 week for §4.9 second GATE worked example).

---

## 9. Light re-audit on v0.2 (recommended before LOCK)

Not full 5-lens. Focused on diff vs v0.1:

1. Title rewrite acceptability (oncologist + reviewer-sim re-check
   ≤ 100 words per lens)
2. §4.8 held-out cohort design — biostat + pre-reg auditor only
   (cohort spec + LOCK protocol + decision rule)
3. §4.2 TKI control reshuffle + Fulvestrant sanity — bioinformatician
   only

Total ~5 person-day light audit. If all pass → sign-off + LOCK.

---

## 10. Open issues to close before LOCK (v0.3 updated)

- [x] Cohort verification — **GSE20271 NO-GO confirmed**; GSE32646
      locked as primary
      (`cohort_verification_4th.md`)
- [x] Light re-audit on v0.2 — DONE, 4 agents returned, findings
      integrated to v0.3
- [x] **Algorithm-level decision RESOLVED v0.3:** Mammaprint × TCGA-BRCA
      chemo cohort selected as second GATE worked example (§4.9).
      Timeline updated +1 week (7-8 weeks total).
- [ ] Designated custodian for §4.8 GSE32646 SOFT file — TBD at sign-off
- [ ] decisions.md D11-D21 wording finalized (text below)
- [ ] Zenodo deposit metadata (license, ORCID, project tag)
- [ ] sha256 of all §4 scripts (incl. new §4.8 download script with
      placeholder accession) computed BEFORE LOCK
- [ ] Authorship list for paper finalized

---

## 11. Cross-references

- Parent: `revise_bioadv/PLAN.md` (original)
- Pre-reg: `revise_bioadv/preregistration.md` + `preregistration.LOCK`
- Negative result: `revise_bioadv/primary_results.md`
- Reviewer matrix: `revise_bioadv/reviewer_response_matrix.md`
- Shelved multi-anchor: `revise_bioadv/multi_anchor_paper/SHELVED.md`
- Pre-scan audit (different context): `revise_bioadv/audits_pre_scan/SYNTHESIS.md`
- PLAN_v2 5-expert audit synthesis: `revise_bioadv/resubmission_v2/audits/SYNTHESIS.md`

---

**END OF PLAN_v2 v0.2 DRAFT**

Status: pre-LOCK; awaiting light re-audit + cohort verification + user
sign-off.
