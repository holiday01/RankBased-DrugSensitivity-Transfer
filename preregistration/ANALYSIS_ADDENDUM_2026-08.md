# Analysis addendum

For the revision of the study deposited at **10.5281/zenodo.20726403**.

**Parent pre-registration**: `preregistration.md`, sha256
`3a9192f3adf0da350d34cf978233799ba571d7e74a5392a2b56201be20f485f7`,
locked `2026-05-28T02:14:24+00:00` (companion `preregistration.LOCK`).

**Status**: DRAFT — not yet locked. Becomes binding only when §12 is completed and the
document is deposited.

**Written**: 2026-07-22, after receipt of the major-revision decision on 2026-07-22 and
**before any revision analysis (A2–A13) was run.** The corrections in §3 were established
first, because several of them determine what the new analyses must be compared against.

**Reference convention.** A bare §N refers to a section of *this* document. Sections of any
other document are always named: "the parent's §6.3", "PLAN_v2 §2.2", "the manuscript's
§2.6". Line numbers in `main.tex` refer to the 943-line 2026-06-28 build, sha256
`1f69058b…`, and every cited line is also quoted verbatim so the reference survives
renumbering.

---

## 1. What this document is, and what it is not

It **is** a pre-specification of the analyses added in response to the reviewers, filed
before those analyses are run, so that their estimators, decision rules and reporting
commitments cannot be selected after the results are seen.

It **is not** a claim that these analyses were pre-registered in the original sense. They
were not. They are new work, specified in advance of execution, and every result arising
from them will be reported as **pre-specified for the revision**, never as *confirmatory*
in the sense of the parent pre-registration. The parent's confirmatory family (T1–T6) is
closed; it was executed and all six tests failed.

It also **corrects the record** (§3). Every one of its subsections records a correction found
by the authors during revision planning, not raised by any reviewer.

**On results that already exist in this document.** **§3, §4.1, §5.2, §6 and §10** — that
enumeration and no other — quote figures the authors had already computed before this document
was drafted: the corrected ΔAUROC table, the scorer-fix consequences and R² ranges, the
de-duplication balance figures, the smallest observed raw LRT p, and the count of failing
minimal subsets. They are stated because withholding a disclosed unfavourable result would be
worse than the awkwardness of disclosing it before the lock.

**All of them are EXPLORATORY and none is re-reported anywhere as pre-specified.** No
pCR-dependent figure outside that enumeration is added to this document before the lock, and
every pCR-dependent figure computed after the lock is governed by §6, §8 and §9. What made
the pre-lock work admissible is that its *load-bearing* findings were adjudicable without
reference to the outcome — a scorer failing a synthetic-data direction test, an identifier
collision between two cohorts — even where reporting their consequences required pCR.

---

## 2. What the parent pre-registration actually specified

Recorded here verbatim because the manuscript currently misstates it in eight places (§3.4).

**PRIMARY** (`preregistration.md` lines 67–70):

> PAM50-basal-adjusted logistic regression `glm(pcr ~ signature_score + pam50_basal_score,
> family=binomial)` yields OR per SD of signature_score ≥ 1.3 with Wald p < 0.0083
> (Bonferroni α = 0.05/6 primary tests, see §6).

**SUPPORTING** (lines 71–72): raw AUROC > 0.55 with 95% bootstrap CI lower bound > 0.50.

**Test family** (parent §6.3, line 264): six tests, T1–T6.

Three facts follow, and all three are load-bearing for the response letter:

1. The pre-specified estimator is a **covariate-adjusted odds ratio on the outcome model**.
   It is not an AUROC of any kind, and not a residual of anything.
2. The string `residual` occurs **zero times** in the hash-verified pre-registration; so
   does the stem `resid`. The terms `nested`, `likelihood-ratio`, `TOST` and `equivalence`
   likewise occur zero times.
3. The cohorts named in the pre-registration are **GSE16446, GSE22226, GSE25066** only.
   **GSE32646 and METABRIC do not appear in it.** No analysis on those two cohorts can be
   pre-specified work under the parent document, independent of the estimator question.

The closest textual hook is §5 line 216, "Both adjusted and unadjusted AUROC are reported."
That is a report-alongside item subordinate to the PRIMARY tier three lines earlier, it
adjusts by one covariate rather than four axes, and it puts the outcome on the left-hand
side rather than the signature. It does not support the manuscript's claim.

---

## 3. Corrections to the record, disclosed before any new analysis

### 3.1 Sign-flipped ΔAUROC in the paired-DeLong analysis (author-found)

`resubmission_v2/scripts/lrt_power_and_delong.py` line 67 called
`delong_var(y, pf, pr)` — full-model predictions in the first slot — while `delong_var`
returns `(AUC(first), AUC(second), …)`, and unpacked the result as `(a_red, a_full)`. The
deposited JSON therefore carries each cohort's **full**-model AUROC under the key
`auc_reduced` and vice versa, and every `delta` is sign-flipped. A second defect of the
same class sits at line 48, where the degenerate `var <= 0` branch returns its two AUROCs
in the opposite order to the normal branch; it never fired on these four cohorts.

Corrected values, recomputed with each AUROC taken from its own model's predictions
(`analysis/A1_delong_corrected.py`, output in `results/A1/`):

| Cohort | n | events | published red / full / Δ | **corrected red / full / Δ** | DeLong P | nested-LRT P |
|---|---|---|---|---|---|---|
| GSE16446 | 114 | 16 | 0.6805 / 0.6652 / −0.0153 | **0.6652 / 0.6805 / +0.0153** | 0.4224 | 0.5672 |
| GSE25066 | 488 | 99 | 0.7867 / 0.7869 / +0.0002 | **0.7869 / 0.7867 / −0.0002** | 0.9411 | 0.4733 |
| GSE22226 | 120 | 32 | 0.7805 / 0.7741 / −0.0064 | **0.7741 / 0.7805 / +0.0064** | 0.3806 | 0.6162 |
| GSE32646 | 115 | 27 | 0.7643 / 0.7656 / +0.0013 | **0.7656 / 0.7643 / −0.0013** | 0.4795 | 0.9715 |

Consequences that must be carried into the manuscript:

- The true ΔAUROC range is **−0.0013 to +0.0153**. The manuscript's "ΔAUROC −0.015 to
  0.000" (lines 353, 419) is that range negated.
- The claim "the model AUROC does not increase in any cohort" (line 351) is **false**.
  Adding the signature raises the model AUROC on GSE16446 (+0.0153) and GSE22226 (+0.0064).
- Table 1's **four-axis column is correct** in all four cohorts. Only the "+signature"
  column is wrong, and only for GSE16446 (0.650 → **0.680**) and GSE22226 (0.768 →
  **0.781**); GSE25066 and GSE32646 were already right.
- **No p-value changes.** DeLong P is invariant to the swap (it is two-sided, and only the
  sign of z flips), and the nested-LRT was never affected. The scientific verdict is
  unchanged: every DeLong P ≥ 0.381, every nested-LRT P ≥ 0.473.
- The LRT power curve reproduces the deposited detection arrays exactly, which is the
  evidence that the rerun changed only the labelling.

The corrected wording must therefore say that the model AUROC **changes by less than 0.016
in absolute value in either direction, with no cohort reaching significance** — not that it
fails to increase.

### 3.2 Bootstrap-Jaccard values (author-found)

The manuscript (lines 132–134) reports bootstrap-Jaccard medians of **0.71 / 0.76** for the
sensitivity and resistance components. Those values match no Doxorubicin signature under
any pool or strategy in `frozen_signatures/MANIFEST.json`.

The canonical values for the deployed signature — pool `solid_only`, strategy
`quantile_30_oncokb_all`, 100 resamples, seed 42, n = 674 cell lines — are
**sens = 0.579 (IQR 0.122), res = 0.714 (IQR 0.109)**. The hash-locked pre-registration
states these same values in four separate places (lines 107, 119, 373, 449) and agrees with
the MANIFEST to three decimals, including the per-signature sha256. **There is no
pre-registration-versus-MANIFEST conflict**; earlier revision notes asserting one were
reading a different drug's entry.

No re-run is required: the deposited MANIFEST already answers this. The manuscript's
downstream gloss "roughly a quarter of genes turn over" must also change — at sens = 0.579
the turnover is closer to 40%. **The correction strengthens the paper's own cautionary
thesis**, since lower membership stability is evidence for, not against, it.

### 3.3 Drug-response data source (author-found)

The manuscript (line 125) states that cell lines "were split at the **GDSC2**
Doxorubicin-AUC quantiles". The actual source is **CTRP / CTD²**. The provenance chain is
unambiguous: the on-disk raw file is `Drug_sensitivity_AUC_(CTD^2)_subsetted.csv`, its
columns carry `CTRP:` identifiers (`DOXORUBICIN (CTRP:36599)`), `freeze_signatures.py`
reads only that file, and all ~186 `drug_ccle_col` entries in the MANIFEST are `CTRP:`.
The GDSC2 asset present on disk is replicate-level dose-viability data, not AUC, so it
could not have supplied the quantities described even in principle.

This is a text correction only; **no analysis code or result is affected.**

### 3.4 Misstatement of the pre-specified estimator (author-found)

The manuscript asserts in **five** places that residual-AUROC was the originally
pre-specified estimator: lines 245–246, 269, 326–327, 430, 745–747. All five are false
against §2 above. The assertion at 326–327 is the most serious, because it sits inside
the manuscript's §2.6 — the section that documents its pre-registration discipline — and directly
contradicts that same section's own admission at lines 320–323 that the axis-decomposition
plan was deposited *after* GSE32646 was scored.

A reviewer's comment on estimator ordering takes this claim at face value and reasons from
it. **That comment must therefore not be answered on its own terms**: doing so would relocate a false statement
to a more prominent position. The response will correct the authors' own manuscript §2.6 and supply a
provenance table, framed throughout as the authors correcting themselves.

**Two further locations assert a pre-specified equivalence margin.** Line 224: *"we therefore
interpret it as an equivalence statement against a pre-specified margin"*. Line 260: *"a
two-one-sided-test (TOST) equivalence (Schuirmann 1987) against the pre-specified margin
(residual AUROC < 0.55…)"*. These are a different kind of error from the five above and are
corrected differently — see §5.1. The number **is** in the deposit, in the unlocked
`PLAN_v2.md`, but as a *falsifier lower* bound rather than an equivalence *upper* bound. The
correction is "mis-attributed and direction-reversed", not "invented".

**An eighth location misstates a different commitment.** Line 335: *"A pre-specified
secondary prediction — a cytotoxic-versus-TKI AUROC gap > 0.3 — was not supported (observed
Δ = −0.136)"*. The word `cytotoxic` occurs **zero times** in the parent, and its H4/T6
prediction is not a gap at all but a ceiling: *max wrong-MoA AUROC < 0.6*. Neither the
threshold 0.3 nor the contrast is pre-registered.

That brings the count of misstated-provenance locations to **eight**: lines 224, 245–246,
260, 269, 326–327, 335, 430, 745–747. If a further location of the same kind is found before
submission it is corrected on the same terms and logged in `decisions_addendum.md`. All line numbers refer to the 943-line `main.tex`,
sha256 `1f69058b…` — the 2026-06-28 build that carries the reviewed title. They do not hold
in the 907-line rival build (`d88d0c4d…`), and each is quoted verbatim in the response letter
so the reference survives any renumbering.

Neither residual-AUROC nor the nested-model LRT was pre-specified, and neither was the
equivalence margin. All three are post-hoc. §5 sets out how they will be ordered and
labelled from here.

### 3.5 Pre-registered commitments not delivered and not disclosed (author-found)

**Ablation grid.** §8b (lines 345–361) the parent commits to a five-axis grid: quantile threshold ×
top-N × gene universe × up:down asymmetry {30:30, 30:10, 10:30} × significance gate {none,
FDR<0.10, FDR<0.05} = **324 cells**, described as the "full grid". What was delivered is
**36 cells** (the first three axes only). The deposited `ablation_grid.tsv` has no
significance-gate column, and every row is symmetric, so the asymmetry and gate arms were
never run. The supplement nonetheless labels the 36-cell table "full". The top-quartile
robustness criterion at line 359 was consequently evaluated against 36 configurations
rather than the 324 it was written for. Analysis **A7** completes the missing arms.

**Decision-Curve Analysis.** The parent's §6.4 (lines 285–300) commits to DCA per primary hypothesis
over a 5–40% threshold range, as a three-state qualifier on H1 (CONFIRMED /
SUPPORTED-ONLY / NOT SUPPORTED). No DCA exists anywhere in the manuscript, the supplement,
or the repository. It is formally moot under the pre-registration's own line 297, since the
PRIMARY criterion failed and that line makes DCA irrelevant in that branch — but §11
requires the deviation to be disclosed, and it was not. It will be disclosed, with the
mootness argument stated rather than assumed.

**A third parent commitment was never delivered and never disclosed: the H4 / T6
wrong-mechanism negative control.** The parent commits to it at lines 136–141 and makes it
test T6 of the locked family at line 273, with the rule that a wrong-MoA AUROC ≥ 0.6
invalidates H1. The manuscript mentions **none** of dasatinib, erlotinib, gefitinib or
"wrong-MoA" anywhere — zero matches. Since the parent's own line 275 makes T6 a constraint
on H1 rather than a hypothesis, and H1's PRIMARY criterion failed, the omission is arguably
moot on the same reasoning as the DCA above — but the parent's §11 still required disclosure, and there
was none. It is disclosed, with the mootness argument stated rather than assumed.

Two further parent commitments are **not yet verified as delivered** and will be checked
before submission: the parent's §6.1 calibration and threshold-metric reporting (Brier, isotonic
slope/intercept, 10-bin reliability, PR-AUC, F1, precision, recall at Youden-J, required in
the primary table rather than the supplement) and the parent §7 CODE-AE head-to-head, which is a
REQUIRED comparator carrying a verbatim failure sentence if it could not be reproduced.

### 3.6 The paired-DeLong companion test is not valid evidence (author-found)

Independent of the label swap, and surviving its correction: both AUROCs in the paired
DeLong comparison come from models fitted on the same data, and the models are nested.
DeLong's test is not valid in that setting (Demler, Pencina & D'Agostino, *Stat Med* 2012).
Measured by permutation on these cohorts, its empirical type-I error at α = 0.05 is 0.004
(GSE16446) and 0.000 (GSE25066); its power against a planted axis-orthogonal effect at
AUROC 0.65 is 0.142, against 0.502 for the nested LRT.

The manuscript's supplementary §S12b presents this test as "the TRIPOD-standard companion
to the nested LRT" showing the signature "does not improve the four-axis model AUROC in any
cohort". Non-significance there is close to guaranteed by construction, so the test cannot
corroborate the null — the reasoning is circular.

**Correction**: ΔAUROC is reported as a descriptive point estimate with a 95% CI; the
DeLong P values are retained only in a footnote carrying this caveat; the inference rests on
the nested LRT, which is valid for nested models. Computed CIs: GSE16446 +0.0153
[−0.0221, +0.0527]; GSE25066 −0.0002 [−0.0050, +0.0046]; GSE22226 +0.0064
[−0.0079, +0.0207]; GSE32646 −0.0013 [−0.0048, +0.0022].

Note that these CIs are wide enough to be worth stating plainly: on GSE16446 the data are
compatible with a ΔAUROC as large as +0.05. That is consistent with, and reinforces, the
bounded conclusion in §10 — large added value is excluded, small increments are not.

### 3.7 The deployed bidirectional score is not the score the Methods describe (author-found)

`resubmission_v2/scripts/phase_6/_common.py` lines 73–96, and thirteen further working
copies plus fourteen deposited copies under `upload/zenodo_v3/`, compute

```
down_score = normalize((1 - rank[down]).mean(), n_down)   # already reversed
raw        = up_score - down_score                        # reversed a second time
```

Subtracting a score already built on reversed ranks is a double negation. With equal set
sizes it reduces algebraically to `raw ∝ r_up + r_dn` — an **undirected** gene-set score,
not the bidirectional singscore described at main.tex lines 127–129 and 138.

Verified two ways. On synthetic data where the answer is fixed by construction
(`analysis/A0_scorer_validation.py`), the sample with up-genes high and down-genes low —
the maximally signature-positive sample — scores **0.0003**, while the sample with *both*
sets high scores **1.9906**. Correlation with the intended directional contrast is
**r = 0.00000000**. On real data, the deposited `score_DOXORUBICIN.tsv` (GSE25066, n = 488)
correlates **r = +1.00000000** with the undirected sum and +0.282 with the contrast.

**Five of the ten grid signatures are affected** — Framework, DoxVariant, MammaPrint, GGI,
RSlike, all of which route through `compute_singscore`. The other five are not: TIS,
cytolytic, TGFBstrom and **HRD** use a plain mean z-score with all genes in one direction,
and DLDA30 uses `mean(UP z) − mean(DOWN z)` with **no rank reversal**, so its subtraction is
correct. That containment is what makes this a code slip propagated by copying rather than
a stated convention, and it is described that way.

A diagnostic that should have caught it earlier is reported rather than omitted: under the
deployed score, MammaPrint (0.407), GGI (0.438) and RSlike (0.315) all had marginal AUROCs
**below chance**. Three of the best-validated breast-cancer signatures scoring
anti-predictively — GGI is a proliferation index, and high proliferation predicts *higher*
pCR — was evidence the scorer was broken.

**Correction**: `raw = up_score + down_score`. All affected cells are rescored, and every
§6 verdict class is computed on the corrected score (deployed figures appear only as a
labelled sensitivity arm). The direction of the consequence is stated plainly: the nested-LRT
null survives on every affected cell, and the redundancy argument **strengthens** — the three
published clinical signatures move from below chance to marginal AUROC 0.72–0.76 with R²
0.766–0.890 on the four axes, which is what redundancy looks like, while the CCLE-derived
signatures fall to R² 0.034–0.302 with no marginal signal. What changes is that the paper now
carries two findings rather than one, and the equivalence margin needs the treatment in §5.

### 3.8 Pre-registration presentation defects (author-found)

Cryptographic integrity is sound: `preregistration.md` re-hashes to the value in
`preregistration.LOCK`, and all four deposited copies are identical. Three cosmetic defects
will be documented in a `PREREG_README.md` accompanying the deposit, **without editing the
hashed file**:

1. Line 5 still reads `Status: DRAFT`, while §9 and §10b both assert "locked 2026-05-28".
2. The parent's §11.1 "Self-lock" contains only forward-references to the LOCK file; its §11 line 469 makes
   filling 11.1 a precondition of being locked, so the document reads as unlocked by its
   own stated criterion while being demonstrably locked by its hash.
3. The parent's §9 cites stale sha256 values for two scripts (`freeze_signatures.py`,
   `verify_frozen.py`) that do not match the deposited files. **All signature-artefact
   hashes match.** This matters because the parent's §9 line 417 states that any artefact whose sha256
   differs from `MANIFEST.json` invalidates the pre-registration — a clause a hostile
   reader could invoke. The README will state that the clause was written for the frozen
   signature artefacts, that those all verify, and that the two script hashes were captured
   before the scripts were extended for later rounds.

---

## 4. Multiplicity, and what counts as a positive

Applies to every analysis in §6.

### 4.1 The analysis unit is a distinct patient

GSE22226 and GSE25066 share **65 patients** — GSE22226's `i-spy_id` and GSE25066's
`sample_id` are the same I-SPY 1 identifier, and every matched GSE25066 row carries
`source == ISPY`. The overlap occupies a different fraction of each analysis set and all
three figures must be quoted, never one alone:

| analysis set | overlapping patients |
|---|---|
| GSE22226, n = 120 | **58** (carrying 11 of its 32 pCR events) |
| GSE25066, n = 488 | **61** |
| GSE25066 HR+, n = 284 (RSlike) | **26** |

**Rule**: in **any** analysis in §6, pooled or per-cohort, no patient contributes more than
once. The de-duplicated GSE25066 (**n = 427, 88 events**) is the analysis set for every use
of GSE25066 in A2, A3, A4, A5, A6, A8, A9, A11 and A13; n = 488 appears only as a labelled
sensitivity arm. The HR+ subset used by RSlike moves 284 → 258 in step.

**The removal set is fixed once and never recomputed per fold.** In every A9
leave-one-cohort-out and discovery-excluded pool — *including the fold that omits
GSE22226* — GSE25066 enters at n = 427. The 61 removed patients are never restored.

Removal is from **GSE25066**, not GSE22226, because GSE25066 enrolled HER2-negative disease
by protocol (485 negative, 6 positive, 4 indeterminate of 495 with a call), so the overlap is
almost entirely GSE22226's HER2-negative part (of the 58: 44 negative, 3 positive, 11
unknown). Removing it from GSE22226 instead would leave a remnant that is 62.5% HER2-positive
against 36.9% in the full analysis set, retains all 11 trastuzumab-treated patients while
removing none, and has nearly double the pCR rate (33.9% vs 19.0%).

Both removal directions are reported. It is disclosed that removal from GSE25066 does shift
clinical T stage — T3 50.8% vs 25.5%, Fisher p = 1.1e-04; T0–T4 omnibus χ² p = 6.4e-04 — a
consequence of I-SPY 1 requiring tumours ≥ 3 cm, even though HER2 status and pCR rate are
unmoved (18.0% vs 20.6%, p = 0.74).

### 4.2 The Benjamini–Hochberg family, enumerated

The family is a fixed list of cells, not a description. Ambiguity here is what lets a family
size be chosen after the fact, so the list is written out and its size fixed **before** the
grid is run:

**Family A — the added-value grid. m = 58, and these are the cells.** Naming a size without
naming the members leaves the members free: thirteen CCLE drugs are scored on all three
discovery cohorts, so "6 drugs" alone admits 1,716 different families that all satisfy
m = 58.

- **Signatures (10)**: Framework, DoxVariant, MammaPrint, GGI, RSlike, DLDA30, TIS,
  cytolytic, HRD, TGFBstrom.
- **pCR cohorts (4)**: GSE16446, GSE22226, GSE25066, GSE32646. → 40 cells.
- **CCLE multi-drug arm (6)**: CARBOPLATIN, CYCLOPHOSPHAMIDE, FLUOROURACIL, PACLITAXEL,
  ETOPOSIDE, TOP2_POISON_CONSENSUS — the deposited artefact names verbatim, all six under
  pool `solid_only`, strategy `quantile_30_oncokb_all` — on discovery cohorts GSE16446, GSE22226, GSE25066. → 18 cells.
- **Excluded, exhaustively.** Wrong-mechanism negative controls: DASATINIB, ERLOTINIB,
  GEFITINIB. Excluded by the parent pre-registration's own ≥ 0.40 bootstrap-Jaccard
  eligibility gate rather than by mechanism: IMATINIB (0.304 / 0.395). Rescalings or composites of a drug already in the grid:
  DOXORUBICIN_fsqn, panDox, DOXORUBICIN_HRpos, DOXORUBICIN_TN, ridge_ridge_rank,
  ridge_ridge_zscore. Agents not administered in any of these regimens: FULVESTRANT,
  LAPATINIB, TRAMETINIB. (DOXORUBICIN itself enters as the Framework signature, not twice.)

**No substitution into or out of these lists after this document is locked.**

- **m = 58 is a fixed denominator.** A cell judged not evaluable enters BH with **p = 1**;
  the denominator stays 58 however many are excluded. Otherwise the realised family size is
  knowable only after the data are seen, and a smaller m is a looser threshold — the
  looseness factor is exactly 58 / m_evaluable.

**Not-evaluable is an exhaustive list, not a judgement.** Left open, it is a one-way silencer:
§4.2 guarantees a not-evaluable cell is never rejected, so an unbounded definition lets an
inconvenient cell be retired. A cell is not evaluable if and only if:

  1. either outcome class has **fewer than 5 events**;
  2. the resolved fraction of the signature's own nominal **up-set or down-set, taken
     separately**, is **< 0.60** — resolution judged against `external_data/hgnc_alias_table.tsv`
     by exact symbol first, then HGNC alias; for a single-direction
     signature, in its one set;
  3. the logistic fit fails to converge or exhibits complete separation, judged
     mechanically: `statsmodels.Logit`, Newton–Raphson, `maxiter` 100, `tol` 1e-8;
     non-convergence is `mle_retvals["converged"] is False`; separation is any
     |coefficient| > 20 or any fitted probability within 1e-8 of 0 or 1. Triggers 1 and 2 are
     settled at lock time, so trigger 3 is the only one an analyst could otherwise reach for
     — hence the numbers.

**No other ground exists.** Every declaration is written to `decisions_addendum.md`, with the
triggering value, **before that cell's p value is computed**. Note that trigger 1 cannot fire
for any Family A cell — the smallest is GSE16446 at 114/16 — so in practice trigger 2 is the
operative one, and the up/down resolved fractions for all 58 cells are recorded in the log at
lock time so the set cannot be revisited later.
**Separate families, each corrected within itself — and each m fixed here, not later.**
Deferring these to `decisions_addendum.md` would leave A4's membership to a gate that needs
pCR (marginal AUROC > 0.55), and §10 condition 3 hangs on A4's BH — so the author would be
deciding after the fact whether condition 3 can fire. Each family is therefore enumerated at
its **maximum** extent now, and a cell that fails its entry gate enters BH at **p = 1** with
the denominator unreduced, exactly as in Family A.

| Family | Maximum cell list | m |
|---|---|---|
| **A4** within-subtype | 10 signatures × 4 pCR cohorts × 3 IHC strata (HR+/HER2−, HER2+, TNBC) | **120** |
| **A7** ablation | quantile {0.20, 0.25, 0.30} × top-N {10, 20, 30, 50, 100} × universe {OncoKB-all, OncoKB-onco-only, protein-coding} × asymmetry {30:30, 30:10, 10:30} × significance gate {none, FDR<0.10, FDR<0.05} | **405** |
| **A11** regimen-matched | 4 pCR cohorts × 2 composites (mean, pre-registered max) | **8** |

PAM50-stratified A4 cells are **exploratory and outside the A4 family**, because which PAM50
strata reach n ≥ 40 is not knowable before the data are read; they are reported with nominal
p only and can never fire condition 3. A cell not evaluable under §4.2 enters its family's BH
at p = 1.
- **The p entering BH is the nested-model LRT p, and only that.** Residual AUROC, marginal
  AUROC and added OR/SD are reported per cell but do not enter the correction.
- **Correction**: Benjamini–Hochberg. **q < 0.10** is the reportable threshold.

### 4.3 Evaluability and reporting

- A cell or stratum meeting a §4.2 not-evaluable trigger is marked *not evaluable*, **enters
  the BH vector at p = 1 per §4.2, and the denominator is never reduced**, and is counted by
  name in the report. This deliberately departs from the parent pre-registration's realised-
  denominator rule at §6.3 line 276; the departure is in the conservative direction and is
  disclosed rather than silent.
- **Commitment**: every cell is reported whatever its result — no cell is dropped for being
  null, and none for being positive. The grid is reported in full as a supplementary table
  even where the main text shows only the heatmap.
- **A positive is a positive.** Any cell reaching BH q < 0.10 is a reportable finding that
  **narrows the abstract's claim**, and will be reported as such rather than explained
  away. This is stated here, before the grid is run, precisely so that it cannot be
  renegotiated afterwards.

---

## 5. Estimator ordering and status

Neither estimator was pre-specified (§3.4). Both are reported, and the ordering is fixed
here in advance:

1. **Residual-AUROC axis-decomposition** — named **first** at every headline location.
2. **Nested-model likelihood-ratio test** — named **second**, co-primary.

Both are labelled **post-hoc** wherever they appear. Every headline location — the
abstract, the first sentence of §3.1, and the first paragraph of §4 — carries both
estimators in the same sentence, plus this disclosure:

> The ordering of estimators was revised after the data were seen; both are reported and
> both return the same verdict.

The parent's literal pre-registered statistic — the PAM50-basal-adjusted logistic OR/SD
Wald test — is additionally run on all five cohorts (analysis A8) so that the one genuinely
pre-specified estimator is reported everywhere, not only where it originally failed.

### 5.1 Equivalence margin

**No equivalence margin was ever pre-specified in the hash-locked pre-registration.** That
document contains no TOST and no equivalence test — `TOST` and `equivalence` occur zero
times, `residual` zero times.

**But the number 0.55 is not invented, and the record must say so.** The same Zenodo record
carries `PLAN_v2.md` (51 KB), whose line 350 reads *"**Pre-specified falsifier (LOCKED):**
residual AUROC 95% CI lower bound > 0.55"*; `residual` occurs 32 times in it and `0.55`
eighteen times. So the manuscript's "pre-specified margin" is **mis-attributed and
direction-reversed**, not fabricated. Three facts fix its true status, and all three are
stated in the revision:

1. **PLAN_v2 was never actually locked.** Its own §2.3 sets out the locking procedure —
   sha256-hash the file, deposit to Zenodo and OSF, obtain a DOI, "NOT a [TBC] placeholder".
   **No `PLAN_v2.LOCK` file exists anywhere in the deposit.** The word "(LOCKED)" on line 350
   therefore describes an intention, not a completed act.
2. **It post-dates the data.** PLAN_v2 §2.2 states that all its §4 analyses are "POST-HOC
   relative to preregistration.LOCK", and the manuscript's own §2.6 (lines 320–323) concedes
   the axis-decomposition plan was deposited *after* GSE32646 was scored.
3. **The estimator was substituted.** PLAN_v2 states the threshold in **both** senses: as a
   falsifier CI-lower bound (L350, §2.6 L264) *and* as a null-confirmation ceiling on the
   point estimate and CI lower bound (§4.8.3 L636; D19 L740, "Residual AUROC ≤ 0.55 → H-b
   confirmed"). The manuscript's TOST upper bound corresponds to the second sense, so the
   residual error is **not** a direction reversal — an earlier draft of this document said it
   was, and that was wrong. The real errors are attribution, to an unlocked post-hoc document,
   and **estimator substitution**: PLAN_v2 defines a per-cohort point-estimate-and-CI rule,
   and the manuscript replaces it with a pooled cross-cohort TOST that PLAN_v2 never specifies.

4. **The α does not carry over either.** PLAN_v2 attaches Bonferroni **α = 0.01** to this
   endpoint (L353, L513, D21 L745). This document uses α = 0.05 because the margin is treated
   as post-hoc and not inherited. That choice is consequential and is stated rather than left
   implicit: the ρ = 1 worst case passes at α = 0.05 (p = 0.0365) and **does not pass at
   α = 0.01**.

The margin is therefore treated as **post-hoc** throughout, and the correction at lines 224
and 260 says "mis-attributed to an unlocked, post-hoc planning document and adjudicated by a
pooled TOST that document does not specify", not "invented" and not "reversed".

The margin is **0.55**, reported across a sensitivity curve of **{0.53, 0.55, 0.575, 0.60}**,
fixed here so it cannot be chosen to fit the result.

**0.53 stays on the curve whatever it returns.** Its pass/fail is adjudicated by the §5.2
cluster bootstrap, which has not been run; no figure currently in hand settles it, and none
is quoted here. If it fails, that failure is reported by the authors as a failure. Deleting the one point on a sensitivity curve that does not pass
— after seeing that it does not pass — leaves a curve on which everything passes, and hides
where the evidence actually stops. That is precisely the degree of freedom this document
exists to remove. **No point may be removed from the curve after the analyses are run;
failing points are reported as failing.** Where the supportable lower bound lies is determined by the §5.2 cluster bootstrap and is
not asserted here; every point on the curve is reported with its p, passing or failing.

An earlier draft of this section justified dropping 0.53 partly by citing a "worst-case
p = 0.13". **That number does not exist in the de-duplicated results.** It comes from
`A0d_pools.json`, a pre-de-duplication pool with a different cell list. It is withdrawn.

Note also that 0.55 appears six times in the parent pre-registration, but **not once as an
equivalence bound**: twice as a threshold to be exceeded (lines 64 and 72, "raw AUROC > 0.55",
the SUPPORTING criterion — both *lower* bounds), once as an interval endpoint (line 13,
"0.55–0.80"), twice as an expected effect size (lines 75 and 338), and once as a false match
inside the string `0.553` (line 146). Using the number as an equivalence *upper* bound is a
new use in a new direction, so the numerical coincidence confers no continuity with the
pre-registration and is not claimed as such.

### 5.2 The pooled estimate is not independent, and this is reported, not corrected away

Even after de-duplication (§4.1), five cells sit on the same GSE25066 patients — RSlike on a
nested HR+ subset of them — and two on the same GSE22226 patients. The naive fixed-effect
pool treats all nine as independent. This is precisely the circularity raised in review, and
de-duplication does **not** fix it.

Every pooled estimate is therefore reported three ways, in this order:

1. **naive fixed-effect** — labelled as assuming independence, which is false;

2. **patient-level cluster bootstrap** — **the primary figure for any equivalence claim.**
   Because §10 condition 7 hangs on this estimator, it is defined here rather than left to
   the analyst: the resampling unit is the distinct patient of §4.1; **B = 2000**; numpy
   `default_rng` seed 42. Each replicate draws patients with replacement from the union of
   the de-duplicated analysis sets; every cell is restricted to the drawn patients and its
   residual AUROC recomputed; a cell with fewer than 5 events in either class is dropped for
   that replicate — the replicate is kept and its realised k recorded, not discarded — and
   the inverse-variance fixed-effect pool is then reassembled. The pooled SE is the standard
   deviation of the 2000 replicate pools, and TOST p = Φ((pooled − margin)/SE_boot).

   Four details that would otherwise be left to the analyst, and change the answer:
   each replicate draws **|union| patients**, the size of the de-duplicated union itself;
   a patient drawn **m times contributes m rows** to each of their cells (the multiplicity
   reading, not the set reading — the two differ by about 26% in the pooled SE); each cell's
   inverse-variance weight is its **observed, un-resampled SE**, held fixed across all
   replicates; and a patient spanning two cohorts is drawn once and enters every cell that
   contains them, never independently per cohort.

   **The same 2000 replicates supply the per-cell CIs.** §6 classes 2, 3 and 4 turn on
   residual-AUROC confidence bounds, so those are the **percentile intervals of that cell's
   2000 replicate values** — not a separate procedure, and not the naive bootstrap used in
   A0b–A0e. Any change to this definition is a §9 deviation.

3. **ρ = 1 worst case** — cells within a cohort treated as perfectly correlated, blocks
   independent: Var = Σ_c (Σ_{i∈c} 1/se_i)² / (Σ_i w_i)², with w_i = 1/se_i². **Not**
   SE_naive × √(k/C), which is the ρ = 1 design effect only under equal cell weights and
   is 26% below the exact value here here; the earlier form is withdrawn. This is an inflation of
   the naive SE, not an independent estimate, and the two must never be presented as
   mutually corroborating figures. The current
   `A0e_dedup.py` computes its "clustered" and "ρ = 1" columns from the same formula and so
   reports one quantity twice; that is a defect in the script, not a second line of evidence,
   and it is corrected before the numbers are used.

And the strictest subset is reported alongside them: **one cell per cohort, all
combinations**, enumerated here so the count cannot drift — Framework/GSE32646; DoxVariant,
GGI, RSlike(HR+), TIS, TGFBstrom/GSE25066; DLDA30/GSE16446; cytolytic, HRD/GSE22226, giving
1 × 5 × 1 × 2 = **10 combinations**. Their TOST is computed on the **naive** fixed-effect
pool, which is appropriate because the four cells in any such combination are already
independent of one another. If A2 changes the pooled cell list, the enumeration is recomputed
and both are shown.

On current evidence **5 of those 10 combinations fail equivalence at 0.55** — TOST one-sided
at **α = 0.05** on the naive fixed-effect pool, the same α and pool named in §10 condition 8
(worst pooled 0.5289, p = 0.186). The α is written here because the count depends on it: at
α = 0.01 six combinations fail, at α = 0.05 five, and at **α = 0.10 only two**, which would
leave condition 8 untriggered. A guard whose threshold can be dissolved by choosing α is not
a guard.

**That count is stated in the manuscript by the authors**, not left for a reviewer to find.

If §10 condition 8 fires, the pre-committed sentence from condition 7 governs the abstract
and every heading. A bounded description — which pools support equivalence at 0.55 and which
do not, stated pool by pool with each p and each assumption, the ρ = 1 pool being the block
form of item 3 — may then appear only as a following sentence, never as the claim itself.
Two pre-specified sentences an author could choose between would be no pre-specification at
all.

---

## 6. Pre-specified revision analyses

Each entry states the question, the data, the model, and what the result would have to look
like to change the paper's conclusion. Analyses are identified by the codes used in
`notes/01_TEN_EXPERT_SYNTHESIS.md` §5.

| ID | Question | Data | Estimator / model | Decision rule |
|---|---|---|---|---|
| **A2** | Does any signature add pCR value beyond the four axes, on any cohort? | 10 signatures × 4 pCR cohorts = 40 cells, plus the CCLE multi-drug arm on the 3 discovery cohorts | Per cell: marginal AUROC; residual-AUROC; nested LRT; added OR/SD with 95% CI; ΔBIC; n, events, gene-resolution fraction | BH q < 0.10 across the family ⇒ reportable positive that narrows the abstract. All 58 cells of Family A reported regardless, together with each NOT-EVALUABLE cell and the trigger that fired. |
| **A3** | Do the conclusions survive a *measured*-clinical baseline, and a stronger one? | GSE25066, GSE22226, GSE32646 (+GSE16446 without ER, uniformly ER−); METABRIC | Reduced model variants: (i) measured IHC ER/PR/HER2; (ii) four axes + grade + clinical T/N stage; (iii) union. Nested LRT vs each | Report base-model AUROC, added OR/SD and LR P for every variant beside the four-axis version. Strengthening the base model can only shrink added value, so a null here is conservative. The grade+stage arm doubles as the real-data positive control for R2-m4. |
| **A4** | Does the signature predict within subtype, or only separate subtypes with different pCR rates? | 4 pCR cohorts; HR+/HER2−, HER2+, TNBC strata; PAM50 strata with n ≥ 40 | Nested LRT within stratum, for the framework and every signature with marginal AUROC > 0.55 | Per stratum: n, events, pCR rate, marginal AUROC with bootstrap CI, added OR/SD, LR P, per-stratum MDE. <5 events ⇒ not evaluable. Report the marginal-minus-within-stratum AUROC gap per cell. |
| **A5** | What are the cohorts actually made of? | All five cohorts' phenotype tables | Tabulation only | New Table 1: n profiled/analysed with exclusion reasons, role, trial and full regimen, platform, verbatim pCR definition, events and rate, IHC ER/PR/HER2, TNBC n, PAM50 distribution, grade, clinical T/N/stage, median age. Includes TRIPOD participant flow and discloses the two near-degenerate axes. |
| **A6** | How much does the unharmonised pCR endpoint matter? | GSE25066 and GSE22226 (RCB class deposited); GSE32646 and GSE16446 unharmonisable | Re-derive RCB-0 and RCB-0/I endpoints; re-run nested LRT and residual-AUROC with all predictors held fixed | Report per cohort: events under each definition, concordance/kappa, added OR/SD and LR P side by side. Plus leave-one-cohort-out of the pooled estimate. State explicitly which cohorts cannot be harmonised and the expected direction of bias. |
| **A7** | Is the mean-difference gene selection stable, and why top 30? | CCLE solid-tumour sensitive vs resistant groups | (a) Per-gene Welch / limma-moderated t-test with BH-FDR; (b) bootstrap selection frequency extended to 1000 resamples; (c) **the two never-delivered pre-registered arms** — significance gate {none, FDR<0.10, FDR<0.05} and asymmetry {30:30, 30:10, 10:30} — plus top-N = 100 | Report how many of the 60 genes survive q<0.10 and q<0.05. Run the nested LRT, not only marginal AUROC, for the best configurations on GSE25066 and GSE32646. Answer "why 30" from the completed grid rather than from the 36-cell subset. |
| **A8** | What does the actually pre-registered statistic say? | All five cohorts | The parent's literal PRIMARY: `glm(pCR ~ signature + PAM50-basal)`, OR/SD, Wald p | **The parent's own two-tier criterion is inherited verbatim, not just its model**: PRIMARY is met only at OR/SD ≥ 1.3 **and** Wald p < 0.0083 (the parent's Bonferroni 0.05/6); SUPPORTING is raw AUROC > 0.55 with bootstrap CI lower bound > 0.50. **Exactly five tests**: the single frozen Framework/DOXORUBICIN signature, once per cohort. A8 sits outside Family A and receives no further correction, and α = 0.0083 is a fixed constant that does not move with the number of evaluable cohorts. Reported on every cohort including GSE32646 and METABRIC, which the parent does not cover — those two are labelled as the pre-registered estimator **extended beyond the parent's scope**, and cannot be described as pre-registered results. METABRIC's endpoint is dichotomised 5-year OS, not pCR; the event is **`dead by 5y`**, patients censored before 5 years are excluded, and the OR is oriented so that a higher signature score predicting *more* death gives OR > 1. Fixing the direction matters: the two codings give reciprocal odds ratios (0.72 vs 1.39), and only one of them clears the 1.3 threshold. The substitution is stated in the same sentence as the result. |
| **A9** | Is the pooled estimate defensible? | The four pCR cohorts | Leave-one-cohort-out and discovery-excluded meta-pools; cohort-clustered (cluster-robust) variance for the residual-AUROC pool; TOST across the §5 margin curve | Report naive and clustered CIs side by side. The pooled row is relabelled "3 of 4 are discovery, not independent — internal consistency, not replication". |
| **A10** | What would it take to exclude a small effect? | Existing spike-in simulation, inverted | Required-n for 80% power to exclude residual AUROC 0.55 and 0.60 at a 23–25% event rate | Report the required n and state whether any public neoadjuvant pCR cohort of that size exists. Converts the power limitation into a study-design recommendation. |
| **A11** | Is the null an artefact of testing the wrong drug? | 4 pCR cohorts; regimens TFAC / AC-T / P-FEC / FEC-docetaxel | Regimen-matched composites (mean, and the pre-registered max over agents actually administered) through the nested LRT | Plus a sensitivity refit excluding the 11 AC-T-Herceptin patients in GSE22226. A positive here would mean the null is drug-attribution, not redundancy — and would be reported as such. |
| **A12** | Why are signatures correlated in patients (Fig 4) but anticorrelated in enrichment (Fig 6)? | Both matrices already on disk; 66 drug pairs | Scatter of NES-profile correlation against deployed-score correlation | Report the association. Label both figures with their space (patient space vs CCLE gene-rank space) and state Figure 6's sign convention explicitly. |
| **A13** | Redundancy, deployment failure, or power? | GSE25066 native algorithmic calls (`ggi_class`, `dlda30_prediction`, `chemosensitivity_prediction`); orthogonal deconvolution immune estimate | (a) Our cohort-free reimplementations vs the original authors' own calls; (b) TIS/cytolytic vs deconvolution | **Stated before running**: high concordance with native calls ⇒ deployment is faithful and the null reflects redundancy; low concordance ⇒ the null is partly a deployment artefact and the abstract must say "as deployed cohort-free". Adds the **six-class** verdict-class column defined below (NOT-EVALUABLE / ADDS-VALUE / INVERSE-DIRECTION / REDUNDANT-WITH-AXES / NO-MARGINAL-SIGNAL / INDETERMINATE), evaluated in that order on the corrected scorer. |

### Verdict-class rule (fixed in advance, used in A2 and A13)

Every clause carries a number, because "collapses toward 0.5" is not a decision rule.
Classes are evaluated **in the order listed** and the first match wins. The order matters:
**ADDS-VALUE is evaluated before NO-MARGINAL-SIGNAL**, because a cell can add value after
axis adjustment while having no marginal signal, and an earlier draft of this rule placed
NO-MARGINAL-SIGNAL first, making ADDS-VALUE unreachable for any cell with marginal AUROC
below 0.60. The reordering is prospective and is not motivated by any cell in hand: under
m = 58 no current cell reaches ADDS-VALUE, since q < 0.10 requires raw p ≤ 0.00172 and the
smallest raw LRT p observed anywhere in A0b/A0c is 0.0695. **m = 58 is not relaxed to make
any cell reachable.**

1. **NOT-EVALUABLE** — fewer than 5 events in either class, **or** fewer than 60% of the
   genes resolved in the up-set or in the down-set *taken separately* (for a
   single-direction signature, in its one set). Reported and counted by name and by which
   clause fired; enters BH at p = 1 per §4.2.
2. **ADDS-VALUE** — nested-LRT BH q < 0.10 **and** the residual-AUROC 95% CI **lower bound
   > 0.5**. The bound is one-sided on purpose: "excludes 0.5" is two-tailed and would file
   an inversely predictive cell as adding value. Narrows the abstract under §4.3.
3. **INVERSE-DIRECTION** — residual-AUROC 95% CI lies entirely **below** 0.5. Counted
   separately, never ADDS-VALUE, and reported — an inverted signature is a finding about
   the signature, not a null.
4. **REDUNDANT-WITH-AXES** — marginal AUROC ≥ 0.60, **and** residual AUROC ≤ 0.55 with its
   95% CI containing 0.5, **and** R²(signature ~ four axes) ≥ 0.50. All three required.
5. **NO-MARGINAL-SIGNAL** — marginal AUROC < 0.60. Redundancy is not testable on that cell,
   and it is **not** evidence for redundancy — it is evidence of nothing.
6. **INDETERMINATE** — anything else.

All CIs here are the patient-level cluster bootstrap percentile intervals defined in §5.2
(B = 2000, seed 42, clustered on the §4.1 patient identifier). R² is computed on that cell's
own cohort, on the de-duplicated analysis set.

**Every verdict class is computed with the corrected bidirectional scorer** (§3.7). Figures
from the deployed scorer appear only as a clearly labelled sensitivity arm and never
determine a class. This matters because the class changes with the scorer — MammaPrint on
GSE25066 is NO-MARGINAL-SIGNAL under the deployed score and REDUNDANT-WITH-AXES under the
corrected one — and the abstract-disclosure trigger above turns on whether
REDUNDANT-WITH-AXES is a majority.

**How BH is applied, stated correctly.** Benjamini–Hochberg is a step-up procedure: the cell
at rank *i* is rejected when p₍ᵢ₎ ≤ i × 0.10/58. The rank-1 threshold is 0.00172, but rank 41
admits p ≤ 0.0707. An earlier draft of this document said q < 0.10 "requires raw p ≤ 0.00172"
and concluded no cell could ever reach ADDS-VALUE; **that was wrong** — it applied the rank-1
threshold to the whole family.

On the nine Family-A cells computed to date the smallest raw LRT p is 0.0695 and no cell is
rejected at its current rank. Whether any of the 58 is rejected is settled by the step-up over
the full family once it is complete, not by this document. The raw LRT p and its nominal-α
verdict are reported for **every** cell alongside the BH result, so a reader preferring a
different correction can apply one. **m = 58 is not relaxed to make any cell reachable.**

Only REDUNDANT-WITH-AXES cells support a redundancy reading. The count in each class is
reported explicitly, **and in the abstract whenever REDUNDANT-WITH-AXES is not a majority of
the 58 cells of Family A**. Where a cell's class differs between the residual-AUROC and nested-LRT
estimators, the LRT governs and the disagreement is flagged and counted.

The R² ≥ 0.50 clause does real work and is not decoration: on the corrected scorer the
published clinical signatures reach R² **0.766–0.890** on the four axes while the
CCLE-derived signatures reach **0.034–0.302**, so 0.50 falls in a genuine gap rather than
being placed to produce a wanted split. The rule separates real redundancy from absence of
signal — a distinction the paper's argument depends on and the previous wording could not
make.

---

## 7. Optional analyses and their gates

**A14 / GSE41998** (n ≈ 279, an additional independent neoadjuvant pCR cohort) —
**OPTIONAL**, with a hard gate fixed here: it must be downloaded, scored and QC-passed by
**2026-08-26**. QC is not a judgement call — it passes only with ≥ 200 patients carrying a
pCR call, ≥ 20 events in each outcome class, and ≥ 0.60 resolution of the signature's up-set
and down-set separately. The QC verdict and its UTC timestamp go into `decisions_addendum.md`
**before any residual AUROC or LRT is computed on the cohort**, so "QC-passed" cannot be
decided after seeing the result. The date and the QC thresholds are the only gates; nothing
about the outcome enters. If that date passes, it is abandoned and **leaves no trace in the
manuscript**. It will not be started and then reported partially, and a null or positive
result cannot be used to decide whether to include it — the date and the QC thresholds are the
only gates.

If it passes, **A14 forms its own BH family**: its cell list and integer m are written to
`decisions_addendum.md` before its first p value is computed, and §10 conditions 1 and 6 apply
to its cells exactly as to Family A cells.

**SCAN-B** — **declined**, as the reviewer explicitly permits, for two reasons only: its
endpoint is survival, so it would reproduce METABRIC's orthogonal-endpoint limitation rather
than remove it; and it cannot narrow a pCR bound because the estimand differs. The
conclusion is instead bounded explicitly by the power floor, with A10 supplying the
required-n figure.

An earlier draft of this section gave a third reason — that at n ≈ 3,000 SCAN-B might yield
a statistically detectable but clinically trivial signal. That reason is **withdrawn**. It
amounts to declining an analysis because it might come out positive, which is not a
defensible ground for omission in a document that is hash-locked and deposited. The first
two reasons stand on their own.

---

## 8. Reporting commitments

1. Every analysis in §6 is reported, whatever it returns.
2. Every Results subsection, table row and figure caption is tagged **[CONFIRMATORY]**,
   **[EXPLORATORY]** or **[PRE-SPECIFIED FOR REVISION]**. Only the parent's T1–T6 may carry
   the first tag.
3. An "Analysis provenance" table appears in the **main text**, one row per analysis, with
   columns: analysis, specifying document, document timestamp, were the relevant data
   already scored at that timestamp, status, verdict.
4. The corrections in §3 appear in a dedicated **"Consistency and corrections"** section of
   the response letter, volunteered rather than defended.
5. The Zenodo deposit is re-minted as a new version with a changelog covering §3, and the
   DOI in the Availability section is updated.
6. A pre-submission QC gate parses every supplementary reference, table and figure
   reference, and every numeric literal in the manuscript, supplement and response letter,
   diffs them against the source JSON/TSV, and **fails the build on mismatch**.

---

## 9. Deviation protocol

Inherited from the parent §11 and binding on this addendum. It covers any deviation from
**§4–§10**, and additionally from each forward-looking commitment made in §3: the
pre-submission delivery check on the parent's §6.1 calibration and threshold metrics and its
§7 CODE-AE comparator; the disclosure of the undelivered DCA together with its mootness
argument; the restriction of the DeLong P values to a caveated footnote; and the
`PREREG_README.md` accompanying the deposit. Those four sit only in §3 and would otherwise
be unprotected.

A deviation requires:

1. **A written rationale appended to `decisions_addendum.md`, deposited alongside this
   document.** The file
   **already exists** — 1102 bytes, sha256
   `4f9409942da6dbb1c09f79c0113f45b6e874de71c1c16a1f83ae8c4a406c467a` — and that value,
   not the hash of an empty file, is recorded verbatim in `ANALYSIS_ADDENDUM_2026-08.LOCK`.
   Its header is the genesis state and is never removed. It is append-only: each entry
   carries a UTC timestamp, the section deviated from, the rationale, and the sha256 of
   every byte preceding that entry's first character. The final file is deposited with the
   revision. **If the file is absent, if its genesis hash does not match, or if a deviation
   has no matching entry, the affected analysis is reported as EXPLORATORY by default and
   the omission is disclosed to the editor.**

   The chain is author-held and therefore self-attesting between the two externally
   witnessed moments — the lock deposit and the revision deposit. That limitation is stated
   rather than papered over; it is what a deviation log can offer without a third-party
   timestamping service.
2. Re-classification of the affected analysis as exploratory.
3. Disclosure in both the manuscript and the response letter.

A deviation discovered after submission is disclosed to the editor rather than left standing.

The range is **§4–§10**, not §4–§7. An earlier draft stopped at §7, which left the reporting
commitments (§8) and the conclusion-change conditions (§10) outside the protocol — the two
sections a motivated author would most want to revise after seeing the results.

---

## 10. What would change the paper's conclusion

Stated in advance so it cannot be adjusted later. Each condition carries a numeric predicate
and a named consequence, because "materially larger", "masked", "flips" and "disagree
substantially" are not decision rules — a condition that cannot be evaluated cannot bind.

The conclusion changes if **any** of the following fires:

| # | Analysis | Trigger (evaluable as written) | Consequence if it fires |
|---|---|---|---|
| 1 | **A2** | Any grid cell is classed **ADDS-VALUE** under §6 — BH q < 0.10 **and** residual-AUROC 95% CI **lower bound > 0.5** | The abstract names that signature and cohort explicitly and drops the unqualified "no signature added detectable value" |
| 2 | **A3** | Added OR/SD against the measured-IHC or grade+stage baseline exceeds the four-axis estimate by ≥ 0.15 on the log-OR scale in ≥ 2 of 4 cohorts, **or** any such cell reaches LRT p < 0.05 where the four-axis version did not | The four-axis construct is reported as insufficient, and the headline is restated against the stronger baseline |
| 3 | **A4** | Any stratum with ≥ 5 events in both classes is classed **ADDS-VALUE** under §6 within its own A4 family — BH q < 0.10 **and** within-stratum residual-AUROC 95% CI **lower bound > 0.5** | The claim is restricted to the marginal analysis and the positive stratum is reported in the abstract |
| 4 | **A6** | Under the harmonised RCB-0/I endpoint, the pooled added OR/SD 95% CI excludes 1.0, **or** ≥ 2 cohorts change verdict class | The endpoint definition is reported as driving the result and the claim is bounded to the deposited call |
| 5 | **A11** | Any regimen-matched composite reaches BH q < 0.10 where the single-agent Doxorubicin signature did not | The null is attributed to drug mis-attribution rather than redundancy, and the title's scope narrows |
| 6 | **A13** | Concordance between our cohort-free reimplementation and the original authors' own deposited calls falls below **κ = 0.60** (or Spearman ρ < 0.60 for continuous calls) on any signature | Every claim for that signature is restricted to "as deployed cohort-free", in the abstract and not only the Discussion |
| 7 | **Equivalence** | The pooled residual AUROC in the PRIMARY (remove-from-GSE25066) pool fails TOST at margin 0.55 under the patient-level cluster bootstrap of §5.2 (one-sided, α = 0.05) | The equivalence claim is withdrawn entirely and replaced by A10's power-bound statement. **Pre-committed fallback sentence**: *"We can exclude a large added value; we cannot establish equivalence at 0.55. The smallest margin the present data support is [value from the curve], reported as such."* **[value from the curve]** is the smallest margin in {0.53, 0.55, 0.575, 0.60} at which the PRIMARY pool passes TOST one-sided at α = 0.05 under the §5.2 cluster bootstrap — that pool and no other. If no point on the curve passes, the sentence says so instead of naming a value. |
| 8 | **Equivalence, minimal subsets** | **≥ 3 of the 10** one-cell-per-cohort combinations enumerated in §5.2 — **or ≥ 30% of them, rounded up, whichever is reached first** if A2 changes the cell list — fail TOST at margin 0.55, **one-sided, α = 0.05, naive fixed-effect pool** | The equivalence claim may not appear in the abstract or in any heading. The condition-7 fallback sentence is used instead, accompanied by the count of failing combinations. |

If A2 changes the pooled cell list, the re-enumeration and its per-combination TOST p values
are written to `decisions_addendum.md` **before any of them is inspected**.

**Condition 8, not condition 7, is the one that has already fired.** On current evidence
5 of the 10 minimal subsets fail (worst p = 0.186), so condition 8 does. An earlier draft had
only condition 7, which left the one result that actually fails with no trigger attached to
it. Writing both here, before the final pools are run, is the point of this section.

**Condition 8 has fired, and that firing is recorded here and not revisited.** The ten
one-cell-per-cohort TOSTs of §5.2 are in hand, the four cells in any such combination are
mutually independent, and they are not a bootstrap quantity — so nothing below defers them.
If A2 changes the cell list, the re-enumeration is reported alongside; it may add failing
combinations but **cannot retract the firing**.

**Condition 7 has not been evaluated.** The §5.2 cluster bootstrap has not been run, and no
figure in this document, in the manuscript, or in `results/A0/` estimates it.
`A0e_dedup.py` computes its ρ = 1 column as SE_naive × √(k/C), which is the ρ = 1 design
effect only under **equal cell weights**; the weights here are unequal (1/se ranges 13.7 to
33.0), and the exact ρ = 1 pooled SE is 36% larger — 0.02770 against 0.02041 on the PRIMARY
pool. Correctly computed, the worst case at margin 0.55 gives p = 0.0365, not 0.0075: close
to α = 0.05, not comfortably past it. **Apart from those ten subset TOSTs, no pass or fail on the §5.1
margin curve is adjudicated by any figure currently in hand.** Every such verdict is recorded
only after the §5.2 bootstrap runs.

Absent all eight, the supportable conclusion is: **large added value is excluded; clinically
small increments are not.**

---

## 11. Estimated scope

Twelve analyses (A0, A1, A2–A13), two of which (A2, A3) are the long pole. Most reuse existing scored
matrices — the expensive cohort-free scoring step is already done for six additional drugs
on three cohorts each, and the within-stratum machinery already exists. A14 is optional
under the §7 date gate.

---

## 12. Lock

- [ ] §3 corrections verified against source artefacts by an independent check
- [ ] §6 analysis list final; no analysis added after this point without a §9 deviation entry
- [ ] sha256 of this document recorded in `ANALYSIS_ADDENDUM_2026-08.LOCK`
- [ ] sha256 of `decisions_addendum.md` at lock time recorded in the same LOCK file
      (currently `4f9409942da6dbb1c09f79c0113f45b6e874de71c1c16a1f83ae8c4a406c467a`,
      1102 bytes — the genesis state, not an empty file)
- [ ] up/down gene-resolution fraction for all 58 Family A cells written to
      `gene_resolution_58.tsv`, its sha256 recorded in the LOCK file, so the not-evaluable
      set of §4.2 cannot be revisited later
- [ ] A4 / A7 / A11 cell lists and integer m final **in this document**, not deferred
- [ ] UTC lock timestamp recorded
- [ ] Deposited to Zenodo as a new version of the existing record
      **10.5281/zenodo.20726403**, alongside the parent pre-registration and its LOCK file,
      **before any A2–A13 analysis is run**
- [ ] DOI recorded in the LOCK file and in the manuscript's Availability section
- [ ] `results/A0/` file manifest with sha256 recorded in the LOCK file, so §1's
      pre-lock exemption can be falsified by a third party

Until every box is checked, this document is a draft and confers no pre-specification
claim. It must not be cited in the manuscript or response letter before it is locked.

**The artefact that is hashed and deposited is this document with the boxes unchecked and no
DOI in it.** The completed checklist, the UTC lock timestamp and the DOI live only in
`ANALYSIS_ADDENDUM_2026-08.LOCK` and the Zenodo metadata — otherwise checking a box would
alter the file whose hash the box records.

---

## 13. Amendments, 2026-08-01

Written **before** the lock and before any A2–A13 analysis is run. Each amendment states what
changed and why. Where an amendment narrows or removes an analysis, that is recorded here rather
than left as a silent non-delivery, per §9.

### 13.1 Gene resolution: the alias table is defective; exact matching is primary

`external_data/hgnc_alias_table.tsv` maps four currently-approved HGNC symbols that are members of
the frozen signature — `HGF`, `SMO`, `CDH1`, `MET` — onto unrelated genes (`IL6`, `SMOX`, `FZR1`,
`SLTM`), and **605 of its rows have an alias field that is itself another gene's canonical
symbol**. `A0f.resolve()`, `A0e.resolve()` and the deployed `phase_6/_common.py:resolve_genes`
all perform an unguarded forward lookup and therefore resolve those four onto the wrong columns.

**Fixed before any analysis runs:** a row `X→Y` is applied only when `X` is not itself a canonical
symbol anywhere in the table. Genuine deprecations that MammaPrint and TGFBstrom require are
retained; the four collisions are blocked. `HGF/SMO/CDH1/MET` are reported as **unresolvable in
these cohorts**, because the matrix builder collapsed their probes into another gene's column, and
`results/A0/gene_resolution_58.tsv` is corrected accordingly.

**Disclosed:** the published Table 2 cells were scored with the four wrong genes in place. This is
a deployment-fidelity correction independent of C1 and is added to the response letter's
corrections section. A second defect in `A0f.resolve()` — it counts query symbols rather than
resolved columns, inflating the recorded fraction when two queries alias to one column — is fixed
at the same time.

### 13.2 A3 is three arms, not one

The A3 specified in §6 is a **robustness** check (does the signature add value over a stronger
baseline). It is not a positive control for the instrument, and §6's sentence calling the
grade+stage arm "the real-data positive control for R2-m4" conflated two different nested tests.
A3 is therefore split, and only A3-r is governed by §10 condition 2:

- **A3-r** — signature added-value against each strengthened baseline. As specified in §6.
- **A3-pc1** — instrument positive control, non-transcriptomic: `LRT[axes + grade + cT + cN]`
  vs `[axes]`. Certifies the machinery on a real covariate, but in an *easy* regime (near-orthogonal
  to the base model). The manuscript must say so in the same sentence as the result.
- **A3-pc2** — instrument positive control, transcriptomic and near-collinear: leave-one-axis-out.
  Reduced model = 3 of {P, H, B, PAM50-basal}; test whether the held-out axis is recovered by the
  identical nested LRT and residual-AUROC decomposition. This is the regime the real test operates
  in and is the more informative control.

**A3-pc1 and A3-pc2 sit outside every BH family** and are reported with nominal p.

**Cohort sets, fixed now.** A3-pc1 and A3-r(ii): GSE25066, GSE32646, GSE16446. **GSE32646 is the
pre-specified primary control cohort** — it is 100 % complete on grade, cT, cN and IHC ER/PR/HER2,
and it is the only genuinely independent cohort. **GSE22226 deposits no nodal field of any kind**,
so the grade+stage arm cannot run there; §6's inclusion of it is a deviation logged under §9.
A3-r(i) (measured IHC): GSE25066, GSE22226, GSE32646 — not GSE16446, which deposits no ER or PR
IHC. On GSE25066 the measured-HER2 term is near-degenerate and the arm is effectively ER/PR-only.

**Covariate coding, fixed now** so the control is the same control in every cohort: grade ordinal
1 df; cT ordinal 1 df; **cN harmonised to binary (N0 vs N+) 1 df in every cohort**, because
GSE32646 records nodal status as binary only. The categorical version is a secondary sensitivity
arm. Family enumeration for A3-r: variants × cohorts as listed above, integer m written to
`decisions_addendum.md` before the first p value.

### 13.3 Interpretation rule for the positive control, fixed before it runs

Each block (grade / cT / cN) is evaluated separately, first match wins. **Grade is
transcriptome-shadowed** — GGI was built to reproduce histological grade from proliferation genes —
so grade failing to add beyond the P axis is predicted by this paper's own thesis and is not
evidence about the instrument. Only the cT/cN block is a clean instrument test.

**Step 1, always reported:** each block's *unadjusted* association with pCR (OR per unit with CI,
block AUROC with CI). This is the planted effect size and it makes the rest interpretable.

**Step 2, the verdict:**
1. **CONTROL UNINFORMATIVE** — unadjusted association null. No conclusion about the instrument.
   A null nested LRT here is arithmetic, not evidence. Never reported as a passed control.
2. **INSTRUMENT VALIDATED** — unadjusted association non-null and nested LRT P < 0.05. Report the
   added OR/SD and ΔAUROC as the **reference scale**: every signature's null is then read as
   "smaller than an effect of this size".
3. **REDUNDANCY CONFIRMED** — association non-null, LRT does not fire, and R²(block ~ axes) ≥ 0.50.
   The base model absorbed the clinical variable. This *extends* the paper's thesis and is reported
   as a finding, not a footnote. Expected for grade; not expected for cT/cN.
4. **INSTRUMENT UNDERPOWERED** — association non-null, LRT does not fire, R² < 0.50. A genuinely
   orthogonal, genuinely predictive block was not detected. **Consequence, fixed now:** "redundant"
   and its cognates may not carry the conclusion; the abstract and Discussion are restricted to
   "no *detectable* added value within a stated power floor", with the per-cohort MDE from A10
   printed next to the headline.

**Step 3, cross-cohort:** GSE32646 has ~27 events, GSE25066 ~99. A fire in the larger and a miss in
the smaller is what an event-count-driven power story predicts and is reported as such — it is not
"the control fired".

### 13.4 A random-signature null was considered and is **not** adopted

A Venet-style empirical null over size- and expression-matched random gene sets was specified in
draft and is declined, for three reasons recorded here so the omission cannot be mistaken for a
result that was run and withheld: (i) no reviewer point requires it; (ii) it would displace A6 and
A11, which answer reviewer points explicitly; (iii) adding an analysis to this pre-registration
after nine of the Family-A cells are already computed is structurally the manoeuvre R2-M2 objects
to, and the disclosure cost exceeds the evidential gain. A pilot also showed the result is not
robust to the null's construction: a signature's percentile moves from the 100th to the 9th
depending on whether random sets are matched on gene–gene coherence, because real gene sets are
co-regulated and independently drawn sets are not. **No random-signature null is run or reported.**

### 13.5 GSE41998 is confirmed in scope under the §7 gate

The optional second independent neoadjuvant pCR cohort is adopted, under §7's existing
**outcome-blind** gate unchanged: ≥200 pCR calls, ≥20 events per class, ≥0.60 up-set and down-set
resolution, verdict and UTC timestamp written to `decisions_addendum.md` **before any residual
AUROC or LRT is computed**, hard date 2026-08-26, abandoned without trace if the date passes.
**Added requirement:** a documented patient-overlap check against all four existing cohorts,
reported whatever it finds, before the gate verdict is written.

### 13.6 A fourth provenance tag, for A8

§8.2's three-way tagging cannot express A8's status. A8 re-runs the pre-registered T1–T6 family on
a corrected predictor. Calling it confirmatory is indefensible; calling it exploratory is too harsh,
because the estimator, the two-tier criterion, α and the family were all fixed pre-data and are
unchanged, and the scorer defect was diagnosed on **outcome-blind synthetic data** before the
re-run — so the correction consumed no researcher degrees of freedom. Fourth tag:

> **[PRE-SPECIFIED ESTIMATOR — CORRECTED EXECUTION]** — the pre-registered test, re-run once on a
> predictor whose defect was established on outcome-blind data before the re-run, with estimator,
> threshold, α and family unchanged.

Both runs are reported side by side. **Pre-committed:** the re-run happens once and the verdict rule
does not move. A change on the OR/SD arm alone, without the Wald-P arm, is a reportable finding and
not a verdict change.

### 13.7 Per-cell TOST is the primary equivalence statement

The pooled TOST is an equivalence statement about the *average across ten different signatures* and
implies nothing about any one of them. It is also one-sided against a 0.55 margin, so an inversely
predictive cell earns full equivalence credit while sitting as far from 0.5 as a predictive one —
the same asymmetry §6's `INVERSE-DIRECTION` class exists to prevent at the verdict level.

**Added:** per-cell TOST at 0.55 for every pooled cell, reported as the primary equivalence
statement, plus a **folded** pool on |residual − 0.5| as a sensitivity arm. The pooled TOST is
retained but relabelled. Values are computed from
`results/A0/A0e_dedup.json → arms["PRIMARY: remove from GSE25066"]`; `A0d_pools.json` is
pre-de-duplication and inflates SEs by √k, and is not used for any reported quantity.

### 13.8 Ordering: the cluster bootstrap runs first

§6's verdict-class rule defines every CI as a patient-level cluster-bootstrap percentile interval
(B = 2000, seed 42). Any A2 or A13 interval produced before that bootstrap exists is
non-conforming, and §10 condition 7's `[value from the curve]` is itself a bootstrap quantity.
Run order: **cluster bootstrap → A13 → A3 → A2 → A8 → A9 → A11 → A10 → A6 → A12 → DCA →
GSE41998 (gate by 2026-08-26)**.

**Pre-committed for the bootstrap:** if `SE_boot > SE_ρ=1`, the ρ=1 column is relabelled "an upper
bound under the block-exchangeable model only", the manuscript states that the bootstrap exceeded
it, and a **stratified** cluster bootstrap is reported alongside as a sensitivity arm. The
unstratified bootstrap of §5.2 governs condition 7 regardless; switching to the stratified version
as primary would be re-choosing the estimator after seeing the answer.

### 13.9 A2 additions

- **Wrong-mechanism negative-control arm.** DASATINIB, ERLOTINIB and GEFITINIB are run on the same
  four cohorts as a labelled negative-control arm **outside Family A**, reported with nominal p.
  Excluding them left the grid with no internal calibration. This also discharges the undisclosed
  T6 non-delivery.
- **Per-cell metadata, printed in the supplementary grid:** scorer type
  (`rank-singscore-bidirectional` / `rank-singscore-unidirectional` / `per-gene-z-mean`),
  `cohort_dependent` yes/no, nominal and resolved set sizes per direction, and
  `deployment_fidelity` (as-published / approximation, with the approximation named). Five of the
  nine published signatures are scored by cohort-wide per-gene z-scoring and are **not** cohort-free
  in the sense §2.1 of the manuscript claims; the grid must not mix the two silently.
- **Trigger scope widened.** §10 condition 1 fires on ADDS-VALUE **and on INVERSE-DIRECTION**.
- **After A2:** the pooled cell list and the one-cell-per-cohort enumeration are re-derived from the
  completed grid, and every per-combination TOST p is written to `decisions_addendum.md`
  **before any of them is inspected**.
- **Validation gate, corrected.** The diagonal cells reproduce `A0h_table2_corrected.json` to ±0.001
  *at Table 2's own n and bootstrap*; only then is the analysis switched to the de-duplicated sets
  and the cluster bootstrap, and the delta reported. The original "reproduce Table 2 exactly" gate
  could not pass: Table 2 uses n = 488, pairs MammaPrint with a survival endpoint, and uses a
  simple bootstrap.
- **RS-like × GSE16446** has no evaluable patients (the score is defined for HR+ only; the cohort is
  uniformly ER-negative) and this ground is not in §4.2's exhaustive not-evaluable list. The cell is
  entered at p = 1 with the reason logged as a stratum-definition vacancy, as a §9 deviation. The
  not-evaluable list is **not** widened.

### 13.10 A13 amended

κ is the wrong instrument for four of the five deposited native calls. `set_class` is 87 % single
class, so κ ≥ 0.60 would require ~91 % agreement and the "unfaithful deployment" verdict would fire
by arithmetic. `chemosensitivity_prediction`, `set_class` and `rcb_0_i_prediction` have no
corresponding re-implementation in the panel at all, so comparing to them measures construct
difference, not deployment fidelity. `dlda30_prediction` is partly in-sample for this cohort.

- **Primary statistic: AUROC of the continuous re-implementation against the native binary call**
  (cut-free and prevalence-free), threshold ≥ 0.80 for "faithful", with cluster-bootstrap CI.
  Secondary: Spearman ρ for ordinal calls, and κ at the prevalence-matched cut, never reported alone.
- **The fidelity verdict is restricted to `ggi_class` and `dlda30_prediction`** (the latter flagged
  as partly in-sample). The other three are reported as descriptive construct comparisons and fire
  no verdict. §6's three-call list is narrowed accordingly; `set_class` and `rcb_0_i_prediction` are
  not added as verdict-bearing tests.
- **A13b, added:** run the identical nested LRT with each native call as the added term against the
  four axes. This is a published positive predictor deployed by its own originators on its own
  patients — no deployment-fidelity confound, transcriptomic, partially collinear with the axes.
  Branches, fixed now: native fires and re-implementation does not ⇒ the null is **deployment**, and
  "as deployed cohort-free" enters the abstract; neither fires ⇒ the null is **redundancy**,
  confirmed against the strongest available challenge; both fire ⇒ condition 1 fires.
- MammaPrint has no native call here, so A13 cannot bear the MammaPrint fidelity claim; that rests
  on the per-cell `deployment_fidelity` column of 13.9.

### 13.11 A4 amendments

`results/A4/A4_within_subtype.json` deviates from §4.2's declared A4 family (10 × 4 × 3 = 120):
it covers three cohorts, and its GSE16446 "TNBC" row is byte-identical to that cohort's marginal
row and is not a stratified analysis. Corrected: **add GSE32646** (100 % complete IHC; the excluded
cohort was the only independent one), **delete the GSE16446 pseudo-stratum**, report **per-stratum
events, EPV, MDE and residual AUROC** in the table itself, apply the §4.1 de-duplication, and print
the simulated type-I error beside any p in the 0.079–0.108 band. Four of the six existing cells have
EPV ≤ 3.2 and one has EPV = 1.0; those nulls are near-zero-power and may not be reported as evidence
of absence.

### 13.12 A6 spec correction

GSE25066's deposited RCB field is **collapsed at RCB-0/I**; RCB-0 is not recoverable. §6's
"re-derive RCB-0 and RCB-0/I" is therefore half-infeasible, and the manuscript's claim that both
cohorts deposit a field from which RCB-0 is derivable is false for GSE25066. Amended: the endpoint
pair is **deposited-pCR vs RCB-0/I**. Two comparisons per cohort, both required — *same-n*
(restricted to the RCB-known subset, isolating the definition effect from attrition) and *full-n*.
Report per cohort: events under each definition, κ between definitions, marginal and residual
AUROC, added OR/SD with CI, LRT P, and whether the six-class verdict changes.

**The two harmonisable cohorts are precisely the two that share 65 patients.** No pooling of them
as independent; any pooled or leave-one-out figure uses the PRIMARY dedup arm; and the manuscript
states that the harmonised analysis is a sensitivity analysis, not a second independent test.
GSE32646 and GSE16446 deposit no RCB and are genuinely unharmonisable — the count of cohorts in
which the question can be asked (2 of 4) is itself the honest answer to R3 point 6.

### 13.13 A7 deferred, disclosed

A7's completed ablation grid (405 cells) is **not run** in this revision. The existing 36-cell grid
is on the pre-correction scorer and must be re-run before any of its numbers are quoted; the
manuscript already concedes the 324-vs-36 non-delivery. Recorded here as a §9 deviation rather than
left silent.

### 13.14 Corrections to this document's own text

- §5's statement that `A0e_dedup.py` computes its clustered and ρ=1 columns from the same formula
  is **stale** — the deposited script fixes that and its docstring says so.
- §3's statement that the manuscript mentions none of dasatinib, erlotinib, gefitinib or
  "wrong-MoA" is **false**: the supplement's T6 row reports a GEFITINIB wrong-MoA AUROC exceeding
  the on-mechanism signature on GSE16446. That is a substantive pre-registered failure and is
  reported, not buried.
- §6's A3 line asserting that "strengthening the base model can only shrink added value" is
  **false**: odds ratios are non-collapsible, so an independent prognostic covariate inflates the
  conditional log-odds-ratio. The direction is not guaranteed a priori, which is why A3 is run
  rather than argued. The same error appears twice in the manuscript and is corrected there.

### 13.15 Binding scope

Every analysis in §6 and in this section runs on the §4.1 de-duplicated analysis sets. n = 488
appears only as a labelled sensitivity arm. This binds A2, A3, A6, A8, A9, A11, A13, the DCA and
the GSE41998 gate — not A2 alone.
