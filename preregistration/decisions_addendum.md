# Deviation log — ANALYSIS_ADDENDUM_2026-08

Required by §9 of `ANALYSIS_ADDENDUM_2026-08.md`. **Append-only.**

Every deviation from §4–§10 of the addendum, or from any forward-looking commitment in §3,
gets one entry below. An analysis that deviates without a matching entry is reported as
EXPLORATORY by default, and the omission is disclosed to the editor.

Each entry records, in this order:

1. UTC timestamp
2. the section deviated from
3. what was done instead, and why
4. the sha256 of this file immediately **before** the entry was appended

This file exists at lock time in its genesis state; that state's sha256, **not the hash of an
empty file**, is recorded in `ANALYSIS_ADDENDUM_2026-08.LOCK`. The header above is never
removed. The file is deposited with the revision.

The 58 per-cell gene-resolution fractions required by §4.2 are **not** written here — they go
to `gene_resolution_58.tsv`, whose own sha256 is recorded in the LOCK file, so that recording
them cannot alter this file's genesis hash.

---

## Entries

Entries appear below in append order; none as of lock.

### 1 — 2026-08-01T23:34:48Z — §6 (A3), §6 (A13), §7, §8.2, §6 (A7), run order

Author decisions taken 2026-08-01, before any A2–A13 analysis was run, and written into the
addendum as §13 before the lock was stamped:

- **§6 A3** split into A3-r (robustness, governed by §10 condition 2), A3-pc1 and A3-pc2
  (instrument positive controls, outside every BH family). §6's description of the grade+stage
  arm as "the real-data positive control" conflated two different nested tests. Cohort sets and
  covariate coding fixed in §13.2; the four-branch interpretation rule fixed in §13.3.
  GSE22226 deposits no nodal field, so the grade+stage arm cannot run there.
- **§7 GSE41998 adopted**, under the existing outcome-blind gate unchanged, plus a documented
  patient-overlap check against all four existing cohorts before the gate verdict is written.
- **A random-signature null was considered and declined** (§13.4). Recorded so the omission
  cannot be mistaken for a result that was run and withheld. No such null is run or reported.
- **§8.2** gains a fourth provenance tag for A8 (§13.6).
- **§6 A13** narrowed: the fidelity verdict is restricted to ggi_class and dlda30_prediction;
  the other three deposited calls fire no verdict (§13.10). A13b added.
- **§6 A7 deferred** and disclosed rather than left silent (§13.13).
- **Run order corrected**: the patient-level cluster bootstrap of §5.2 runs first, because §6's
  verdict-class rule defines every CI as a bootstrap percentile interval (§13.8).
- **Gene resolution**: the HGNC alias table is defective (four frozen-signature symbols aliased
  onto unrelated genes; 605 rows whose alias is itself canonical). Exact matching is primary;
  the alias table is filtered (§13.1). The published Table 2 cells were scored with four wrong
  genes in place — disclosed as a deployment-fidelity correction independent of C1.

Zenodo deposit and DOI remain outstanding author actions; the LOCK carries the timestamp and
the content hashes, and `zenodo_doi` is still `<FILL AT DEPOSIT>`.

sha256 of this file immediately before this entry: `4f9409942da6dbb1c09f79c0113f45b6e874de71c1c16a1f83ae8c4a406c467a`

### 2 — 2026-08-02 — §10 condition 2 FIRED (A3-r), framing deferred to the writing stage

**What fired.** Condition 2, arm B ("any such cell reaches LRT p < 0.05 where the four-axis
version did not"). One cell: **GSE22226, A3-r variant (i), measured-IHC baseline** —
LRT p = 0.0084, added OR/SD 1.919, against the four-axis version's p = 0.1160 on the same cohort.

Arm A did **not** fire: ΔlogOR ≥ 0.15 versus the four-axis estimate occurs in 1 of 4 cohorts
(GSE22226 +0.282; GSE32646 +0.138, below threshold; GSE25066 +0.067), against the "≥ 2 of 4"
requirement.

**The prescribed consequence is not executable as written, and that is disclosed rather than
resolved by weakening the condition.** Condition 2 says "the four-axis construct is reported as
insufficient, and the headline is restated against the stronger baseline". The baseline that
fired is the **weaker** one: measured IHC reaches AUROC 0.710 against the four axes' 0.774 on
the same cohort. It is also a two-marker baseline — GSE22226 deposits ER and HER2 only, no PR —
fitted on n = 101 after 19 patients are lost to missing IHC.

**Context that bears on it, recorded now so it is not assembled after the framing decision:**
- GSE25066's measured-IHC arm is flatly null (p = 0.9950, OR/SD 1.001) on n = 423.
- GSE32646's measured-IHC arm is null (p = 0.2929) on complete three-marker IHC, n = 115.
- The fired cell is in the cohort that **retains** the 65 shared I-SPY 1 patients under the
  PRIMARY dedup arm, and it is the same patients whose GSE25066 counterpart is null.
- Every stronger-baseline arm is null in every cohort: axes+grade+cT+cN p = 0.9089 / 0.8869 /
  0.4351, union p = 0.7059 / 0.9298.

**Reading, stated before the framing is chosen:** the signature adds value over a two-marker IHC
baseline and not over the four transcriptomic axes, because the axes carry proliferation and
basal information that ER/HER2 IHC does not. That supports the design choice, and it is also a
clinically live point a reviewer may press. **The framing is the author's decision and is taken
at the writing stage, not here.** The trigger, the cell and this context go into the manuscript
whichever framing is chosen.

sha256 of this file immediately before this entry: `19377122e5cfb5ddd800a2903e19e203381e81887a6ce164c94e654d1c71e819`

### 3 — 2026-08-02 — A2 complete. Condition 1 does NOT fire; the class-count rule DOES.

Family A ran at its fixed m = 58. **No cell reaches ADDS-VALUE**, so §10 condition 1 does not
fire and the abstract's unqualified "no signature added detectable pCR value" stands.

**But the §6 class-count rule fires**: class counts go into the abstract whenever
REDUNDANT-WITH-AXES is not a majority of the 58 cells. It is not — it is **8 of 58 (14%)**.
The full distribution:

| class | n | % |
|---|---|---|
| NO-MARGINAL-SIGNAL | 33 | 57 |
| INDETERMINATE | 15 | 26 |
| REDUNDANT-WITH-AXES | 8 | 14 |
| NOT-EVALUABLE | 1 | 2 |
| INVERSE-DIRECTION | 1 | 2 |
| ADDS-VALUE | 0 | 0 |

This bears directly on how the paper may argue. The rule itself says NO-MARGINAL-SIGNAL is
**not** evidence for redundancy — "it is evidence of nothing" — so the redundancy narrative is
carried by 8 cells, not 58, and the abstract must say so. The one INVERSE-DIRECTION Family-A
cell is FLUOROURACIL on GSE16446. The one NOT-EVALUABLE cell is RSlike on GSE16446 (stratum
vacancy: the score is defined for HR-positive disease, the cohort is uniformly ER-negative),
entered at p = 1 per §4.2 without widening the not-evaluable list.

**Wrong-mechanism negative controls (outside Family A, §13.9).** GEFITINIB and ERLOTINIB on
GSE16446 reach nominal p = 0.0025 and 0.0046 — but strongly INVERSE (marginal 0.280 / 0.309,
residual 0.286 / 0.305), not adding value. Reported as what they are: a wrong-mechanism control
that is inversely associated with pCR in a uniformly triple-negative cohort. Including this arm
is what makes the grid interpretable at all; a reader can now see the test is not inert and can
see which way a spurious cell goes.

Framing of the class counts in the abstract is the author's decision, taken at the writing
stage. The counts themselves are not negotiable and go in whichever framing is chosen.

sha256 of this file immediately before this entry: `b67f1dc82a5c1d65a3e6338fd23999b7b0f4f05a67bbc192d27bd986e1486633`

### 4 — 2026-08-02T03:38:46.733380+00:00 — §7 GSE41998 gate: **PASS**. Written before any outcome-dependent computation.

The optional second independent neoadjuvant pCR cohort was adopted by author decision
(entry 1) and its outcome-blind gate is evaluated here. **This entry is written before any
residual AUROC, nested LRT or other outcome-dependent quantity is computed on GSE41998.** The
ordering is enforced structurally, not by discipline: `A14a_gse41998_gate.py` contains no code
that computes such a quantity, and `A14b` is gated on this entry existing.

| criterion | value | threshold | verdict |
|---|---|---|---|
| pCR calls | 253 | >= 200 | PASS |
| events, pCR | 69 | >= 20 | PASS |
| events, no pCR | 184 | >= 20 | PASS |
| up-set gene resolution | 0.900 | >= 0.60 | PASS |
| down-set gene resolution | 0.933 | >= 0.60 | PASS |
| patient overlap vs all four cohorts (§13.5) | 0 | 0 | PASS |

**Provenance.** NCBI/GEO returns HTTP 403 to this host — consistent with the server-side flood
lockout recorded earlier — so the series was retrieved from the EBI BioStudies mirror
**E-GEOD-41998** (same submission, 279 samples, verified public). 279/279 per-sample tables
fetched serially with a fixed inter-request delay, zero failures. Platform is A-AFFY-37 =
Affymetrix HG-U133A (GPL96), the same platform as GSE25066, so the already-deposited GPL96
annotation was reused and no further external request was made.

**Overlap check (§13.5), documented rather than asserted.** GSM accession and specimen /
patient identifier set intersection against GSE16446, GSE22226, GSE25066 and GSE32646:
**0 in every pairwise comparison**. Provenance is also disjoint — GSE41998 is a Bristol-Myers
Squibb trial of AC followed by ixabepilone or paclitaxel, against EORTC 10994, I-SPY 1, the
MDACC/USO/LBJ TFAC series and Osaka. This is the same class of check that detected the
65-patient GSE22226/GSE25066 overlap, applied prospectively this time.

**One deliberate departure from how the existing cohorts were built, disclosed here.** The
defective HGNC alias table of §13.1 was **not** applied at matrix-build time. Applying it is
what collapsed HGF, SMO, CDH1 and MET — four approved symbols in the frozen signature — into
unrelated genes in the existing cohorts. All four are present and correctly labelled in
GSE41998, which is why its resolution (0.900 / 0.933) exceeds the existing cohorts'. The
asymmetry is a consequence of correcting the defect and is reported with the result rather
than smoothed over.

Hard date 2026-08-26 was not reached; the gate is settled on 2026-08-02.

sha256 of this file immediately before this entry: `1c1f9358ac04968035a3eb7f1137db4733a1b8f9ad95cb627c2fbaf5e35e3f8d`

### 5 — 2026-08-02 — A14b: GSE41998 result. Three signatures add value on the second independent cohort.

The gate passed (entry 4) and the outcome-dependent analysis was then run. GSE41998 sits
**outside Family A** — section 4.2 fixes that family at 58 cells and it may not be enlarged
after the fact — so these cells carry nominal p and no BH correction. Section 10 condition 1
therefore does not fire. That is a statement about bookkeeping, not about importance.

**The framework signature is null**, which is the paper's headline claim and it survives:
residual AUROC 0.548 [0.469, 0.625], LRT p = 0.321, on n = 253 with 69 events.

**Three other signatures are not null**, with residual-AUROC 95% CI lower bounds above 0.5:

| signature | marginal | residual [95% CI] | R² on axes | LRT p |
|---|---|---|---|---|
| cytolytic (GZMA/PRF1) | 0.671 | 0.623 [0.546, 0.699] | 0.158 | **0.0035** |
| TIS (IFN-γ) | 0.675 | 0.623 [0.547, 0.700] | 0.192 | **0.0046** |
| DoxVariant | 0.620 | 0.611 [0.528, 0.687] | 0.091 | **0.0147** |

Applying BH within this cohort's own 10 tests, the smallest q is ≈ 0.035, so the result is not
an artefact of multiplicity at the within-cohort level either.

**This bears directly on R2-m5.** The reviewer's point was that TILs are established
subtype-independent pCR predictors yet fail here, and that the nulls might reflect deployment or
power rather than redundancy. On a cohort with 253 calls and 69 events, the two immune
signatures do **not** fail. The reviewer's suspicion is supported by our own data.

**The caveat that must travel with it, stated now rather than after the framing is chosen.** The
four-axis base model is markedly weaker on GSE41998 than anywhere else — AUROC **0.673**, against
0.766 (GSE32646), 0.774 (GSE22226) and 0.797 (GSE25066). A weaker base model leaves more
outcome variance for any added term to explain. The direction is not guaranteed a priori
(§13.14: odds ratios are non-collapsible), but it is the first thing a reviewer will ask, and the
low R² on axes for these three cells (0.09–0.19) says the same thing from the other side: these
signatures are close to orthogonal to the axes in this cohort, in a way they are not elsewhere
(GGI 0.923, MammaPrint 0.786, DLDA30 0.890 on the same patients).

**Secondary, and unique to this cohort.** GSE41998 randomises AC followed by ixabepilone versus
paclitaxel — the only randomised within-cohort treatment contrast in the study. The framework
signature is null in both arms (ixabepilone n=132, p=0.985; paclitaxel n=121, p=0.094).

**No result is withdrawn, reordered or reweighted on the basis of this outcome.** The framing —
what enters the abstract, what enters the Discussion, how the immune result is bounded — is the
author's decision at the writing stage. The numbers and this caveat go in whichever framing is
chosen.

sha256 of this file immediately before entry 5: `4961c5a8140c2231576d7a8f70fc230f19f7b0adfbd4b0e55ed7c441fe099075`

### 6 — 2026-08-02 — §13.1 was specified but not implemented in the first execution. Corrected and re-run.

**What happened.** §13.1 requires that an HGNC alias row `X→Y` be applied only when `X` is not
itself a canonical symbol elsewhere in the table, because 605 rows fail that test and four of
them collide with members of the frozen signature (`HGF→IL6`, `SMO→SMOX`, `CDH1→FZR1`,
`MET→SLTM`). That filter was written into the addendum before the lock, but
`A0e_dedup.py:resolve()` was **not** patched to apply it, and A0i, A2, A3, A13 and A14b were
first executed on the unfiltered resolver — two wrong genes per direction in the framework
signature.

**How it was found.** Not by the process. `fig2.pdf`'s residual annotation disagreed with
Table 2 for the same cell, and tracing that disagreement led to the resolver. The drift gate,
the ledger and the five-expert audits would not have caught it; only the numeric cross-check
did, and only because a figure happened to print the value.

**Magnitude.** Under the filter, resolution equals exact matching in all four cohorts — the
filter removes exactly the collided genes and nothing else — and the per-cohort results
reproduce the deposited Table 1 (`A0g`) exactly, which the unfiltered run did not. So the
filtered rule is also the rule Table 1 was always built on; the unfiltered run was the outlier.

**Every conclusion is unchanged. Several numbers moved.**

| quantity | unfiltered (withdrawn) | filtered (of record) |
|---|---|---|
| Framework residual, GSE32646 | 0.5105 | **0.4958** (now agrees with Table 2's 0.496) |
| pooled residual | 0.5003 | **0.4995** |
| pooled TOST at 0.55 | 3.03×10⁻⁴ | **2.53×10⁻⁴** |
| folded pool / its TOST at 0.55 | 0.5302 / 0.0733 | **0.5296 / 0.0668** |
| per-cell equivalence at 0.55 | 2 of 9, both anti-predictive | **unchanged** |
| A2 NO-MARGINAL-SIGNAL / INDETERMINATE / INVERSE-DIRECTION | 33 / 15 / 1 | **28 / 19 / 2** |
| A2 REDUNDANT-WITH-AXES, ADDS-VALUE | 8, none | **unchanged** |
| A3 GSE22226 measured-IHC arm | p = 0.0084 | **p = 0.0162** |
| A14b DoxVariant | p = 0.0147 | **p = 0.0184** |
| A13 concordance, branch verdicts | 0.952 / 0.962, REDUNDANCY | **unchanged** |
| A14b TIS, cytolytic | p = 0.0046, 0.0035 | **unchanged** |

**Trigger status is unchanged.** §10 condition 1 still does not fire (no ADDS-VALUE cell);
condition 2 still fires on the GSE22226 measured-IHC arm; the class-count rule still fires
(REDUNDANT-WITH-AXES is 8 of 58, not a majority). GSE41998 is unaffected because its matrix was
deliberately built without the alias table (§13.5), which is why it was the only cohort whose
numbers did not move.

**Consequence.** `analysis/A0e_dedup.py` is modified, so its sha256 changes and
`ANALYSIS_ADDENDUM_2026-08.LOCK` is regenerated. The addendum text itself is unchanged — §13.1
already said this; the deviation was in execution, not in specification. The unfiltered figures
above are withdrawn and are not reported anywhere.

sha256 of this file immediately before entry 6: `4158abc7a2afd8d56197b8a802a27a15896dbce4b11db03ebf7c636815f9a254`

### 7 — 2026-08-02 — A4 extended per §13.11. One within-stratum cell now reaches nominal significance.

§13.11 requires A4 to add GSE32646, drop GSE16446's pseudo-stratum, report per-stratum EPV and
MDE, and apply the §4.1 de-duplication. Done. GSE22226's cells reproduce the deposited analysis
exactly; GSE25066's shift only by the mandated de-duplication (TNBC 199→167, HR+/HER2− 284→258).

**The manuscript's current sentence — "no within-stratum cell reaches even nominal significance
(six evaluable cells, all P > 0.10)" — is no longer true once the independent cohort is
included.** Of eight evaluable cells, one reaches p < 0.05: **GSE32646 TNBC, residual AUROC
0.762, LR p = 0.0121**.

Everything needed to read it correctly, reported alongside and not in a footnote:

- n = 26 with 10 events, **EPV = 2.0**. Five of the eight evaluable cells have EPV ≤ 3.
- The **minimum detectable added OR/SD at 80 % power in that stratum is 4.6** — only an
  enormous effect is detectable at this size, so a nominal p of 0.012 there is fragile.
- Simulated type-I error at this EPV runs **0.079–0.108** against a nominal 0.05, so the cell is
  not correctly sized either.
- A4's pre-specified family is 120 cells (§4.2). Under BH at m = 120 this p does **not** survive.

So the honest revision of the sentence is that one of eight cells reaches nominal significance,
in the smallest stratum in the study, at an event count where the test is neither adequately
powered nor correctly sized, and it does not survive the pre-specified correction. It is
reported, not suppressed, and the "does not predict within them either" phrasing is withdrawn.

**A methodological error caught before it entered anything.** A first draft of the extension
script wrote its own stratum definition — it used PR in addition to ER and dropped patients with
missing HER2 rather than treating them as negative, as the deposited `A4.strata_of` does. Under
that partition GSE25066's TNBC stratum moved from 199 to 148 and **two** cells reached p < 0.05.
Changing a stratum definition and reporting the resulting significance is precisely the
researcher degree of freedom this paper exists to warn about. The deposited definition is now
imported verbatim rather than re-implemented, and the discarded variant is recorded here so the
choice is visible.

**A5** extended to all five cohorts (R3 point 3 asked for each): GSE32646 n = 115, 27 events,
TNBC 26 / HER2+ 34 / HR+/HER2− 55; METABRIC n = 412 with 151 five-year events, ER+ 158,
HER2+ 103, endpoint dichotomised OS and not pCR.

sha256 of this file immediately before entry 7: `897218d4f969c1dbd6c9c15bda28b273fec718c8d6a5b40738cd2ca961a442a9`

### 8 — 2026-08-04 — §6 A6/A8/A9/A10/A11/A12 and the DCA re-run under the §13.1 resolver; one code defect corrected

Entry 6 recorded that the §13.1 alias filter was implemented and that A0i, A2, A3, A13 and A14b
were re-run under it. **Six further §6 analyses and the decision curve were not**, and their
deposited results remained on the withdrawn unfiltered resolver. This was found by a five-expert
audit of the manuscript, not by the pipeline: `A11_regimen.json` gave LR P = 0.116 for the
GSE22226 framework test where the manuscript's own Table 1 prints 0.16, and
`A9_pooled_defensibility.json` carried the 0.5003 pool that entry 6 itself names as withdrawn.

All seven were re-run once, with no change to any script other than the defect below, and with the
verdict rules untouched. Measured effect of the resolver on each:

| analysis | withdrawn | re-run | verdict change |
|---|---|---|---|
| A6 GSE25066 same-n deposited LR P | 0.5366 | 0.2907 | none |
| A6 GSE22226 full-n LR P | 0.1160 | 0.1553 | none |
| A8 GSE16446 OR/SD (corrected scorer) | 1.339 | 1.484 | none — still fails Wald |
| A8 GSE22226 OR/SD (corrected scorer) | 1.279 | 1.256 | **now below the 1.3 arm**; still fails overall |
| A9 all-cells pooled residual | 0.5003 | 0.4995 | none; now identical to A0i |
| A9 drop-GSE25066 pooled / TOST@0.55 | 0.5201 / 0.161 | 0.5168 / 0.136 | none |
| A10 required n at 0.55 / 0.60 | 560 / 146 | 560 / 146 | unchanged — synthetic simulation, resolver-independent |
| A11 GSE22226 full / excl. trastuzumab | 0.1160 / 0.1049 | 0.1553 / 0.1482 | none |
| A12 r(NES corr, singscore corr) | +0.619 | +0.619 | unchanged |
| DCA max delta net benefit | +0.0065 / +0.0010 | +0.0065 / +0.0010 | none |

**Six of the six pre-registered T1–T6 tests still fail. No firing condition changes state.**

**Code defect corrected, disclosed rather than silently fixed.** `analysis/A_remaining.py:dca()`
drew a bootstrap replicate index `i` and then never used it — each of the 500 replicates recomputed
the point estimate on the full data, so the deposited `delta_ci_lo` and `delta_ci_hi` were exactly
equal to `delta_net_benefit` and the interval was degenerate. The replicate now resamples patients
and refits both the reduced and full models. The point estimates are unaffected. Under the repaired
interval the added net benefit's 95% bootstrap CI **contains zero at all 36 thresholds on both
cohorts**, which is the first time the decision curve carried usable uncertainty at all.

Both changes are corrections of execution, not of specification. Nothing in §6, §8 or §10 is
amended and no verdict rule moved.

### 9 — 2026-08-04 — §13.13: the 36-cell ablation grid re-run on the corrected scorer. A pre-registered check flips PASS to FAIL.

§13.13 states that the existing 36-cell ablation grid is on the pre-correction scorer and
**must be re-run before any of its numbers are quoted**. The manuscript quoted it anyway, in §2.1,
as the justification for choosing 30 genes per direction. The grid is now re-run.

**The original script carried the same C1 defect.** `resubmission_v2/scripts/ablation_grid_completion.py`
computes `raw = norm(up) - norm(1 - down)`, which expands to `r_up + r_dn` — the undirected sum,
not the bidirectional contrast. This is the same sign error found in the deployed scorer, in a
second independent script. The re-run copy, `analysis/A7b_ablation_corrected.py`, is verbatim apart
from that sign and the output path.

| quantity | pre-correction | corrected | 
|---|---|---|
| primary config (q=0.30, N=30, OncoKB-all) marginal AUROC on GSE25066 | 0.5954 | **0.4874** |
| rank within the 36 | 5 | **30** |
| quartile | 0.139 | **0.833** |
| parent pre-registration §8b top-quartile check | **PASS** | **FAIL** |
| best AUROC per N (10 / 20 / 30 / 50) | 0.580 / 0.583 / 0.599 / 0.593 | 0.590 / 0.558 / 0.499 / 0.541 |
| grid range | — | 0.457 to 0.590 |

**Consequence, executed rather than deliberated.** A pre-registered robustness check that fails is
reported as failed. §2.1's three claims built on the old grid — "ranks 6th of 36", "the five above
it differ by ΔAUROC < 0.004", and "30 sits on a flat optimum" — are false on the corrected score
and are withdrawn. The claim that survives is the one that bears on the paper's conclusion rather
than on the choice of 30: no configuration in the grid reaches an AUROC that would change it.

**Direction of the effect on the paper's thesis.** The corrected primary configuration is *weaker*
(0.487, below chance) than the buggy one (0.595), so the correction moves the evidence toward the
null the paper reports, not away from it. That does not make the failed check less of a failure,
and it is reported as a failure regardless of which way it points.

**Scope note that follows from this.** The C1 sign defect has now been found in two independent
scripts — the deployed `phase_6/_common.py` scorer and this ablation script. Any other deposited
artefact computed with a hand-written singscore should be checked for the same expansion before
its numbers are quoted. `SUPPLEMENTARY.md` Table S2 reproduces all 36 pre-correction cells and
must be refreshed from `results/A7b/ablation_grid.tsv`.

**Verification of entry 9, added 2026-08-04 after the author asked whether the new run, not the old one, might be the wrong one.** Three independent checks, none of which relies on the argument that the corrected sign is theoretically right:

1. **Cross-script agreement.** The primary ablation configuration (q=0.30, N=30, OncoKB-all) *is* the frozen framework signature. The re-run gives its marginal AUROC on GSE25066 (n=488) as **0.487393**. `results/A4/A4_within_subtype.json`, produced by a different script under the corrected scorer, gives **0.4874** for the same quantity. The pre-correction value, 0.5954, is reproduced by no corrected analysis anywhere in the deposit.
2. **Diff scope.** `diff` of `analysis/A7b_ablation_corrected.py` against the original returns exactly two changed lines: the singscore sign and `OUT_DIR`. No selector, universe, cohort, threshold or AUROC code was touched.
3. **The sign identity, tested rather than asserted.** On simulated data the original `norm(up) - norm(1 - down)` form correlates with the undirected rank sum `r_up + r_dn` at **Pearson r = 1.0000000000**; the corrected `norm(up) + norm(1 - down)` form correlates with it at 0.026. The original was computing an undirected gene-set score. This is the same diagnostic that identified the C1 defect in the deployed scorer, applied to the ablation script.

The failed top-quartile check therefore stands.

### 10 — 2026-08-04 — A new analysis at revision stage: an outcome-blind criterion for the gene count

Entry 9 withdrew the manuscript's justification for N = 30 (it rested on an ablation rank that
the corrected scorer reverses). That left the gene count unjustified. This entry adds a
replacement and labels it honestly.

**This analysis is NOT in the addendum and is NOT pre-registered.** It is specified and run at
revision stage, and it is labelled as such in §2.1 and in Supplementary §S13. The addendum's
ordering rule — add, lock, then run — is not satisfied for it, and that is disclosed rather than
worked around.

**Why running it anyway is defensible, stated so a reader can disagree.** The criterion is
**outcome-blind by construction**: `analysis/A7c_gene_count_stability.py` reads only the CCLE
expression matrix, the CCLE drug-response table, the OncoKB gene list and the CCLE model
annotation. It never opens a patient cohort, a pCR call or any clinical endpoint, and it cannot
be run against one — there is no code path from it to an outcome. It therefore consumes no
researcher degrees of freedom on the reported result, which is the reason the ordering rule
exists. This is the same basis on which §13.6 treats the A8 scorer correction as consuming none.

**What it does.** Cell lines are resampled with replacement from the CCLE solid-only pool
(674 lines; OncoKB-all universe of 1,202 genes), the DOXORUBICIN signature is re-derived on each
resample under the frozen selection rule, and the bootstrap set is compared with the full-data
set by Jaccard index (100 resamples, seed 42), for N in {5, 10, 15, 20, 25, 30, 40, 50, 75, 100}.

**Validation.** At N = 30 the scan returns Jaccard 0.579 (sensitivity) and 0.714 (resistance),
which are exactly the values deposited for the frozen signature. The machinery reproduces the
frozen pipeline rather than approximating it.

**Result.** Stability rises to roughly N = 15–30 and is flat beyond: resistance component 0.548
(N=5), 0.538 (N=10), 0.667 (N=15–25), 0.714 (N=30), 0.705–0.739 (N=40–100). The random-selection
baseline E[J] = N/(2M−N) runs 0.002–0.043 and is reported alongside, so the curve is not an
artefact of set size. N = 30 is on the plateau; it is distinguishable from N = 5 and N = 10 and
not from N = 40–100.

**The finding worth reporting is the disagreement.** The corrected ablation grid's best cell uses
N = 10, at which this selection rule is the least stable setting tested apart from N = 5. Choosing
the gene count by AUROC on the evaluation cohort would have selected the configuration least
reproducible under resampling. That is the paper's own cautionary thesis applied to its own
hyper-parameter, and it is reported as such rather than as a tuning success.

**What was NOT done, deliberately.** N was not re-chosen. The frozen signature keeps N = 30, which
was committed before any patient data were scored. Re-selecting N now — by any criterion — would
replace a pre-committed value with a post-hoc one and would be a worse methodological position
than the one this entry repairs.

### 11 — 2026-08-04 — §7: the A14 BH family was never enumerated, and entry 5 stated the opposite of §7. Condition 1 fires.

**§7, verbatim.** "If it passes, **A14 forms its own BH family**: its cell list and integer m are
written to `decisions_addendum.md` before its first p value is computed, and §10 conditions 1 and 6
apply to its cells exactly as to Family A cells."

So a family was pre-specified. Entry 5 said the opposite — "these cells carry nominal p and no BH
correction. Section 10 condition 1 therefore does not fire" — and `results/A14/A14b_analysis.json`
records `family = "outside Family A (section 7 optional cohort); nominal p, no BH"`. Half of that is
right: these cells are indeed outside Family A. The inference drawn from it is not. §7 gives them
their own family rather than no family.

**The deviation, stated plainly.** §7 required the cell list and integer m in this log before the
first p value. No entry contains them. Entry 4 recorded the QC gate and stopped there; entry 5 went
straight to results. The ordering §7 imposes was therefore not satisfied, and this entry is the
disclosure rather than a repair — the recording cannot be back-dated.

**What the deviation does not do is make m adjustable.** §4.2 fixes the signature list at ten and
confines the six-drug CCLE arm to the three discovery cohorts, which GSE41998 is not. Ten
signatures on one cohort is the family's maximum extent under §4.2's own enumeration rule, so
**m = 10**, determined by lists locked before the cohort was adopted. There is no researcher
latitude in it, which is the only reason the correction below is reportable at all.

**m = 68 is the one value the locked document forbids.** Folding these cells into Family A
contradicts §4.2 — "m = 58 is a fixed denominator", "No substitution into or out of these lists
after this document is locked" — and entry 5 says so itself ("it may not be enlarged after the
fact") in the same paragraph in which it enlarges nothing but reasons as though the alternative
were nominal p. An intermediate correction computed at m = 68 and concluded that no cell survives;
that computation is withdrawn here as contrary to §4.2.

**BH step-up at m = 10** (`analysis/A14c_bh_family.py` → `results/A14/A14c_bh_family.json`). The p
entering BH is the nested-model LRT p and only that, per §4.2. CIs are carried through unmodified
from A14b (B = 2000, seed 42, §5.2 cluster bootstrap).

| signature | LRT p | BH q | residual [95% CI] | verdict (§6) |
|---|---|---|---|---|
| cytolytic | 0.0035 | **0.0230** | 0.623 [0.546, 0.699] | **ADDS-VALUE** |
| TIS | 0.0046 | **0.0230** | 0.623 [0.547, 0.700] | **ADDS-VALUE** |
| DoxVariant | 0.0184 | **0.0614** | 0.607 [0.525, 0.683] | **ADDS-VALUE** |
| DLDA30 | 0.1200 | 0.2999 | 0.567 [0.489, 0.639] | INDETERMINATE |
| Framework | 0.3210 | 0.5384 | 0.548 [0.469, 0.625] | NO-MARGINAL-SIGNAL |
| RSlike | 0.3265 | 0.5384 | 0.395 [0.205, 0.625] | NO-MARGINAL-SIGNAL |
| HRD | 0.4134 | 0.5384 | 0.478 [0.397, 0.553] | NO-MARGINAL-SIGNAL |
| TGFBstrom | 0.4308 | 0.5384 | 0.534 [0.455, 0.615] | NO-MARGINAL-SIGNAL |
| GGI | 0.6886 | 0.7651 | 0.532 [0.453, 0.615] | REDUNDANT-WITH-AXES |
| MammaPrint | 0.9895 | 0.9895 | 0.500 [0.422, 0.580] | REDUNDANT-WITH-AXES |

**§10 condition 1 fires.** Three cells meet ADDS-VALUE — BH q < 0.10 and residual-AUROC 95% CI
lower bound > 0.5. Its consequence is fixed in advance: the abstract names those signatures and
that cohort explicitly and drops the unqualified "no signature added detectable value".

**The manuscript already discharges it; this entry is the record catching up, not the manuscript.**
The abstract already reports q = 0.023 for cytolytic and TIS and q = 0.061 for the Dox-variant
"under the family the analysis plan specifies for it", and §4 already states m = 10, classes the
three cells ADDS-VALUE, and rejects m = 68 as "a third family that the plan does not license".
What was inconsistent with the manuscript was this log — entry 5 — and the `family` field of
`results/A14/A14b_analysis.json`, both of which still asserted nominal p with no BH. An earlier
draft of this entry claimed the reverse, that the manuscript rested on the withdrawn m = 68
figure; that claim was wrong, is withdrawn, and is recorded here rather than deleted.

**The headline claim is unaffected.** The Framework signature returns residual AUROC 0.548
[0.469, 0.625], LRT p = 0.321, q = 0.538 — null on the second independent cohort, as on the first
four. Condition 1 constrains how the abstract may generalise from that to the field, not the
cell-line-signature result itself.

**Files.** `analysis/A14b_gse41998_analysis.py` is pinned by sha256 in the LOCK's
`deposit_manifest`, so neither it nor its output was edited; a third party re-hashing both still
gets a match. The correction is a new artefact plus this entry. Entry 5 is likewise left standing —
this log is append-only, and rewriting a past entry to agree with a later finding is the failure
mode the append-only rule exists to prevent.

sha256 of this file immediately before entry 11: `58939a534ea6db3a99300250069699c204b95b828274a29a579fe8d63cacdb1f`

### 12 — 2026-08-04 — The LOCK's deposit_manifest no longer verifies: 4 of 41 files mismatch. Observation; no file altered to make it pass.

Re-hashing every path in `ANALYSIS_ADDENDUM_2026-08.LOCK` gives **37/41 match, 4 mismatch**. The
LOCK's own `instructions` field says the manifest exists so "a third party can check that no
pre-lock result was added or altered after the lock", so a silent mismatch defeats its only
purpose. Recording it is the point of this log; it is recorded before any decision about what to
do, so that the decision cannot be made by quietly re-hashing first.

| path | changed by | already disclosed? |
|---|---|---|
| `analysis/A_remaining.py` | `033ca7f` batch6c — re-run of the six analyses the §13.1 resolver fix missed, plus the DCA bootstrap defect | yes, entry 8 |
| `analysis/_fig2_patched.py` | `fd1e1ec` — fig2 stops printing `nan` | as a figure fix; not as a manifest change |
| `analysis/_fig5_patched.py` | `fd1e1ec`, then `102ce33` | as a figure fix; not as a manifest change |
| `results/A0/A0e_dedup.json` | `636ea1d` batch6d — A0e re-run under the corrected resolver | the re-run is disclosed in entries 6 and 8; **this file's manifest entry is not** |

**The last row is the one that matters.** The three scripts are analysis code whose corrections are
described in entries 8 and 9. `A0e_dedup.json` is a **pre-lock result**, and pre-lock results are
exactly what the manifest is there to protect. Its numbers changed for a documented and defensible
reason — entry 6's resolver filter — but a reader verifying the deposit sees a bare hash failure on
a pre-lock artefact with nothing in the manifest to explain it. That is indistinguishable, from
outside, from the thing the manifest is designed to catch.

**Entry 6 promised a regeneration that did not fully happen.** It states that `A0e_dedup.py` "is
modified, so its sha256 changes and `ANALYSIS_ADDENDUM_2026-08.LOCK` is regenerated". The script's
hash does match the manifest, so a regeneration did occur at lock-stamp time; the later re-runs of
08-04 were never folded in. The commitment was kept once and then overtaken.

**No action taken, and deliberately so.** Regenerating the LOCK now would make the manifest verify,
and would also be an author decision about a hash-locked deposit taken while its discrepancy is the
open question — the same shape as re-choosing a threshold after seeing that it fails.
`ANALYSIS_ADDENDUM_2026-08.md` itself still hashes to the locked `6fe26424…`, so the
pre-registration text is intact and nothing here touches it. The options, for the author:
regenerate the LOCK and cite this entry as the reason; or deposit the LOCK unchanged with this
entry as the explanation of the four failures. Either is defensible written down; neither is
defensible silent.

sha256 of this file immediately before entry 12: `037aa405870e744348b7dbdff132d2da8e906fc98f81dbcbd21f861f8b48a981`

### 13 — 2026-08-07 — Two analyses specified and run at revision stage: the endpoint check on GSE41998 (R3-6) and its cohort-characteristics row (R3-3)

**Neither is in the addendum and neither is pre-registered.** Both are labelled
[SPECIFIED AT REVISION STAGE] wherever they appear, the same tag entry 10 introduced for the
gene-count criterion. They are recorded here before their results enter the manuscript.

**A6b — endpoint substitution on GSE41998 (`analysis/A6b_gse41998_endpoint.py`).** The
addendum's A6 harmonises pCR to RCB-0/I on the cohorts that deposit an RCB class; GSE41998 was
not in the study when A6 was written. The manuscript stated that GSE41998 carries an RCB-0/I
call on the same patients and then declined to run the check, calling it "the natural next
step". Reviewer 3's point 6 asked for exactly this and the data were already parsed in the
deposit, so declining cost more than running it.

It is also the cleanest place in the study to ask the question. The other cohorts with an RCB
field, GSE22226 and GSE25066, share 65 I-SPY 1 patients, so an endpoint comparison there is
confounded with the overlap. GSE41998 has zero overlap with any cohort in the study and carries
both calls on the **same** patients, so the endpoint definition is the only thing that varies.

| endpoint | n | events | marginal | residual [95% CI] | LR P |
|---|---:|---:|---:|---|---:|
| pCR (deposited call, as used throughout) | 253 | 69 | 0.500 | 0.548 [0.469, 0.625] | 0.3210 |
| RCB-0/I (harmonised, R3-6) | 253 | 86 | 0.480 | 0.523 [0.448, 0.598] | 0.6690 |

The pCR arm reproduces `A14b_analysis.json` exactly — same n, same events, same residual and
interval — which is the check that the two scripts agree before the new arm is believed.
**Both are null, so the verdict does not depend on how response is defined**, and the sentence
declining the analysis is withdrawn. Under §10 no condition fires: condition 4 requires the
pooled added OR/SD interval to exclude 1.0 or two cohorts to change verdict class, and neither
happens here.

**A5b — the GSE41998 row of the cohort-characteristics table
(`analysis/A5b_gse41998_characteristics.py`).** A5's deposited table covers five cohorts and
not this one, so R3-3 was answered for five sixths of the study. The row is computed from the
phenotype file rather than transcribed from §2.2, so the table and the text have one source;
the two agree (253 analysed, 69 events, 152 ER-negative, 24 HER2-positive, 125 TNBC). HER2
missing is counted negative, following the deposited `A4.strata_of`.

**Neither analysis touches a pre-lock artefact.** Both write new files under `results/`, and
`ANALYSIS_ADDENDUM_2026-08.LOCK` is not regenerated, so the manifest discrepancy recorded as
entry 12 is neither worsened nor quietly repaired.

sha256 of this file immediately before entry 13: `debbbb16e9e0ff84c6b5c2295e2b64ea6ec38bcbea353ca73217bf366331bfaf`
