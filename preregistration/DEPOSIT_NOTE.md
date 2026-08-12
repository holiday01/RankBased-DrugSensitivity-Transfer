# Deposit note — what the hashed files cannot say about themselves

Written 2026-08-12, at the deposit of the pre-registration artefacts.

This file is **not hashed by any lock**. That is deliberate: everything below is a
correction or a clarification that a reader needs, and every file it concerns is either
hash-locked or hashed *by* a lock, so editing any of them to say these things would
break the integrity they exist to provide. The pattern is the one `PREREG_README.md`
already uses one level down — *"Nothing here edits the hashed file"*.

Read `PREREG_README.md` first. It explains the parent pre-registration. This file
explains the deposit.

---

## 1. The DOI. Cite `10.5281/zenodo.20726402`.

**The article, `PREREG_README.md` and `ANALYSIS_ADDENDUM_2026-08.LOCK` all cite
`10.5281/zenodo.20726403`, and all three are wrong in the same way.**

Zenodo issues two DOIs per record. `20726403` is the **version** DOI of v1.0.0 — a
2026-06-17 GitHub release archive of `holiday01/RankBased-DrugSensitivity-Transfer`
containing three notebooks and a `scripts/` tree. It is a real record, but it is not
the deposit this study describes and never was; the curated tree was staged in June and
first published as v2.0.0.

The **concept** DOI is `10.5281/zenodo.20726402`. It resolves to whichever version is
current and is what everything should cite.

| where | what it says | what it should say |
|---|---|---|
| `PREREG_README.md` line 4 | `10.5281/zenodo.20726403 (concept DOI — always resolves to the latest version)` | `10.5281/zenodo.20726402` — and 20726403 is a *version* DOI, not a concept DOI |
| `ANALYSIS_ADDENDUM_2026-08.LOCK` → `parent_concept_doi` | `10.5281/zenodo.20726403` | `10.5281/zenodo.20726402` |
| `ANALYSIS_ADDENDUM_2026-08.LOCK` → `zenodo_doi` | `<FILL AT DEPOSIT>` | this deposit's concept DOI, `10.5281/zenodo.20726402` |

**None of the three was edited.** `PREREG_README.md`'s sha256 is recorded in the
addendum's lock and currently *matches*; editing it to fix a typo would convert a
passing check into a failing one. The lock's own two fields were left alone because
altering a hash-locked artefact after the fact is the manoeuvre the lock exists to make
visible, and it is not worth doing for a field this note can carry instead.

## 2. Which hashes verify, as of this deposit

Recomputed 2026-08-12 against `ANALYSIS_ADDENDUM_2026-08.LOCK` and
`preregistration.LOCK`:

| file | recorded in | result |
|---|---|---|
| `preregistration.md` | `preregistration.LOCK` → `3a9192f3…` | **matches** |
| `ANALYSIS_ADDENDUM_2026-08.md` | addendum lock → `6fe26424…` | **matches** |
| `PREREG_README.md` | addendum lock → `9013ebb8…` | **matches** |
| `decisions_addendum.md` | addendum lock → `897218d4…` | **does not match — by design, see §3** |
| the 41-path `deposit_manifest` | addendum lock | **32 match, 4 differ, 5 are not here — see §4** |

The load-bearing one is the first: the 2026-05-28 lock on the parent pre-registration
verifies. That is the claim the study's confirmatory arm rests on, and it is now
checkable by anyone.

## 3. `decisions_addendum.md` is an append-only log, so its recorded hash is a
   snapshot of an empty file

The addendum's lock records `decisions_addendum_md_sha256` as of 2026-08-02, when the
log had **no entries** — `PREREG_README.md` says so in its own reading order: *"the
deviation log. Empty at lock time."* The log now carries **thirteen** entries, so it
necessarily hashes to something else.

This is the file behaving as specified, not evidence of alteration. The log carries its
own chain instead: several entries record the sha256 of the file immediately before the
entry was appended, so the sequence of states is checkable from inside the document.

## 4. The `deposit_manifest`: 32 of 41 verify here. Two separate reasons.

Re-hashing the manifest's 41 paths **against this deposit** gives **32 match, 4 differ,
5 absent**. Both non-matching groups are deliberate and each has a different cause, so
they are set out separately.

**Note the number you will find quoted elsewhere.** `decisions_addendum.md` entry 12 and
the manuscript's appendix both say *37 of 41*. That is the count in the **authoring
working area**, where all 41 paths exist. In the deposit five of them are absent by the
selection decision below, so the same four differences appear against a smaller
denominator. Both counts are correct about different file sets; neither is a revision of
the other.

### The 4 that differ — entry 12

| path | changed by |
|---|---|
| `analysis/A_remaining.py` | the §13.1 resolver fix and the DCA bootstrap defect, disclosed in entry 8 |
| `analysis/_fig2_patched.py` | a figure fix (fig2 stopped printing `nan`) |
| `analysis/_fig5_patched.py` | two figure fixes |
| `results/A0/A0e_dedup.json` | re-run under the corrected resolver of entry 6 |

**Entry 12 of `decisions_addendum.md` is the full account**, written 2026-08-04 *before*
any decision about what to do — precisely so the decision could not be made by quietly
re-hashing first. The first three are analysis code whose corrections are disclosed in
entries 8 and 9. The fourth, `A0e_dedup.json`, is a **pre-lock result**, the class of
object the manifest exists to protect; its numbers changed for the documented reason in
entry 6, but a third party sees only a bare hash failure. Entry 12 states the two
defensible options — regenerate the lock and cite the entry, or deposit the lock
unchanged with the entry as the explanation.

**This deposit takes the second option.** The lock is deposited exactly as written on
2026-08-02. Nothing was re-hashed to make a check pass.

### The 5 that are absent — a selection decision, not a hash problem

`results/A0/score_deployed_GSE16446.tsv`, `…GSE22226.tsv`, `…GSE25066.tsv`,
`…GSE32646.tsv` and `results/A0/A0b_rescore_comparison.json`.

These are the intermediate records of correcting the scorer defect. No reviewer asked
for them, no reported result reads them, and the scripts that produced them **are**
deposited, so the correction stays auditable. They were dropped on 2026-08-11 under the
rule that the deposit carries what the reviewers' analyses need; the selection was
computed from the path literals in `analysis/*.py` and a grep of the four submitted
documents, not chosen case by case.

## 5. What is here, and what each thing is

| file | what it is | locked |
|---|---|---|
| `preregistration.md` | the parent pre-registration: the confirmatory family T1–T6, the estimator, and the success threshold. All six tests were executed and **all six failed**. | yes, 2026-05-28 |
| `preregistration.LOCK` | its sha256 and lock timestamp | — |
| `PLAN_v2.md` | post-hoc planning. **Not a pre-registration**; its own §2.2 says its §4 analyses are post-hoc relative to `preregistration.LOCK`. The residual-AUROC estimator and the 0.55 margin trace to this document, not to the locked one. | **no** |
| `ANALYSIS_ADDENDUM_2026-08.md` | this revision's analysis plan, written before the revision analyses ran | yes, 2026-08-02 |
| `ANALYSIS_ADDENDUM_2026-08.LOCK` | its sha256, the deposit manifest, and the two DOI fields §1 corrects | — |
| `decisions_addendum.md` | the deviation log, thirteen entries, append-only | no, by design |
| `PREREG_README.md` | four presentation defects in the parent pre-registration, disclosed by the authors rather than left to be found | — |
| `DEPOSIT_NOTE.md` | this file | no |

**The revision plan was not publicly deposited before its analyses ran.** It carries a
written sha256 from before execution, which fixes the analysis list and the decision
rules in advance; it does not carry third-party attestation of that ordering, and the
manuscript says so. No revision analysis is described as pre-registered.
