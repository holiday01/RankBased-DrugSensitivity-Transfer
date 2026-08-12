# Reading the pre-registration deposit

This note accompanies `preregistration.md` and `preregistration.LOCK` in Zenodo record
**10.5281/zenodo.20726403** (concept DOI — always resolves to the latest version). It exists because three features of the deposited pre-registration look
like defects on a first read, and the authors would rather explain them than have a reader
discover them unexplained. **Nothing here edits the hashed file**, which must not change.

Committed in `ANALYSIS_ADDENDUM_2026-08.md` §3.8.

---

## 1. The cryptographic lock is sound

```
sha256(preregistration.md) = 3a9192f3adf0da350d34cf978233799ba571d7e74a5392a2b56201be20f485f7
preregistration.LOCK       = the same value, plus manifest_json_sha256 and
                             lock_timestamp_utc 2026-05-28T02:14:24+00:00
```

All four deposited copies of `preregistration.md` are byte-identical. Re-hash the file and
compare against the LOCK file to verify.

## 2. Line 5 still reads `Status: DRAFT`

It does, and it should not. The header was written before the lock and was never updated,
while §9 and §10b of the same document both assert "locked 2026-05-28". The document was
locked; the header simply was not revised before hashing, and revising it afterwards would
have broken the hash the lock depends on.

Filesystem times are consistent with a genuine lock: `preregistration.md` last modified
02:09, `preregistration.LOCK` written 02:14, same day, and both files are mode `r--r--r--`.

## 3. §11.1 "Self-lock" contains only forward-references

§11.1 says the lock timestamp and self-sha256 are "captured in companion file
`preregistration.LOCK`" rather than containing them inline, while §11 line 469 makes filling
§11.1 a *precondition* of being locked. Read literally, the document declares itself unlocked
while being demonstrably locked by its hash.

This is a drafting circularity, not an integrity problem: a document cannot contain its own
sha256, which is why the value lives in the companion file. The companion file exists,
verifies, and carries the timestamp.

## 4. §9 cites stale hashes for two scripts

§9 records sha256 values for `freeze_signatures.py` and `verify_frozen.py` that do not match
the deposited files:

| Script | §9 records | Actual and MANIFEST-recorded |
|---|---|---|
| `freeze_signatures.py` | `cc6cde39b4976fd1…` | `d1619e203b341b18…` |
| `verify_frozen.py` | `c4c11b7a2224de91…` | `475e2d8958e64e13…` |

**Every signature-artefact hash in §9 matches.** Only these two script hashes differ,
consistent with the scripts being extended for the Round 2/3 additions after §9 was first
filled — the script mtime is 2026-05-28 01:47, the MANIFEST timestamp 02:14 the same day.
The file that hashed to `cc6cde39…` is not on disk and the directory is not a git repository,
so that intermediate state cannot be recovered.

**Why this needs saying.** §9 line 417 states: *"Any artifact whose full sha256 differs from
`MANIFEST.json` invalidates the pre-registration."* Read at its widest, that clause could be
invoked here. It was written to protect the **frozen signature artefacts** — the objects
whose stability the pre-registration's claims depend on — and all of those verify. The
authors flag the tension rather than rely on a reader not noticing it.

## 5. `PLAN_v2.md` is in this record and is NOT locked

The record also contains `PLAN_v2.md`. It is not part of the hash-locked pre-registration and
**no `PLAN_v2.LOCK` exists**, although PLAN_v2's own §2.3 specifies the locking procedure it
would need (sha256 hash, Zenodo and OSF deposit, a DOI and "NOT a [TBC] placeholder"). Its
own §2.2 states that all its §4 analyses are "POST-HOC relative to preregistration.LOCK".

This matters for reading the manuscript. Several manuscript passages describe a
"pre-specified" residual-AUROC estimator and a "pre-specified margin" of 0.55. Those terms
trace to PLAN_v2, not to the hash-locked pre-registration, in which `residual`, `TOST` and
`equivalence` occur **zero times**. PLAN_v2 states the 0.55 threshold in two senses — a
falsifier CI-lower bound (L350) and a null-confirmation ceiling (§4.8.3 L636, D19 L740) — and
attaches Bonferroni α = 0.01 to the endpoint (L353, L513, L745).

`ANALYSIS_ADDENDUM_2026-08.md` §3.4 and §5.1 correct the manuscript on all of this and set
out the full provenance. In short: the estimators are **post-hoc**, and the addendum says so
in the manuscript rather than leaving the reader to reconstruct it from the deposit.

---

## What to read, in order

1. `preregistration.md` + `.LOCK` — the confirmatory family T1–T6. All six were executed and
   all six failed; that is reported in the manuscript.
2. `PLAN_v2.md` — post-hoc planning, unlocked. Not a pre-registration.
3. `ANALYSIS_ADDENDUM_2026-08.md` + `.LOCK` — the revision's pre-specification, together with
   eight author-found corrections to the record, filed before any revision analysis was run.
4. `decisions_addendum.md` — the deviation log. Empty at lock time.
