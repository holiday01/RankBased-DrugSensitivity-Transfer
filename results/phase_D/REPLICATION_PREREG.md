# Track D — Bortezomib signature replication: PRE-REGISTERED cohort list & plan

**Locked 2026-06-11, BEFORE touching any GSE19784/GSE2658 expression data.**
User-directed: keep GSE9782, add two external cohorts. To avoid cohort-shopping,
the cohort list is frozen here and ALL results are reported regardless of outcome.

## Signatures (frozen; sha256 in build logs)
- `BORTEZOMIB.tsv` — pan-cancer pool (framework verbatim)
- `BORTEZOMIB_heme.tsv` — heme-lineage-matched

## Cohorts (FROZEN — no additions to rescue a negative)
| Cohort | Setting | Drug | Endpoint | Status |
|---|---|---|---|---|
| GSE9782 | relapsed/refractory | bortezomib MONO | binary response (R/NR) | tested POSITIVE (AUROC 0.60–0.66) — retained |
| GSE68871 | frontline | VTD combo | response (≥VGPR etc.) | tested NULL (~0.48) — retained in record |
| **GSE19784** | frontline, transplant-elig | **PAD (bortezomib) vs VAD (control)** | PFS / OS | NEW — to test |
| **GSE2658** | frontline | UAMS TT (bortezomib in TT3) | OS / EFS | NEW — to test |

## Pre-specified analysis (same singscore deployment as GSE9782)
1. Score both signatures via full-transcriptome singscore on each cohort.
2. **Survival association (primary for the new cohorts):** Cox PH, HR per +1 SD of
   signature score, on OS and PFS/EFS. Direction expectation: higher score =
   more bortezomib-sensitive ⇒ HR < 1 (better outcome).
3. **Binary surrogate (for AUROC comparability with GSE9782):** event
   (progression/death) by cohort-median follow-up ⇒ AUROC + bootstrap CI.
4. **GSE19784 specificity test (the key one):** fit signature×arm interaction.
   A bortezomib-PREDICTIVE (not merely prognostic) signature should associate
   with outcome MORE strongly in the PAD (bortezomib) arm than the VAD arm.
5. Baselines: MKI67 (proliferation), UAMS-70 / proliferation index if available.

## Reporting commitment
- Report HR, CI, p, C-index, AUROC for BOTH signatures on BOTH new cohorts and
  BOTH arms of GSE19784, regardless of result.
- The verdict integrates ALL 4 cohorts (GSE9782 + GSE68871 + GSE19784 + GSE2658).
- No additional cohorts will be added after seeing these results.

## Setting caveat (honest)
GSE19784/GSE2658 are frontline, bortezomib-in-COMBINATION, survival endpoints —
NOT matched to GSE9782 (relapsed, mono, response). They test cross-setting
generalization, not strict replication. Interpreted accordingly.
