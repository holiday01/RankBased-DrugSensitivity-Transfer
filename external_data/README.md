# external_data — cohort data, fetched by you

The analyses resolve their inputs through the `DATA_ROOT` environment variable, and
this directory is where the cohort data goes. No cohort data is distributed here;
every cohort is public and belongs to the group that published it. This file lists
what each script needs and where it comes from.

## Downloaded by a script in this deposit

| cohort | script |
|---|---|
| GSE16446, GSE25066, GSE22226 | `code/framework/download_gse.py` (series matrix + platform annotation, each logged with its sha256) |
| GSE32646 | `code/axis_analysis/phase_2/gse32646_download.py` |
| METABRIC | `code/axis_analysis/phase_4/metabric_locate_download.py` (cBioPortal datahub, `brca_metabric`) |
| GSE9782, GSE68871, GSE136337 | `code/axis_analysis/phase_D/download_*.py` |

Parsing to a matrix and extracting pCR labels: `code/framework/parse_gse_to_matrix.py`
and `code/framework/extract_pcr_labels.py`.

## GSE41998

No download script is provided. `analysis/A14a_gse41998_gate.py` builds
`GSE41998_expression.tsv`, `GSE41998_axes.tsv` and `GSE41998_phenotype.tsv` from
per-sample GEO tables it reads from:

```
external_data/GSE41998/samples/GSM*_sample_table.txt
```

`GSE41998_sample_accessions.json` in this directory lists the 279 accessions the
analysis used; each is retrievable from GEO. The script also reads
`external_data/GSE25066/platform_annot.annot.gz` (GSE41998 is on GPL96, the same
platform, and `download_gse.py` fetches that annotation) and
`external_data/pam50_centroids.tsv`.

## Reference files

`pam50_centroids.tsv` (published PAM50 centroids), `hgnc_alias_table.tsv` and
`hgnc_complete_set.tsv` (HGNC symbol tables), `oncokb_all_universe.tsv` (the OncoKB
gene universe). All public.
