# CLAUDE.md — oncoanalyser

Oncoanalyser WGS/WTS + clinical CSVs → cBioPortal. HPC-only (SLURM + Apptainer).

## Key Rules

- `--type both` (default) or `--type dna-only` — dna-only skips all RNA-derived processing (Isofox expression, Isofox fusion, SAGE RNA-append integration into mutations) and omits the expression profile/meta file, even if RNA files exist on disk. Validated in `subworkflows/local/utils/main.nf`. Mutational signature fitting (Sigs) is unaffected either way — it's DNA-only already, but the VCF it reads changes (see below).
- `--type dna-only` also changes the somatic/germline VCF source: oncoanalyser only runs PAVE (annotation) and `sage_append` (RNA-count append) when RNA is present, so dna-only reads the raw SAGE output — `sage/somatic/<subject>-T.sage.somatic.vcf.gz` and `sage/germline/<subject>-T.sage.germline.vcf.gz` — instead of `pave/<subject>-T.pave.somatic.vcf.gz` / `pave/<subject>-T.pave.germline.vcf.gz`. Applies everywhere a somatic VCF is read: mutations (`ch_sage_vcf`), DBS/ID signature fitting (`ch_sigs_dbs`, `ch_sigs_id`) in `workflows/genomic.nf`.
- R scripts → `container_r`; Python → `container_python`; SigProfiler → `container_sigprofiler`
- `data_sv.txt` rows require Hugo symbols at both sites — filter unannotated rows
- SV classification (`gen_esvee_sv_to_cbioportal.R`): BND ALT strand → `(+,-)` DEL, `(-,+)` DUP, `(+,+)/(−,−)` INV, diff chr TRANSLOC. DNA SVs: `DNA_Support=Yes, RNA_Support=No`
- RNA fusions (`gen_isofox_fusion_to_cbioportal.R`): `Class=FUSION, DNA_Support=No, RNA_Support=Yes`. Both merge into `data_sv.txt`
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`
- No internet on compute nodes — `NXF_OFFLINE=true`; pre-pull containers on login nodes
- VEP/PCGR data must be pre-staged

## Timeline Generation (`gen_timeline.R`)

Generates a single combined `data_timeline.txt` from ARGO clinical CSVs. All event types are merged into one file, distinguished by EVENT_TYPE column:

- `SURGERY` — surgical treatments (subtype, site, intent)
- `TREATMENT` — systemic therapies (drug, type, intent)
- `STATUS` — follow-up disease status
- `SPECIMEN` — specimen collection (site, type, sample type)
- `LAB_TEST` — biomarker results (PSA, ER/PR/HER2, CEA, etc.)

Columns are the union of all event types; columns not applicable to a given EVENT_TYPE are empty. Common columns first (PATIENT_ID, START_DATE, STOP_DATE, EVENT_TYPE), then the rest alphabetically.

Key behaviors:

- **START_DATE defaults to 0 when missing** — any NA/empty date becomes day 0 (diagnosis day)
- Filters patients from `sample_registrations` (Tumour + Total DNA + Solid tissue + regex `-\d+[A-Z]*[DR]T$`)
- When `genomic_subjects` TSV is provided (both mode), restricts timeline output to genomic subjects only
- Nearest follow-up visit determines specimen SAMPLE_TYPE (Primary / Recurrence / Metastasis)

## MOHCCN Mapping Tables

Three mapping CSVs in `assets/` translate MOHCCN plain-text values to ontology/ICD-O codes:

- `mohccn_*_primary_site.csv` — ICD-O topography codes (reversed: code → label for surgery/specimen site)
- `mohccn_*_specimen_tissue_source.csv` — tissue source codes
- `mohccn_*_treatment_intent.csv` — treatment intent codes

Params: `mohccn_primary_site_map`, `mohccn_specimen_tissue_source_map`, `mohccn_treatment_intent_map`. Used by both `clin_format.R` and `gen_timeline.R`.

## Incremental Processing

The genomic workflow checks for pre-existing output files per subject. If all expected outputs (CNV, SV, expression, mutations) already exist, processing is skipped. This allows adding new subjects to the samplesheet and re-running without reprocessing the entire cohort. Delete a subject's output directory to force reprocessing.

## Process Labels (`conf/base.config`)

Default: 1 CPU, 1 GB, 4 min (scaled by `task.attempt`). See `conf/base.config` for full table.

## Mutational Signatures

See `docs/mutational_signatures.md` for detailed SBS/DBS/ID documentation.

## Testing

nf-test gotchas:

- Use `path(f.toString())` — channel file outputs are `String`, not `Path`
- Sort snapshots: `.sort { it.toString().split('/').last() }`
- `collectFile` with `storeDir` won't create dirs — call `file("${params.outdir}/GROUP").mkdirs()` in test setup
- `genomic_ml` uses `options "-stub-run"` to skip `DOWNLOAD_KNOWN_FUSIONS`
