# CLAUDE.md — oncoanalyser

Oncoanalyser WGS/WTS + clinical CSVs → cBioPortal. HPC-only (SLURM + Apptainer).

## Key Rules

- `--type both` (default) or `--type dna-only` — dna-only skips all RNA-derived processing (Isofox expression, Isofox fusion, SAGE RNA-append integration into mutations) and omits the expression profile/meta file, even if RNA files exist on disk. Validated in `subworkflows/local/utils/main.nf`. Mutational signature fitting (Sigs) is unaffected either way — it's DNA-only already, but the VCF it reads changes (see below).
- `--type dna-only` also changes the somatic/germline VCF source: oncoanalyser only runs PAVE (annotation) and `sage_append` (RNA-count append) when RNA is present, so dna-only reads the raw SAGE output — `sage/somatic/<subject>-T.sage.somatic.vcf.gz` and `sage/germline/<subject>-T.sage.germline.vcf.gz` — instead of `pave/<subject>-T.pave.somatic.vcf.gz` / `pave/<subject>-T.pave.germline.vcf.gz`. Applies everywhere a somatic VCF is read: mutations (`ch_sage_vcf`), DBS/ID signature fitting (`ch_sigs_dbs`, `ch_sigs_id`) in `workflows/genomic.nf`.
- `--type dna-only` skips steps per oncoanalyser tool: a subject whose `<folder>/<tool>/` directory is absent (`sage/somatic`, `sage/germline`, `purple`, `esvee`, `sigs`) is dropped from that tool's channels (`hasToolOutput()` in `workflows/genomic.nf`, one summary warning per tool, not per subject). In `GENOMIC_AGGREGATE_OUTPUT`, meta files and case lists are emitted only for data types that produced output (e.g. WES without `sigs/` → no SBS data or SBS meta files). The case list label travels with its sample list, because `GENERATE_CASE_LIST`'s three inputs are otherwise paired by position. "both" mode keeps the fixed upstream meta and case lists
- `PACKAGE_CBIOPORTAL` must receive every `GENOMIC_AGGREGATE_OUTPUT` output (including `sigs_counts_id`); a missing `.mix()` leaves the file in the output folder but not in the `.tar.gz`
- `cancer_type` must be an OncoTree code (https://oncotree.mskcc.org/) already in the target portal, unless `generate_cancer_type = true`
- R scripts → `container_r`; Python → `container_python`; SigProfiler → `container_sigprofiler`
- `data_sv.txt` rows require Hugo symbols at both sites — filter unannotated rows
- SV classification (`gen_esvee_sv_to_cbioportal.R`): BND ALT strand → `(+,-)` DEL, `(-,+)` DUP, `(+,+)/(−,−)` INV, diff chr TRANSLOC. DNA SVs: `DNA_Support=Yes, RNA_Support=No`
- RNA fusions (`gen_isofox_fusion_to_cbioportal.R`): `Class=FUSION, DNA_Support=No, RNA_Support=Yes`. Both merge into `data_sv.txt`
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`
- No internet on compute nodes — `NXF_OFFLINE=true`; pre-pull containers on login nodes
- `container_vcf2maf` must bundle **both** `vcf2maf.pl` and `vep` (e.g. `vcf2maf_ensembl-vep`) whenever `vep_data` is set — a vcf2maf-only image makes VCF2MAF fail
- Keep `executor.queueSize` well under the cluster's per-account `MaxSubmitJobs` (1000 on Narval/Rorqual) or large cohorts die on `AssocMaxSubmitJobLimit`
- VEP/PCGR data must be pre-staged
- SAGE/PAVE germline VCFs carry no genotypes (`GT=./.`), and CPSR's `vcf2tsvpy` step drops `./.` rows (there is no CPSR flag to change this), so CPSR silently classifies nothing. `FILTER_GERMLINE_DNA` therefore sets GT from the normal's `FMT/AF` (`≥0.9` → `1/1`, `>0` → `0/1`) with `bcftools +setGT`, and removes `INFO/IMPACT` (PAVE), which clashes with a CPSR 2.3.1 reserved tag and aborts CPSR
- `Mutation_Status` in the per-sample MAF: rows from the somatic vcf2maf MAF → `Somatic` (set in `gen_convert_cpsr_to_maf.R` where blank; in both mode the `RNA(d=,v=,filt=)` text from INTEGRATE_RNA_VARIANTS is kept), rows from CPSR → `Germline`. Germline has priority: a somatic row that is the same variant as a CPSR call (chrom ignoring `chr`, Start_Position, Reference_Allele, Tumor_Seq_Allele2) is dropped, so only the CPSR row remains. CPSR rows get the `chr` prefix like vcf2maf rows. Never leave it blank: `ml_format_mutation.R` splits on it, and dplyr `filter()` drops NA

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

The genomic workflow checks for pre-existing output files per subject. If all expected outputs (CNV, SV, expression, mutations) already exist, processing is skipped. Both the skip decision and the reuse of cached files go through `isSubjectComplete()` in `workflows/genomic.nf` — a subject that is re-run must never also be read from cache, or every file reaches the merge steps twice. In `--type dna-only` no `tpm.tsv` is produced, so every subject is re-run (rely on `-resume`). This allows adding new subjects to the samplesheet and re-running without reprocessing the entire cohort. Delete a subject's output directory to force reprocessing.

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
