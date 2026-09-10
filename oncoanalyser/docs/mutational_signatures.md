# Mutational Signatures — oncoanalyser

Three signature types (SBS, DBS, ID) are fitted against COSMIC v3.6 GRCh38 reference (`assets/COSMIC_Human_v3.6.zip`), producing 6 cBioPortal GENERIC_ASSAY data files per group. Metadata files in `assets/` provide `category`, `etiology`, `main_effect` for NAME/DESCRIPTION columns.

## SBS — Single Base Substitutions (101 signatures)

**Container:** `container_sigprofiler` (SigProfilerAssignment + pandas + scipy + pysam)

**Contribution** (`data_mutational_signatures_contribution_SBS.txt`):

- Input: `{subject}-T.sig.snv_counts.csv` (96-channel trinucleotide counts)
- Script: `bin/run_sigprofiler_sbs.py` → `SigProfilerAssignment.Analyzer.cosmic_fit(context_type="96")`
- Metadata: `assets/cosmic_sbs_metadata.tsv`
- Columns: `ENTITY_STABLE_ID` (`mutational_signatures_contribution_{SBS_ID}`), `NAME`, `DESCRIPTION`, `{sample}`
- Merge: `gen_merge_sigs_to_cbioportal.R` (R, outer join by ENTITY_STABLE_ID)

**Counts** (`data_mutational_signatures_counts_SBS.txt`):

- Input: same `sig.snv_counts.csv`
- Script: `bin/gen_sigs_counts_to_cbioportal.R` (context transform: `C>A_ACA` → `A[C>A]A`)
- Columns: `ENTITY_STABLE_ID` (`mutational_signatures_matrix_{5'}_{from}-{to}_{3'}`), `NAME`, `{sample}`
- Merge: `gen_merge_sigs_counts_to_cbioportal.R`

## DBS — Doublet Base Substitutions (22 signatures)

**Contribution** (`data_mutational_signatures_contribution_DBS.txt`):

- Input: `{subject}-T.pave.somatic.vcf.gz` (`{subject}-T.sage.somatic.vcf.gz` under `--type dna-only`, since PAVE/sage_append don't run without RNA) — extracts adjacent SNV pairs + 2bp MNVs → 78-channel DBS matrix
- Script: `bin/run_sigprofiler_dbs.py` → `Analyzer.cosmic_fit(context_type="DINUC")`
- Metadata: `assets/cosmic_dbs_metadata.tsv`
- Strand normalization to 10 canonical ref dinucleotides (AC, AT, CC, CG, CT, GC, TA, TC, TG, TT)

**Counts** (`data_mutational_signatures_counts_DBS.txt`):

- Same script produces both contribution and counts
- `ENTITY_STABLE_ID` = `mutational_signatures_matrix_{ref}-{alt}`, `NAME` = `AC>CA` format

## ID — Indels (17 signatures)

**Contribution** (`data_mutational_signatures_contribution_ID.txt`):

- Input: `{subject}-T.pave.somatic.vcf.gz` (`{subject}-T.sage.somatic.vcf.gz` under `--type dna-only`) — extracts indels → 83-channel ID matrix via pysam + reference FASTA
- Script: `bin/run_sigprofiler_id.py` → `Analyzer.cosmic_fit(context_type="ID")`
- Metadata: `assets/cosmic_id_metadata.tsv`
- 83-channel classification: 1bp C/T del/ins at homopolymers (24) + 2-5bp repeat-mediated del/ins (48) + 2-5bp microhomology del (11)

**Counts** (`data_mutational_signatures_counts_ID.txt`):

- Same script produces both contribution and counts
- `ENTITY_STABLE_ID` = `mutational_signatures_matrix_{size}_{type}_{base}_{count}`, `NAME` = `1:Del:C:0` format

## Common

- All meta files use `generic_assay_type: MUTATIONAL_SIGNATURE`, `datatype: LIMIT-VALUE`, `pivot_threshold_value: 0.0`.
- Merge modules reuse `gen_merge_sigs_to_cbioportal.R` / `gen_merge_sigs_counts_to_cbioportal.R` with different output filenames.
- Zero-mutation edge case: scripts write header-only files; <50 mutations: stderr warning but continues.
