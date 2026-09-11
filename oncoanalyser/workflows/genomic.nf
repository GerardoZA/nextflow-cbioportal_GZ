/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { GENOMIC_CNV                  } from '../subworkflows/local/genomic_cnv'
include { GENOMIC_SV                   } from '../subworkflows/local/genomic_sv'
include { GENOMIC_EXPRESSION           } from '../subworkflows/local/genomic_expression'
include { GENOMIC_MUTATIONS            } from '../subworkflows/local/genomic_mutations'
include { GENOMIC_ML                   } from '../subworkflows/local/genomic_ml'
include { GENOMIC_AGGREGATE_OUTPUT     } from '../subworkflows/local/genomic_aggregate_output'
include { GENERATE_META_FILE           } from '../modules/local/generate_meta_file'
include { GENERATE_CANCER_TYPE        } from '../modules/local/generate_cancer_type'
include { ISOFOX_FUSION_TO_CBIOPORTAL  } from '../modules/local/isofox_fusion_to_cbioportal'
include { MERGE_SAMPLE_SV              } from '../modules/local/merge_sample_sv'
include { SIGPROFILER_SBS               } from '../modules/local/sigprofiler_sbs'
include { SIGPROFILER_DBS              } from '../modules/local/sigprofiler_dbs'
include { SIGPROFILER_ID               } from '../modules/local/sigprofiler_id'
include { SIGS_COUNTS_TO_CBIOPORTAL    } from '../modules/local/sigs_counts_to_cbioportal'
include { PACKAGE_CBIOPORTAL           } from '../modules/local/package_cbioportal'


// True when running without RNA data. Compares case-insensitively because params values
// given on the command line are immutable — "--type DNA-only" cannot be normalised in
// place by the validator, so a plain `params.type == 'dna-only'` would silently be false.
def isDnaOnly() {
    return (params.type ?: 'both').toString().toLowerCase() == 'dna-only'
}

// Resolve an oncoanalyser output file path; log a warning and return null if absent.
def findOncoFile(meta, path_str, label) {
    def f = file(path_str, checkIfExists: false)
    if (!f.exists() || f.isEmpty()) {
        log.warn "File not found for ${meta.subject} (${label}): ${path_str}"
        return null
    }
    return f
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMIC {

    take:
        samplesheet_list        // from nf-schema — one record per samplesheet row
        ensembl_annotations     // path — BioMart TSV for CNV gene mapping + ESVEE gene overlap
        ensembl_annotations_expr// path — BioMart TSV for Isofox expression Entrez ID mapping
        vep_data                // channel<path> — pre-staged VEP cache (may be empty)
        pcgr_data               // channel<path> — pre-staged PCGR reference data (may be empty)
        needs_vep               // boolean — true when vep_data is not supplied
        needs_pcgr              // boolean — true when pcgr_data is not supplied
        fasta                   // path — GRCh38 reference FASTA (for vcf2maf)
        cosmic_data             // channel<path> — COSMIC/ChimerKB fusion data for ML step
        chimer_data

    main:

        ch_versions = channel.empty()

        // ── Parse samplesheet ────────────────────────────────────────────────
        // One row per subject: group, subject_id, sample_id, folder.
        // All modality files are resolved relative to folder.

        ch_samples = samplesheet_list
            .map { rec ->
                [
                    group  : "${params.study_id}",
                    subject: "${rec[0].subject}",
                    sample : "${rec[0].sample}",
                    folder : "${rec[0].folder}",
                ]
            }

        // ── Incremental processing: detect already-processed subjects ────────
        // Check output directory for existing per-subject files.
        // If all expected outputs exist, skip processing and reuse cached files.

        ch_files_ran = ch_samples
            .filter { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                baseDir.exists() && baseDir.isDirectory()
            }
            .flatMap { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                def files = []

                // CNV files (use meta.sample for filename)
                def seg = file("${baseDir}/${meta.sample}_data_cna_hg38.seg", checkIfExists: false)
                def long_cnv = file("${baseDir}/${meta.sample}_data_cna_long.txt", checkIfExists: false)
                if (seg.exists() && long_cnv.exists()) {
                    files.add([meta + [pipeline: 'cnv'], [seg: seg, longfile: long_cnv]])
                }

                // SV file
                def sv = file("${baseDir}/${meta.sample}.data_sv.txt", checkIfExists: false)
                if (sv.exists()) {
                    files.add([meta + [pipeline: 'sv'], sv])
                }

                // Expression file
                def tpm = file("${baseDir}/${meta.sample}.tpm.tsv", checkIfExists: false)
                if (tpm.exists()) {
                    files.add([meta + [pipeline: 'expression'], tpm])
                }

                // Mutation file
                def maf = file("${baseDir}/${meta.sample}.somatic_rna_germline.maf", checkIfExists: false)
                if (maf.exists()) {
                    files.add([meta + [pipeline: 'mutation'], maf])
                }

                // RNA fusion file (optional — only present for samples with RNA data)
                def rna_fusion = file("${baseDir}/${meta.sample}.isofox_fusion.data_sv.txt", checkIfExists: false)
                if (rna_fusion.exists()) {
                    files.add([meta + [pipeline: 'sv_rna_fusion'], rna_fusion])
                }

                // Signatures contribution file (optional — only present when sufficient mutations exist)
                def sigs = file("${baseDir}/${meta.sample}.data_mutational_signatures_contribution_SBS.txt", checkIfExists: false)
                if (sigs.exists()) {
                    files.add([meta + [pipeline: 'sigs'], sigs])
                }

                // Signatures counts file (optional — only present when sufficient mutations exist)
                def sigs_counts = file("${baseDir}/${meta.sample}.data_mutational_signatures_counts_SBS.txt", checkIfExists: false)
                if (sigs_counts.exists()) {
                    files.add([meta + [pipeline: 'sigs_counts'], sigs_counts])
                }

                // DBS signature contribution file
                def sigs_dbs = file("${baseDir}/${meta.sample}.data_mutational_signatures_contribution_DBS.txt", checkIfExists: false)
                if (sigs_dbs.exists()) {
                    files.add([meta + [pipeline: 'sigs_dbs'], sigs_dbs])
                }

                // DBS signature counts file
                def sigs_counts_dbs = file("${baseDir}/${meta.sample}.data_mutational_signatures_counts_DBS.txt", checkIfExists: false)
                if (sigs_counts_dbs.exists()) {
                    files.add([meta + [pipeline: 'sigs_counts_dbs'], sigs_counts_dbs])
                }


                return files
            }

        // Get set of already-processed subject names (those with all 4 outputs)
        existing_subjects = ch_samples
            .filter { meta ->
                def baseDir = file("${params.outdir}/${meta.group}/${meta.subject}")
                if (!baseDir.exists()) return false
                def seg = file("${baseDir}/${meta.sample}_data_cna_hg38.seg", checkIfExists: false)
                def long_cnv = file("${baseDir}/${meta.sample}_data_cna_long.txt", checkIfExists: false)
                def sv = file("${baseDir}/${meta.sample}.data_sv.txt", checkIfExists: false)
                def tpm = file("${baseDir}/${meta.sample}.tpm.tsv", checkIfExists: false)
                def maf = file("${baseDir}/${meta.sample}.somatic_rna_germline.maf", checkIfExists: false)
                return seg.exists() && long_cnv.exists() && sv.exists() && tpm.exists() && maf.exists()
            }
            .map { meta -> meta.subject }
            .collect()
            .map { it.toSet() }
            .ifEmpty([] as Set)

        // Filter to only new samples that need processing
        ch_samples_to_run = ch_samples
            .combine(existing_subjects)
            .filter { meta, existing_set -> meta.subject !in existing_set }
            .map { meta, existing_set -> meta }

        // Log skipped samples
        ch_samples
            .combine(existing_subjects)
            .filter { meta, existing_set -> meta.subject in existing_set }
            .subscribe { meta, set ->
                log.info "Skipping already-processed subject: ${meta.subject}"
            }

        // Warn if all samples already processed
        ch_samples_to_run
            .collect()
            .filter { list ->
                if (list.isEmpty()) {
                    log.warn "All subjects in samplesheet already processed. Skipping genomic processing."
                }
                return true
            }

        // ── Build per-modality input channels ─────────────────────────────────
        // File naming conventions (relative to folder):
        //   sage/somatic/${sample_id}-T.sage.somatic.vcf.gz   — SAGE somatic DNA
        //   sage/germline/${sample_id}-N.sage.germline.vcf.gz — SAGE germline DNA
        //   sage/append/${sample_id}-T.sage.append.vcf.gz     — SAGE somatic RNA append
        //   esvee/caller/${sample_id}-T.esvee.unfiltered.vcf.gz
        //   purple/${sample_id}-T.purple.cnv.somatic.tsv
        //   purple/${sample_id}-T.purple.cnv.gene.tsv
        //   isofox/${sample_id}-T.isf.gene_data.csv

        // SAGE somatic VCF → mutations
        // dna-only: oncoanalyser never runs PAVE/sage_append without RNA, so the somatic
        // VCF is the raw SAGE output (sage/somatic/) rather than the PAVE-annotated one
        ch_sage_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = isDnaOnly()
                    ? findOncoFile(meta,
                        "${meta.folder}/sage/somatic/${meta.subject}-T.sage.somatic.vcf.gz",
                        'mutation (SAGE somatic, dna-only)')
                    : findOncoFile(meta,
                        "${meta.folder}/pave/${meta.subject}-T.pave.somatic.vcf.gz",
                        'mutation (PAVE somatic)')
                vcf ? [meta + [pipeline: 'mutation'], vcf] : null
            }
            .filter { it != null }

        // SAGE germline VCF → germline mutations (same dna-only caveat as above)
        ch_sage_germline_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = isDnaOnly()
                    ? findOncoFile(meta,
                        "${meta.folder}/sage/germline/${meta.subject}-T.sage.germline.vcf.gz",
                        'germline mutation (SAGE, dna-only)')
                    : findOncoFile(meta,
                        "${meta.folder}/pave/${meta.subject}-T.pave.germline.vcf.gz",
                        'germline mutation (PAVE)')
                vcf ? [meta + [pipeline: 'mutation_germline'], vcf] : null
            }
            .filter { it != null }

        // SAGE RNA-append VCF → somatic RNA mutations (skipped entirely in dna-only mode)
        ch_sage_rna_vcf = isDnaOnly()
            ? Channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def vcf = findOncoFile(meta,
                        "${meta.folder}/sage_append/somatic/${meta.subject}-T.sage.append.vcf.gz",
                        'mutation (SAGE RNA append)')
                    vcf ? [meta + [pipeline: 'mutation_rna'], vcf] : null
                }
                .filter { it != null }

        // PURPLE CNV somatic + gene TSV → copy-number
        ch_purple_cnv = ch_samples_to_run
            .map { meta ->
                def somatic = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.somatic.tsv",
                    'cnv (PURPLE somatic)')
                def gene = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.gene.tsv",
                    'cnv (PURPLE gene)')
                (somatic && gene) ? [meta + [pipeline: 'cnv'], somatic, gene] : null
            }
            .filter { it != null }

        // ESVEE unfiltered VCF (tumor only) → structural variants
        ch_esvee_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/esvee/${meta.subject}-T.esvee.somatic.vcf.gz",
                    'sv (ESVEE tumor)')
                vcf ? [meta + [pipeline: 'sv'], vcf] : null
            }
            .filter { it != null }

        // Isofox gene expression CSV → TPM (skipped entirely in dna-only mode)
        ch_isofox_exp = isDnaOnly()
            ? Channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def exp = findOncoFile(meta,
                        "${meta.folder}/isofox/${meta.subject}-T-RNA.isf.gene_data.csv",
                        'expression (Isofox)')
                    exp ? [meta + [pipeline: 'expression'], exp] : null
                }
                .filter { it != null }

        // Isofox pass_fusions CSV (tumor RNA) → RNA fusions for data_sv.txt (skipped entirely in dna-only mode)
        ch_isofox_fusion = isDnaOnly()
            ? Channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def fusions = findOncoFile(meta,
                        "${meta.folder}/isofox/${meta.subject}-T-RNA.isf.pass_fusions.csv",
                        'rna fusion (Isofox)')
                    fusions ? [meta + [pipeline: 'sv_rna_fusion'], fusions] : null
                }
                .filter { it != null }

        // SBS signature fitting: reuse snv_counts.csv with pipeline:'sigs' for cache compatibility
        ch_sigs_for_assignment = ch_samples_to_run
            .map { meta ->
                def f = findOncoFile(meta,
                    "${meta.folder}/sigs/${meta.subject}-T.sig.snv_counts.csv",
                    'SNV counts (SigProfiler SBS)')
                f ? [meta + [pipeline: 'sigs'], f] : null
            }
            .filter { it != null }

        // DBS signature fitting: extract from PAVE somatic VCF (raw SAGE somatic in dna-only)
        ch_sigs_dbs = ch_samples_to_run
            .map { meta ->
                def f = isDnaOnly()
                    ? findOncoFile(meta,
                        "${meta.folder}/sage/somatic/${meta.subject}-T.sage.somatic.vcf.gz",
                        'somatic VCF (SigProfiler DBS, dna-only)')
                    : findOncoFile(meta,
                        "${meta.folder}/pave/${meta.subject}-T.pave.somatic.vcf.gz",
                        'somatic VCF (SigProfiler DBS)')
                f ? [meta + [pipeline: 'sigs_dbs'], f] : null
            }
            .filter { it != null }

        // ID signature fitting: extract indels from PAVE somatic VCF (raw SAGE somatic in dna-only)
        ch_sigs_id = ch_samples_to_run
            .map { meta ->
                def f = isDnaOnly()
                    ? findOncoFile(meta,
                        "${meta.folder}/sage/somatic/${meta.subject}-T.sage.somatic.vcf.gz",
                        'somatic VCF (SigProfiler ID, dna-only)')
                    : findOncoFile(meta,
                        "${meta.folder}/pave/${meta.subject}-T.pave.somatic.vcf.gz",
                        'somatic VCF (SigProfiler ID)')
                f ? [meta + [pipeline: 'sigs_id'], f] : null
            }
            .filter { it != null }

        // SIGS SNV counts CSV → mutational signature trinucleotide counts
        ch_sigs_counts = ch_samples_to_run
            .map { meta ->
                def snv_counts = findOncoFile(meta,
                    "${meta.folder}/sigs/${meta.subject}-T.sig.snv_counts.csv",
                    'mutational signature SNV counts')
                snv_counts ? [meta + [pipeline: 'sigs_counts'], snv_counts] : null
            }
            .filter { it != null }

        // ── Run subworkflows ──────────────────────────────────────────────────

        GENOMIC_CNV(ch_purple_cnv, ensembl_annotations)

        GENOMIC_SV(ch_esvee_vcf, ensembl_annotations)

        ISOFOX_FUSION_TO_CBIOPORTAL(ch_isofox_fusion)

        SIGPROFILER_SBS(ch_sigs_for_assignment, file(params.cosmic_reference), file(params.sbs_metadata))

        SIGPROFILER_DBS(ch_sigs_dbs, file(params.cosmic_reference), file(params.dbs_metadata))

        SIGPROFILER_ID(ch_sigs_id, file(params.cosmic_reference), file(params.id_metadata), fasta)

        SIGS_COUNTS_TO_CBIOPORTAL(ch_sigs_counts)

        GENOMIC_EXPRESSION(ch_isofox_exp, ensembl_annotations_expr)

        GENOMIC_MUTATIONS(
            ch_sage_germline_vcf,
            ch_sage_vcf,
            ch_sage_rna_vcf,
            fasta,
            vep_data,
            pcgr_data,
            needs_vep,
            needs_pcgr
        )

        // ── Mix new results with pre-existing cached results ──────────────────

        all_cnv_seg = GENOMIC_CNV.out.segfile
            .mix(ch_files_ran
                .filter { meta, files -> meta.pipeline == 'cnv' }
                .map { meta, files -> [meta, files.seg] })

        all_cnv_long = GENOMIC_CNV.out.longfile
            .mix(ch_files_ran
                .filter { meta, files -> meta.pipeline == 'cnv' }
                .map { meta, files -> [meta, files.longfile] })

        // Group all SV files (fresh + cached) by sample, then merge into one per-sample file
        ch_sv_all_raw = GENOMIC_SV.out.sv_out
            .mix(ISOFOX_FUSION_TO_CBIOPORTAL.out.sv)
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sv' })
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sv_rna_fusion' })

        ch_sv_per_sample = ch_sv_all_raw
            .map { meta, file -> [meta.sample, meta, file] }
            .groupTuple()
            .map { sample, metas, files ->
                [metas[0], files instanceof List ? files : [files]]
            }

        MERGE_SAMPLE_SV(ch_sv_per_sample)

        all_sv = MERGE_SAMPLE_SV.out.sv

        all_expression = GENOMIC_EXPRESSION.out.out
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'expression' })

        all_mutations = GENOMIC_MUTATIONS.out.out
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'mutation' })

        all_sigs = SIGPROFILER_SBS.out.sigs
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs' })

        all_sigs_counts = SIGS_COUNTS_TO_CBIOPORTAL.out.sigs_counts
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs_counts' })

        all_sigs_dbs = SIGPROFILER_DBS.out.sigs_dbs
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs_dbs' })

        all_sigs_counts_dbs = SIGPROFILER_DBS.out.sigs_counts_dbs
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs_counts_dbs' })

        all_sigs_id = SIGPROFILER_ID.out.sigs_id
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs_id' })

        all_sigs_counts_id = SIGPROFILER_ID.out.sigs_counts_id
            .mix(ch_files_ran.filter { meta, f -> meta.pipeline == 'sigs_counts_id' })

        // ── Aggregate per-group outputs ───────────────────────────────────────

        GENOMIC_AGGREGATE_OUTPUT(
            all_cnv_seg,
            all_cnv_long,
            all_sv,
            all_expression,
            all_mutations,
            all_sigs,
            all_sigs_counts,
            all_sigs_dbs,
            all_sigs_counts_dbs,
            all_sigs_id,
            all_sigs_counts_id,
        )

        // ── ML formatting ─────────────────────────────────────────────────────

        GENOMIC_ML(
            GENOMIC_AGGREGATE_OUTPUT.out.cnv,
            GENOMIC_AGGREGATE_OUTPUT.out.expression,
            GENOMIC_AGGREGATE_OUTPUT.out.mutation,
            GENOMIC_AGGREGATE_OUTPUT.out.sv,
            cosmic_data,
            chimer_data,
        )

        // ── Study-level metadata ──────────────────────────────────────────────

        all_groups = ch_samples.map { meta -> meta.group }.unique()

        meta_text = """type_of_cancer: ${params.cancer_type}
cancer_study_identifier: add_text
name: ${params.study_id}
description: ${params.project_description}
add_global_case_list: true
reference_genome: hg38
        """

        GENERATE_META_FILE(all_groups, "study", meta_text)

        // ── Cancer type file (optional) ──────────────────────────────────────

        ch_cancer_type = Channel.empty()
        if (params.generate_cancer_type) {
            GENERATE_CANCER_TYPE(all_groups)
            ch_cancer_type = GENERATE_CANCER_TYPE.out.flatMap { group, ct, meta_ct ->
                [tuple(group, ct), tuple(group, meta_ct)]
            }
        }

        // ── Package all cBioPortal files into a tar.gz per group ──────────────

        all_package_files = GENOMIC_AGGREGATE_OUTPUT.out.cnv
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.cnv_seg)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sv)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.expression)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.mutation)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_dbs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts_dbs)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_id)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.meta_files)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.case_files)
            .mix(GENERATE_META_FILE.out)
            .mix(ch_cancer_type)
            .groupTuple()

        PACKAGE_CBIOPORTAL(all_package_files)

        // ── Subject → tumor sample linking file ───────────────────────────────

        ch_linking_file = ch_samples
            .map    { meta -> tuple(meta.group, "${meta.subject}\t${meta.sample}") }
            .unique()
            .groupTuple()
            .map { group, lines ->
                def file_content = "subject_id\tsample_id\n" + lines.join("\n")
                def output_file  = file("${params.outdir}/${group}/util_linking_file.txt")
                output_file.parent.mkdirs()
                output_file.text = file_content
                return tuple(group, output_file)
            }

        // ── Software versions ─────────────────────────────────────────────────

        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:     'software_versions.yml',
                sort:     true,
                newLine:  true,
            )

    emit:
        linking_file = ch_linking_file  // channel<tuple(group, file)> — subject→sample linking file per group
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
