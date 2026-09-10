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
include { resource                     } from '../subworkflows/local/utils'


// Resolve an oncoanalyser output file path; log a warning and return null if absent.
def findOncoFile(meta, path_str, label) {
    def f = file(path_str, checkIfExists: false)
    if (!f.exists() || f.isEmpty()) {
        log.warn "File not found for ${meta.subject} (${label}): ${path_str}"
        return null
    }
    return f
}


// Somatic or germline (kind) variant VCF for a subject. "both" reads PAVE. dna-only prefers
// PAVE too, and falls back to the raw SAGE VCF only when oncoanalyser wrote no PAVE file for
// the subject (it can skip PAVE/sage_append without RNA). Raw SAGE carries no PAVE INFO, so
// gnomAD_AF and the HMF.* MAF columns stay empty for a subject read from it.
def findVariantVcf(meta, kind, label) {
    def pave = "${meta.folder}/pave/${kind}/${meta.subject}-T.pave.${kind}.vcf.gz"
    def pave_f = file(pave, checkIfExists: false)
    if (params.type != 'dna-only' || (pave_f.exists() && !pave_f.isEmpty())) {
        return findOncoFile(meta, pave, "${label} (PAVE ${kind})")
    }
    return findOncoFile(meta,
        "${meta.folder}/sage/${kind}/${meta.subject}-T.sage.${kind}.vcf.gz",
        "${label} (SAGE ${kind}, dna-only fallback)")
}


// Every per-subject file this workflow publishes to <outdir>/<group>/<subject>/.
//
// SINGLE SOURCE OF TRUTH: both the "is this subject already done" gate and the
// cache loader read this table. Keeping them as two hand-written lists is how
// the ID signatures ended up gated-but-never-loaded — add a modality here and
// it is picked up by both sides at once.
//
//   tag      — pipeline tag the mix sites below select on
//   suffix   — filename, appended to meta.sample
//   required — must exist for the subject to count as processed. The signature
//              files are optional: SigProfiler legitimately emits nothing for a
//              subject with too few mutations, and gating on them would
//              reprocess such a subject on every single run.
def subjectOutputs() {
    return [
        [tag: 'cnv_seg',         suffix: '_data_cna_hg38.seg',                               required: true ],
        [tag: 'cnv_long',        suffix: '_data_cna_long.txt',                               required: true ],
        [tag: 'sv',              suffix: '.data_sv.txt',                                     required: true ],
        [tag: 'expression',      suffix: '.tpm.tsv',                                         required: true ],
        [tag: 'mutation',        suffix: '.somatic_rna_germline.maf',                        required: true ],
        [tag: 'sigs',            suffix: '.data_mutational_signatures_contribution_SBS.txt', required: false],
        [tag: 'sigs_counts',     suffix: '.data_mutational_signatures_counts_SBS.txt',       required: false],
        [tag: 'sigs_dbs',        suffix: '.data_mutational_signatures_contribution_DBS.txt', required: false],
        [tag: 'sigs_counts_dbs', suffix: '.data_mutational_signatures_counts_DBS.txt',       required: false],
        [tag: 'sigs_id',         suffix: '.data_mutational_signatures_contribution_ID.txt',  required: false],
        [tag: 'sigs_counts_id',  suffix: '.data_mutational_signatures_counts_ID.txt',        required: false],
    ]
}


// Resolve a subject's previously-published outputs.
// Returns null when the subject still has work to do, otherwise a map of
// tag -> file. All-or-nothing on purpose: a half-written subject is reprocessed
// from scratch, and none of its stale files are merged into the study.
def resolveSubjectCache(meta) {
    def subjectDir = file("${params.outdir}/${meta.group}/${meta.subject}")
    if (!subjectDir.exists() || !subjectDir.isDirectory()) {
        return null
    }
    def found = [:]
    def complete = subjectOutputs().every { output ->
        def f = file("${subjectDir}/${meta.sample}${output.suffix}", checkIfExists: false)
        if (f.exists() && !f.isEmpty()) {
            found[output.tag] = f
            return true
        }
        return !output.required
    }
    return complete ? found : null
}


// Cached per-subject files carrying the given tag, as [meta, file].
// The tag must exist in subjectOutputs(), so a typo at a mix site fails loudly
// instead of quietly yielding an empty channel — which is exactly how the ID
// signatures went missing for cached subjects in the first place.
def cachedFiles(ch_cache, String tag) {
    assert tag in subjectOutputs()*.tag : "Unknown cached-output tag '${tag}'"
    return ch_cache.filter { meta, _f -> meta.pipeline == tag }
}


// Subjects the previous run put in the study, read from the linking file it left
// behind. That file is an exact record of the last cohort, which beats scanning the
// group directory for subject folders (that directory also holds case_lists/,
// machine_learning/ and the flat data_*/meta_* files).
//
// Must be called before ch_linking_file overwrites it for this run.
def previousCohort(group) {
    def linking = file("${params.outdir}/${group}/util_linking_file.txt", checkIfExists: false)
    if (!linking.exists()) {
        return [] as Set
    }
    return linking.readLines()
        .drop(1)
        .findAll { line -> line?.trim() }
        .collect { line -> line.split('\t')[0] }
        .toSet()
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
        vep_data                // channel<path> — pre-staged VEP cache (may be empty)
        pcgr_data               // channel<path> — pre-staged PCGR reference data (may be empty)
        mafsmith_data           // channel<path> — pre-staged mafsmith home (may be empty)
        needs_vep               // boolean — true when vep_data is not supplied
        needs_pcgr              // boolean — true when pcgr_data is not supplied
        needs_mafsmith          // boolean — true when mafsmith_data is not supplied
        fasta                   // path — GRCh38 reference FASTA (for SigProfiler indel context)
        cosmic_data             // channel<path> — COSMIC/ChimerKB fusion data for ML step
        chimer_data
        hotspots_data           // channel<path> — pre-staged cancerhotspots.org data (may be empty)

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
        // A subject counts as done only when every required per-subject output is
        // already on disk (see subjectOutputs()). Done subjects contribute their
        // cached files and run no tasks; every other subject is reprocessed from
        // source and contributes nothing from disk, so a partially written subject
        // can never be merged in twice.

        ch_subject_cache = ch_samples
            .map { meta -> [meta, resolveSubjectCache(meta)] }

        // Cached per-subject files, tagged so the mix sites below can select them.
        ch_files_ran = ch_subject_cache
            .filter { _meta, cache -> cache != null }
            .flatMap { meta, cache ->
                cache.collect { tag, cached_file -> [meta + [pipeline: tag], cached_file] }
            }

        // Subjects that still need processing
        ch_samples_to_run = ch_subject_cache
            .filter { _meta, cache -> cache == null }
            .map { meta, _cache -> meta }

        // Log skipped samples
        ch_subject_cache
            .filter { _meta, cache -> cache != null }
            .subscribe { meta, cache ->
                log.info "Skipping already-processed subject: ${meta.subject} (${cache.size()} cached file(s))"
            }

        // One-line run summary, and a warning when there is nothing left to do
        ch_subject_cache
            .map { _meta, cache -> cache == null ? 'run' : 'cached' }
            .toList()
            .subscribe { states ->
                def to_run = states.count('run')
                log.info "Incremental processing: ${states.count('cached')} subject(s) cached, ${to_run} to process"
                if (to_run == 0) {
                    log.warn "All subjects in samplesheet already processed. Only the group-level merge will run."
                }
            }

        // The samplesheet defines the cohort: the group-level study files are rebuilt
        // from it on every run, so a subject dropped from the samplesheet disappears
        // from the study even though its per-subject folder is still on disk. Read the
        // previous cohort now, before ch_linking_file rewrites it further down.
        previous_cohort = previousCohort("${params.study_id}")

        ch_samples
            .map { meta -> meta.subject }
            .toList()
            .subscribe { subjects ->
                def dropped = (previous_cohort - (subjects as Set)).sort()
                if (dropped) {
                    log.warn "Subject(s) ${dropped.join(', ')} were in the previous run of study " +
                             "${params.study_id} but are not in this samplesheet. The samplesheet defines the " +
                             "cohort, so they will NOT appear in the merged study files, even though their " +
                             "output folders remain on disk."
                }
            }

        // ── Build per-modality input channels ─────────────────────────────────
        // File naming conventions (relative to folder), oncoanalyser 3.0 layout:
        //   pave/somatic/${subject}-T.pave.somatic.vcf.gz          — PAVE somatic DNA
        //   pave/germline/${subject}-T.pave.germline.vcf.gz        — PAVE germline DNA
        //   sage_append/${subject}-T/${subject}-T.sage.append.vcf.gz — SAGE somatic RNA append
        //   esvee/${subject}-T.esvee.somatic.vcf.gz
        //   purple/${subject}-T.purple.cnv.somatic.tsv
        //   purple/${subject}-T.purple.cnv.gene.tsv
        //   isofox/${subject}-T.isf.gene_data.tsv
        //   isofox/${subject}-T.isf.pass_fusions.tsv
        //   sigs/${subject}-T.sig.snv_counts.csv

        // PAVE somatic VCF → mutations (dna-only: raw SAGE when pave/ is absent, see findVariantVcf)
        ch_sage_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = findVariantVcf(meta, 'somatic', 'mutation')
                vcf ? [meta + [pipeline: 'mutation'], vcf] : null
            }
            .filter { it != null }

        // PAVE germline VCF → germline mutations (same dna-only fallback as above)
        ch_sage_germline_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = findVariantVcf(meta, 'germline', 'germline mutation')
                vcf ? [meta + [pipeline: 'mutation_germline'], vcf] : null
            }
            .filter { it != null }

        // SAGE RNA-append VCF → somatic RNA mutations (skipped entirely in dna-only mode)
        ch_sage_rna_vcf = params.type == 'dna-only'
            ? channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def vcf = findOncoFile(meta,
                        "${meta.folder}/sage_append/${meta.subject}-T/${meta.subject}-T.sage.append.vcf.gz",
                        'mutation (SAGE RNA append)')
                    vcf ? [meta + [pipeline: 'mutation_rna'], vcf] : null
                }
                .filter { it != null }

        // PURPLE CNV somatic + gene TSV → copy-number
        // The purity file is optional: it only supplies the sample's sex for the
        // chrX/chrY baseline, so a subject without one still gets CNV output —
        // just with a diploid baseline everywhere, same as before this file existed.
        ch_purple_cnv = ch_samples_to_run
            .map { meta ->
                def somatic = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.somatic.tsv",
                    'cnv (PURPLE somatic)')
                def gene = findOncoFile(meta,
                    "${meta.folder}/purple/${meta.subject}-T.purple.cnv.gene.tsv",
                    'cnv (PURPLE gene)')
                def purity = file("${meta.folder}/purple/${meta.subject}-T.purple.purity.tsv", checkIfExists: false)
                if (!purity.exists() || purity.isEmpty()) {
                    log.warn "No PURPLE purity file for ${meta.subject}; chrX/chrY will use a diploid baseline"
                    purity = []
                }
                (somatic && gene) ? [meta + [pipeline: 'cnv'], somatic, gene, purity] : null
            }
            .filter { it != null }

        // ESVEE somatic VCF (tumor only) → structural variants
        ch_esvee_vcf = ch_samples_to_run
            .map { meta ->
                def vcf = findOncoFile(meta,
                    "${meta.folder}/esvee/${meta.subject}-T.esvee.somatic.vcf.gz",
                    'sv (ESVEE tumor)')
                vcf ? [meta + [pipeline: 'sv'], vcf] : null
            }
            .filter { it != null }

        // Isofox gene expression CSV → TPM (skipped entirely in dna-only mode)
        ch_isofox_exp = params.type == 'dna-only'
            ? channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def exp = findOncoFile(meta,
                        "${meta.folder}/isofox/${meta.subject}-T.isf.gene_data.tsv",
                        'expression (Isofox)')
                    exp ? [meta + [pipeline: 'expression'], exp] : null
                }
                .filter { it != null }

        // Isofox pass_fusions CSV (tumor RNA) → RNA fusions for data_sv.txt (skipped entirely in dna-only mode)
        ch_isofox_fusion = params.type == 'dna-only'
            ? channel.empty()
            : ch_samples_to_run
                .map { meta ->
                    def fusions = findOncoFile(meta,
                        "${meta.folder}/isofox/${meta.subject}-T.isf.pass_fusions.tsv",
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

        // DBS signature fitting: extract from PAVE somatic VCF (dna-only: raw SAGE fallback)
        ch_sigs_dbs = ch_samples_to_run
            .map { meta ->
                def f = findVariantVcf(meta, 'somatic', 'SigProfiler DBS')
                f ? [meta + [pipeline: 'sigs_dbs'], f] : null
            }
            .filter { it != null }

        // ID signature fitting: extract indels from PAVE somatic VCF (dna-only: raw SAGE fallback)
        ch_sigs_id = ch_samples_to_run
            .map { meta ->
                def f = findVariantVcf(meta, 'somatic', 'SigProfiler ID')
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

        SIGPROFILER_SBS(ch_sigs_for_assignment, file(resource('sbs_signatures')), file(resource('sbs_metadata')))

        SIGPROFILER_DBS(ch_sigs_dbs, file(resource('dbs_signatures')), file(resource('dbs_metadata')))

        SIGPROFILER_ID(ch_sigs_id, file(resource('id_signatures')), file(resource('id_metadata')), fasta)

        SIGS_COUNTS_TO_CBIOPORTAL(ch_sigs_counts)

        GENOMIC_EXPRESSION(ch_isofox_exp, ensembl_annotations)

        GENOMIC_MUTATIONS(
            ch_sage_germline_vcf,
            ch_sage_vcf,
            ch_sage_rna_vcf,
            vep_data,
            pcgr_data,
            mafsmith_data,
            needs_vep,
            needs_pcgr,
            needs_mafsmith
        )

        // ── Mix new results with pre-existing cached results ──────────────────
        // cachedFiles() rejects a tag that is not in subjectOutputs(), so a cached
        // subject can never go missing from one of these outputs unnoticed.

        all_cnv_seg = GENOMIC_CNV.out.segfile
            .mix(cachedFiles(ch_files_ran, 'cnv_seg'))

        all_cnv_long = GENOMIC_CNV.out.longfile
            .mix(cachedFiles(ch_files_ran, 'cnv_long'))

        // Group the freshly produced SV files (DNA + RNA fusion) by sample and merge
        // them into one file per sample. A cached subject already has that merged
        // file on disk, so it bypasses MERGE_SAMPLE_SV entirely — re-merging it would
        // feed the process an input with the same name as its own output.
        ch_sv_per_sample = GENOMIC_SV.out.sv_out
            .mix(ISOFOX_FUSION_TO_CBIOPORTAL.out.sv)
            .map { meta, sv_file -> [meta.sample, meta, sv_file] }
            .groupTuple()
            .map { _sample, metas, files ->
                [metas[0], files instanceof List ? files : [files]]
            }

        MERGE_SAMPLE_SV(ch_sv_per_sample)

        all_sv = MERGE_SAMPLE_SV.out.sv
            .mix(cachedFiles(ch_files_ran, 'sv'))

        all_expression = GENOMIC_EXPRESSION.out.out
            .mix(cachedFiles(ch_files_ran, 'expression'))

        all_mutations = GENOMIC_MUTATIONS.out.out
            .mix(cachedFiles(ch_files_ran, 'mutation'))

        all_sigs = SIGPROFILER_SBS.out.sigs
            .mix(cachedFiles(ch_files_ran, 'sigs'))

        all_sigs_counts = SIGS_COUNTS_TO_CBIOPORTAL.out.sigs_counts
            .mix(cachedFiles(ch_files_ran, 'sigs_counts'))

        all_sigs_dbs = SIGPROFILER_DBS.out.sigs_dbs
            .mix(cachedFiles(ch_files_ran, 'sigs_dbs'))

        all_sigs_counts_dbs = SIGPROFILER_DBS.out.sigs_counts_dbs
            .mix(cachedFiles(ch_files_ran, 'sigs_counts_dbs'))

        all_sigs_id = SIGPROFILER_ID.out.sigs_id
            .mix(cachedFiles(ch_files_ran, 'sigs_id'))

        all_sigs_counts_id = SIGPROFILER_ID.out.sigs_counts_id
            .mix(cachedFiles(ch_files_ran, 'sigs_counts_id'))

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

        // ── ML formatting (opt-in: --run_ml) ──────────────────────────────────

        if (params.run_ml) {
            GENOMIC_ML(
                GENOMIC_AGGREGATE_OUTPUT.out.cnv,
                GENOMIC_AGGREGATE_OUTPUT.out.expression,
                GENOMIC_AGGREGATE_OUTPUT.out.mutation,
                GENOMIC_AGGREGATE_OUTPUT.out.sv,
                cosmic_data,
                chimer_data,
                hotspots_data,
            )
        }

        // ── Study-level metadata + cancer type file ──────────────────────────
        // Skipped on a follow-up load (--incremental): cBioPortal already holds
        // the study and cancer type from the first load.

        all_groups = ch_samples.map { meta -> meta.group }.unique()

        ch_study_files = Channel.empty()
        if (!params.incremental) {
            meta_text = """type_of_cancer: ${params.cancer_type}
cancer_study_identifier: add_text
name: ${params.study_id}
description: ${params.project_description}
add_global_case_list: true
reference_genome: hg38
        """

            GENERATE_META_FILE(all_groups, "study", meta_text)
            GENERATE_CANCER_TYPE(all_groups)
            ch_study_files = GENERATE_META_FILE.out
                .mix(GENERATE_CANCER_TYPE.out.flatMap { group, ct, meta_ct ->
                    [tuple(group, ct), tuple(group, meta_ct)]
                })
        }

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

        // ── Files that make up the cBioPortal study package ───────────────────
        // Packaging itself happens in main.nf, so that in "both" mode the archive
        // also picks up the clinical and timeline files produced by CLINICAL.

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
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.sigs_counts_id)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.meta_files)
            .mix(GENOMIC_AGGREGATE_OUTPUT.out.case_files)
            .mix(ch_study_files)
            .mix(ch_linking_file)

        // ── Software versions ─────────────────────────────────────────────────

        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name:     'software_versions.yml',
                sort:     true,
                newLine:  true,
            )

    emit:
        linking_file  = ch_linking_file     // channel<tuple(group, file)> — subject→sample linking file per group
        package_files = all_package_files   // channel<tuple(group, file)> — files to put in the study archive
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
