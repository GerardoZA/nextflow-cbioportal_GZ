include { GENERATE_CASE_LIST } from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'
include { MERGE_EXPRESSION_FILES_TO_CBIOPORTAL } from '../../../modules/local/merge_expression_files_to_cbioportal'
include { MERGE_SIGS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_to_cbioportal'
include { MERGE_SIGS_COUNTS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_to_cbioportal'
include { MERGE_SIGS_DBS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_dbs_to_cbioportal'
include { MERGE_SIGS_COUNTS_DBS_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_dbs_to_cbioportal'
include { MERGE_SIGS_ID_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_id_to_cbioportal'
include { MERGE_SIGS_COUNTS_ID_TO_CBIOPORTAL } from '../../../modules/local/merge_sigs_counts_id_to_cbioportal'

    workflow GENOMIC_AGGREGATE_OUTPUT {

        take:
            cnv_results_seg // channel [long: [meta, seg_file], seg [meta, seg_file]]
            cnv_results_long
            sv_results
            expression_results
            mutation_results
            sigs_results
            sigs_counts_results
            sigs_dbs_results
            sigs_counts_dbs_results
            sigs_id_results
            sigs_counts_id_results

        main:

            ch_versions = channel.empty()

            // to get all groups, just take .seg files (we assume seg and long are the same)
            all_groups = cnv_results_seg.map {meta, sample -> meta.group}.unique()

            // merge cnv ----------------------------------------
            cnv_seg_output = cnv_results_seg
                 .map {meta, filepath -> [meta.group, filepath]}
                 .groupTuple()
                 .flatMap { group, files ->
                    files.collect { filepath -> [group, filepath]}
                 }
                 .collectFile(storeDir: "${params.outdir}",
                             keepHeader : true,
                             skip: 1,
                             sort: 'deep') { group, filepath ->
                                 ["${group}/data_cna_hg38.seg", filepath.text]
                             }
                 .map { filepath -> tuple(filepath.parent.name, filepath) }

             cnv_long_output = cnv_results_long
                 .map {meta, filepath -> [meta.group, filepath]}
                 .groupTuple()
                 .flatMap {group, files ->
                     files.collect { filepath -> [group, filepath]}
                 }
                 .collectFile(storeDir : "${params.outdir}",
                              keepHeader : true,
                              skip : 1,
                              sort: 'deep') { group, filepath ->
                                 ["${group}/data_cna_long.txt", filepath.text]
                              }
                .map {filepath -> tuple(filepath.parent.name, filepath) }

            // merge sv -------------------------------------------------------------
            sv_output = sv_results
                .map {meta, filepath -> [meta.group, filepath]}
                .groupTuple()
                .flatMap { group, files ->
                    files.collect { filepath -> [group, filepath]}
                }
                .collectFile(storeDir: "${params.outdir}",
                            keepHeader : true,
                            skip: 1,
                            sort : 'deep') { group, filepath ->
                                ["${group}/data_sv.txt", filepath.text]
                            }
                .map {filepath -> tuple(filepath.parent.name, filepath) }

            // merge expression with merging not possible within pure nextflow. This will publish the file to output via the module
            tpm_file_list = expression_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]  // Take first meta since they share the same group
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    if (sortedFiles.size() < 2) {
                        log.warn "GENOMIC_EXPRESSION: Found ${sortedFiles.size()} TPM file(s) for group ${group}. Need at least 2 files to merge. Skipping merge step."
                        return null
                    }
                    return tuple(meta, sortedFiles)
                }
                .filter { it != null }

            expression_output = MERGE_EXPRESSION_FILES_TO_CBIOPORTAL(
                tpm_file_list
            )

            // merge mutations
            mutation_output = mutation_results
                .map {meta, filepath -> [meta.group, filepath]}
                .groupTuple()
                .flatMap { group, files ->
                    files.collect { filepath -> [group, filepath]}
                }
                .collectFile(storeDir: "${params.outdir}",
                           keepHeader : true,
                           skip: 2,
                           sort: 'deep') { group, filepath ->
                              ["${group}/data_mutations_dna_rna_germline.txt", filepath.text]
                           }
                .map {filepath -> tuple(filepath.parent.name, filepath) }


            // merge signatures -------------------------------------------------------------
            sigs_file_list = sigs_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_output = MERGE_SIGS_TO_CBIOPORTAL(sigs_file_list)

            // merge signature counts -------------------------------------------------------
            sigs_counts_file_list = sigs_counts_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_counts_output = MERGE_SIGS_COUNTS_TO_CBIOPORTAL(sigs_counts_file_list)

            // merge DBS signatures --------------------------------------------------------
            sigs_dbs_file_list = sigs_dbs_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_dbs_output = MERGE_SIGS_DBS_TO_CBIOPORTAL(sigs_dbs_file_list)

            // merge DBS signature counts --------------------------------------------------
            sigs_counts_dbs_file_list = sigs_counts_dbs_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_counts_dbs_output = MERGE_SIGS_COUNTS_DBS_TO_CBIOPORTAL(sigs_counts_dbs_file_list)

            // merge ID signatures ---------------------------------------------------------
            sigs_id_file_list = sigs_id_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_id_output = MERGE_SIGS_ID_TO_CBIOPORTAL(sigs_id_file_list)

            // merge ID signature counts ---------------------------------------------------
            sigs_counts_id_file_list = sigs_counts_id_results
                .map { meta, filepath -> tuple(meta.group, meta, filepath) }
                .groupTuple()
                .map { group, metas, files ->
                    def meta = metas[0]
                    def sortedFiles = files.sort { a, b ->
                        def na = a.toString().split(/[\/\\]/).last()
                        def nb = b.toString().split(/[\/\\]/).last()
                        na <=> nb
                    }
                    return tuple(meta, sortedFiles)
                }

            sigs_counts_id_output = MERGE_SIGS_COUNTS_ID_TO_CBIOPORTAL(sigs_counts_id_file_list)

            // create meta files and case lists ---------------------------------------------------------

            cnv_sample_list = cnv_results_seg
                .map {meta, filepath -> meta.sample}
                .collect()
                .map { it.sort(false).join('\t') }

            mutation_sample_list = mutation_results
                .map {meta, filepath -> meta.sample}
                .collect()
                .map { it.sort(false).join('\t') }

            sv_sample_list = sv_results
                .map {meta, filepath -> meta.sample}
                .collect()
                .map { it.sort(false).join('\t') }

            case_name_all = channel.of("cnv", "sequenced", "sv")
            case_sample_lists = cnv_sample_list.concat(mutation_sample_list).concat(sv_sample_list)
            all_groups_cases = all_groups.combine(case_name_all).map { g, _c -> g }

            GENERATE_CASE_LIST(
                all_groups_cases,
                case_name_all,
                case_sample_lists
            )

            // to get all groups, just take .seg files (we assume seg and long are the same)
            // add_text is the key word to replac with the group.
            meta_text_seg = """cancer_study_identifier: add_text
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: SEG
reference_genome_id: hg38
description: Somatic CNA data (copy number segment file)
data_filename: data_cna_hg38.seg
            """

            meta_text_long = """cancer_study_identifier: add_text
genetic_alteration_type: COPY_NUMBER_ALTERATION
datatype: DISCRETE_LONG
stable_id: cna
show_profile_in_analysis_tab: TRUE
profile_name: Copy-number alterations
profile_description: ADD TEXT
data_filename: data_cna_long.txt
            """

            meta_text_sv = """cancer_study_identifier: add_text
genetic_alteration_type: STRUCTURAL_VARIANT
datatype: SV
stable_id: structural_variants
show_profile_in_analysis_tab: true
profile_name: Structural variants from DNA
profile_description: Structural Variant Data DNA
data_filename: data_sv.txt
            """

            meta_text_expression = """cancer_study_identifier: add_text
genetic_alteration_type: MRNA_EXPRESSION
datatype: CONTINUOUS
stable_id: rna_seq_v2_mrna
show_profile_in_analysis_tab: true
profile_name: mRNA expression (RNA-Seq TPM)
profile_description: Expression levels (RNA-Seq TPM values)
data_filename: data_expression.txt
            """

            meta_text_mutations = """cancer_study_identifier: add_text
genetic_alteration_type: MUTATION_EXTENDED
stable_id: mutations
datatype: MAF
show_profile_in_analysis_tab: true
profile_description: ADD TEXT
profile_name: Mutations
data_filename: data_mutations_dna_rna_germline.txt
    """

            meta_text_sigs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_SBS
show_profile_in_analysis_tab: true
profile_name: Mutational Signatures
profile_description: Mutational signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_SBS.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
    """

            meta_text_counts = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_SBS
show_profile_in_analysis_tab: false
profile_name: Mutational Signatures Counts
profile_description: Mutational signature trinucleotide counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_SBS.txt
generic_entity_meta_properties: NAME
    """

            meta_text_sigs_dbs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_DBS
show_profile_in_analysis_tab: true
profile_name: DBS Mutational Signatures
profile_description: Doublet base substitution signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_DBS.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
"""

            meta_text_counts_dbs = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_DBS
show_profile_in_analysis_tab: false
profile_name: DBS Mutational Signatures Counts
profile_description: Doublet base substitution counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_DBS.txt
generic_entity_meta_properties: NAME
    """

            meta_text_sigs_id = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_contribution_ID
show_profile_in_analysis_tab: true
profile_name: ID Mutational Signatures
profile_description: Indel signatures based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_contribution_ID.txt
generic_entity_meta_properties: NAME,DESCRIPTION
pivot_threshold_value: 0.0
    """

            meta_text_counts_id = """cancer_study_identifier: add_text
genetic_alteration_type: GENERIC_ASSAY
generic_assay_type: MUTATIONAL_SIGNATURE
datatype: LIMIT-VALUE
stable_id: mutational_signatures_counts_ID
show_profile_in_analysis_tab: false
profile_name: ID Mutational Signatures Counts
profile_description: Indel classification counts based on COSMIC v3.6 (WGS)
data_filename: data_mutational_signatures_counts_ID.txt
generic_entity_meta_properties: NAME
    """

        // expression meta file is only emitted in DNA+RNA ("both") mode — a dna-only study
        // never produces data_expression.txt, so a meta file pointing at it would be an
        // orphan entry that breaks cBioPortal import
        meta_entries = [
            ["cna_hg38", meta_text_seg],
            ["cna_long", meta_text_long],
            ["sv", meta_text_sv],
        ]
        if (params.type != 'dna-only') {
            meta_entries << ["expression", meta_text_expression]
        }
        meta_entries += [
            ["sequenced", meta_text_mutations],
            ["mutational_signatures_contribution_SBS", meta_text_sigs],
            ["mutational_signatures_counts_SBS", meta_text_counts],
            ["mutational_signatures_contribution_DBS", meta_text_sigs_dbs],
            ["mutational_signatures_counts_DBS", meta_text_counts_dbs],
            ["mutational_signatures_contribution_ID", meta_text_sigs_id],
            ["mutational_signatures_counts_ID", meta_text_counts_id],
        ]

        file_name_all = channel.fromList(meta_entries.collect { it[0] })
        meta_text_all  = channel.fromList(meta_entries.collect { it[1] })
        all_groups_meta = all_groups.combine(file_name_all).map { g, _f -> g }

        GENERATE_META_FILE(
            all_groups_meta,
            file_name_all,
            meta_text_all
        )

    emit:
        cnv         = cnv_long_output
        cnv_seg     = cnv_seg_output
        expression  = expression_output
        mutation    = mutation_output
        sv          = sv_output
        sigs            = sigs_output
        sigs_counts     = sigs_counts_output
        sigs_dbs        = sigs_dbs_output
        sigs_counts_dbs = sigs_counts_dbs_output
        sigs_id         = sigs_id_output
        sigs_counts_id  = sigs_counts_id_output
        meta_files      = GENERATE_META_FILE.out
        case_files  = GENERATE_CASE_LIST.out
}
