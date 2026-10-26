//
// Subworkflow with functionality specific to the CRCHUM-CITADEL/nextflow-cbioportal pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           	// boolean: Display version and exit
    monochrome_logs   	// boolean: Do not use coloured log outputs
    nextflow_cli_args 	// array: List of positional nextflow CLI args
    mode              	// string: pipeline mode [clinical, genomic, both]
    outdir            	// string: The output directory where the results will be saved
    genomic_input     	// string: Path to input samplesheet
    clinical_input    	// string: Path to input samplesheet
    study_id          	// string: Study identifier for metadata and top-level output folder
	project_description // string : String of project description to be put in metadata

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //

	// custom function validation
    validateInputParameters()

	// nf-schema validation
    validateParameters()

    //
    // Create channel from input file provided through params.input
    //
    ch_genomic_samplesheet  = Channel.empty()
    ch_clinical_samplesheet = Channel.empty()

    if (mode in ['genomic', 'both']){
        if (!resource('ensembl_annotations')){
            error("ERROR: Missing ensembl_annotations (BioMart TSV). Set --resources_dir to the resource kit, pass --ensembl_annotations, or use -profile citadel at CRCHUM.")
        }

        ch_genomic_samplesheet = Channel.fromList(samplesheetToList(genomic_input, "assets/schema_genomic_input.json"))
    }

    if (mode in ['clinical', 'both']){
        ch_clinical_samplesheet = clinical_input
            ? Channel.fromList(samplesheetToList(clinical_input, "assets/schema_clinical_input.json"))
            : Channel.empty()
    }

	if (study_id == "" ) {
		log.warn "study_id not set. 'test_name' will be used."
		study_id = "test_name"
	}

	if (project_description == "") {
		log.warn "project description not set. 'test_description' will be used"
		project_description = 'test_description'
	}


    emit:
    genomic_samplesheet  = ch_genomic_samplesheet
    clinical_samplesheet = ch_clinical_samplesheet
	name 		= study_id
	description = project_description
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = [:]

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            // TODO: wait for HPC access
            // completionEmail(
            //     summary_params,
            //     email,
            //     email_on_fail,
            //     plaintext_email,
            //     outdir,
            //     monochrome_logs,
            //     []
            // )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // check modes and input
    if (!params.mode){
        error("ERROR: Pipeline mode not chosen in configuration file. Choices : 'genomic', 'clinical', or 'both'")
    }
    params.mode = params.mode.toLowerCase()
    if ( !params.mode in ['genomic','clinical','both'] ) {
        error("Error: Invalid pipeline mode chosen. Choices : 'genomic', 'clinical', or 'both'")
    }

    // check sequencing data type
    // NOTE: params values supplied on the command line are immutable — assigning to
    // params.type here is silently ignored, so the value is NOT normalised in place.
    // Every consumer must compare case-insensitively (see isDnaOnly() in the workflows).
    def type_value = (params.type ?: 'both').toString().toLowerCase()
    if ( !(type_value in ['both', 'dna-only']) ) {
        error("ERROR: Invalid --type value '${params.type}'. Choices: 'both' or 'dna-only'")
    }

    // make sure there's input
    if (params.mode in ["genomic", "both"] && !params.genomic_samplesheet){
        error("ERROR: Could not find genomic samplesheet. Not running any tests. Check input in nextflow.config")
    }

    // The resource kit must match the version this release was validated against.
    if (params.resources_dir) {
        def kit = file(params.resources_dir)
        if (!kit.isDirectory()) {
            error("ERROR: resources_dir does not exist or is not a directory: ${params.resources_dir}")
        }
        def version_file = file("${params.resources_dir}/VERSION")
        def kit_version  = version_file.exists() ? version_file.text.trim() : 'none'
        if (kit_version != params.resources_version) {
            error("ERROR: resources_dir holds resource kit version '${kit_version}', but this pipeline release needs '${params.resources_version}'. Download nextflow-cbioportal-resources-v${params.resources_version} (see docs/resources.md).")
        }

        // Warn about kit files that are missing and not overridden by a param.
        def absent = resourceLayout().findAll { name, rel ->
            !explicitParam(name) && !file("${params.resources_dir}/${rel}").exists()
        }.keySet()
        if (absent) {
            log.warn "Resource kit ${params.resources_dir} is missing, and no param overrides: ${absent.join(', ')}"
        }
    }

    if (params.mode in ["genomic", "both"]) {
        def genome_reference = resource('genome_reference')
        if (!genome_reference) {
            error("ERROR: genome_reference is required for genomic mode. Set --resources_dir to the resource kit, or pass --genome_reference.")
        }

        if (!file(genome_reference).exists()) {
            error("ERROR: Genome reference file does not exist: ${genome_reference}")
        }

        requireResources(['sbs_signatures', 'dbs_signatures', 'id_signatures', 'sbs_metadata', 'dbs_metadata', 'id_metadata'], 'genomic mode')

        // FORMAT_PROCESS_ML_SV unions the COSMIC and ChimerKB fusion lists, so the two
        // are all-or-nothing. Setting only one used to reach the process with an empty
        // path input and abort the run with an unattributable "Path must not be empty".
        // Only checked when ML runs; otherwise the fusion files are never read.
        if (params.run_ml && params.cosmic_data) {
            def cosmic_data = file(params.cosmic_data)
            if (!cosmic_data.exists()) {
                error("ERROR: Cosmic data file does not exist: ${params.cosmic_data}")
            }

            def chimer_data = resource('chimer_data')
            if (!chimer_data) {
                error("ERROR: cosmic_data is set but chimer_data is empty. The ML SV step needs both fusion databases — set chimer_data (or --resources_dir), or clear cosmic_data to skip FORMAT_PROCESS_ML_SV.")
            }

            if (!file(chimer_data).exists()) {
                error("ERROR: Chimer data file does not exist: ${chimer_data}")
            }
        }

    }

    if (params.mode in ["clinical", "both"] && !params.clinical_samplesheet){
        log.warn "No clinical samplesheet provided. Template clinical files will be generated from the linking file."
    }

    // MOHCCN maps are only read when a clinical samplesheet is given.
    if (params.mode in ["clinical", "both"] && params.clinical_samplesheet) {
        requireResources(['mohccn_primary_site_map', 'mohccn_specimen_tissue_source_map', 'mohccn_treatment_intent_map'], 'clinical output')
    }

    if (params.mode == "clinical" && !params.clinical_samplesheet && params.sample_registrations) {
        def regs_file = file(params.sample_registrations)
        if (!regs_file.exists()) {
            error("ERROR: sample_registrations file does not exist: ${params.sample_registrations}")
        }
    }


}


// Path of every reference file inside the resource kit (cosmic_data is licence-gated, never in the kit).
def resourceLayout() {
    def sigs = 'genomic/cosmic_mutational_signatures'
    def moh  = 'clinical/MoH/dictionary'
    return [
        ensembl_annotations               : 'genomic/annotations/biomart_grch38_ensembl_113_with_entrez_id.tsv',
        chimer_data                       : 'genomic/annotations/ChimerKB4.xlsx',
        cancer_hotspots_data              : 'genomic/annotations/cancerhotspots_single.json',
        vep_data                          : 'genomic/vep/cache',
        pcgr_data                         : 'genomic/pcgr',
        mafsmith_data                     : 'genomic/mafsmith/mafsmith_0.1.0',
        genome_reference                  : 'genomic/reference/Homo_sapiens_assembly38.fasta',
        sbs_signatures                    : "${sigs}/COSMIC_Human_SBS-96_GRCh38_v3.6.csv",
        dbs_signatures                    : "${sigs}/COSMIC_Human_DBS-78_GRCh38_v3.6.csv",
        id_signatures                     : "${sigs}/COSMIC_Human_ID-83_GRCh38_v3.6.csv",
        sbs_metadata                      : "${sigs}/cosmic_sbs_metadata.tsv",
        dbs_metadata                      : "${sigs}/cosmic_dbs_metadata.tsv",
        id_metadata                       : "${sigs}/cosmic_id_metadata.tsv",
        mohccn_primary_site_map           : "${moh}/mohccn_clinical_data_modelv3-1_sep2024_primary_site.csv",
        mohccn_specimen_tissue_source_map : "${moh}/mohccn_clinical_data_modelv3-1_sep2024_specimen_tissue_source.csv",
        mohccn_treatment_intent_map       : "${moh}/mohccn_clinical_data_modelv3-1_sep2024_treatment_intent.csv",
    ]
}


// Explicit param, else the kit file if it exists, else ''; resolved at run time so a profile-set resources_dir is seen.
def resource(String name) {
    assert name in resourceLayout() : "Unknown resource '${name}'"
    if (explicitParam(name)) {
        return explicitParam(name).toString()
    }
    if (!params.resources_dir) {
        return ''
    }
    def path = "${params.resources_dir}/${resourceLayout()[name]}".toString()
    return file(path).exists() ? path : ''
}


// Kit-backed params are not declared in nextflow.config; containsKey avoids the "undefined parameter" warning.
def explicitParam(String name) {
    return params.containsKey(name) ? params.get(name) : null
}


def helpRequested() {
    return ['help', 'helpFull'].any { name -> explicitParam(name) }
}


// Fail up front, naming every missing reference file.
def requireResources(List names, String purpose) {
    def missing = names.findAll { name ->
        def path = resource(name)
        !path || !file(path).exists()
    }
    if (missing) {
        def detail = missing.collect { name -> "  ${name}: ${resource(name) ?: '(not set)'}" }.join('\n')
        error("ERROR: Reference files needed for ${purpose} are missing. Set --resources_dir to the resource kit, or pass each one explicitly:\n${detail}")
    }
}
