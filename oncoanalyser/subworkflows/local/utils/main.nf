//
// Subworkflow with functionality specific to the CRCHUM-CITADEL/nextflow-sante-precision pipeline
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
        if (!params.ensembl_annotations){
            error("ERROR: Missing ensembl_annotations file (tsv format). Check input in nextflow.config")
        }

        if (!params.ensembl_annotations_expr){
            error("ERROR: Missing ensembl_annotations_expr file (tsv format). Check input in nextflow.config")
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

    if (params.mode in ["genomic", "both"]) {
        if (!params.genome_reference) {
            error("ERROR: genome_reference parameter is required for genomic mode")
        }

        def genome_file = file(params.genome_reference)

        if (!genome_file.exists()) {
            error("ERROR: Genome reference file does not exist: ${params.genome_reference}")
        }

        if (params.cosmic_data) {
            def cosmic_data = file(params.cosmic_data)
            if (!cosmic_data.exists()) {
                error("ERROR: Cosmic data file does not exist: ${params.cosmic_data}")
            }
        }

        if (params.hotspots_data) {
            def hotspots_data = file(params.hotspots_data)
            if (!hotspots_data.exists()) {
                error("ERROR: Hotspots data file does not exist: ${params.hotspots_data}")
            }
        }

    }

    if (params.mode in ["clinical", "both"] && !params.clinical_samplesheet){
        log.warn "No clinical samplesheet provided. Template clinical files will be generated from the linking file."
    }

    if (params.mode == "clinical" && !params.clinical_samplesheet && params.sample_registrations) {
        def regs_file = file(params.sample_registrations)
        if (!regs_file.exists()) {
            error("ERROR: sample_registrations file does not exist: ${params.sample_registrations}")
        }
    }


}
