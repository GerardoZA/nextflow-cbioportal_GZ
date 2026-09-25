include { FORMAT_ML_CNV } from '../../../modules/local/format_ml_cnv'
include { FORMAT_ML_EXPRESSION } from '../../../modules/local/format_ml_expression'
include { FORMAT_ML_MUTATION } from '../../../modules/local/format_ml_mutation'
include { FORMAT_PROCESS_ML_SV } from '../../../modules/local/format_process_ml_sv'
include { PROCESS_ML_CNV } from '../../../modules/local/process_ml_cnv'
include { PROCESS_ML_EXPRESSION } from '../../../modules/local/process_ml_expression'
include { PROCESS_ML_MUTATION } from '../../../modules/local/process_ml_mutation'
include { DOWNLOAD_KNOWN_FUSIONS } from '../../../modules/local/download_known_fusions'

workflow GENOMIC_ML {
    take:
        cnv_result_long     // tuple (group, filepath)
        expression_result   // tuple (group, filepath)
        mutation_result     // tuple (group, filepath)
        sv_result           // tuple (group, filepath)
        cosmic_data         // file
        chimer_data

    main:

        // known_fusions_data = DOWNLOAD_KNOWN_FUSIONS()

        FORMAT_ML_CNV(
            cnv_result_long
        )

        PROCESS_ML_CNV(
           FORMAT_ML_CNV.out
        )

        FORMAT_ML_EXPRESSION(
            expression_result
        )

        PROCESS_ML_EXPRESSION(
           FORMAT_ML_EXPRESSION.out
        )

        FORMAT_ML_MUTATION(
            mutation_result
        )

        // Pre-staged cancerhotspots.org JSON — compute nodes may have no internet
        hotspots_json = params.hotspots_data ? file(params.hotspots_data, checkIfExists: true) : []

        PROCESS_ML_MUTATION(
            FORMAT_ML_MUTATION.out,
            hotspots_json
        )

        if (cosmic_data) {
            FORMAT_PROCESS_ML_SV(
                cosmic_data,
                chimer_data,
                sv_result
            )
        }

}
