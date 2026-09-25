process PROCESS_ML_MUTATION {
    publishDir { "${params.outdir}/${group}/machine_learning/processed" }, mode: 'copy'
    label "process_medium"

    container params.container_r

    input:
        tuple val(group), path(mutation_results)
        path hotspots_json   // optional: [] = fetch from cancerhotspots.org (needs internet)

    output:
        path "mutations_processed_*.tsv"

    script:
    def hotspots_arg = hotspots_json ? "${hotspots_json}" : ""
    """
    ml_mutation_processor.R $mutation_results $hotspots_arg
    """

    stub:
    """
    touch mutations_processed_somatic.tsv
    """
}
