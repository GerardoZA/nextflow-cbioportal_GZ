process GENERATE_CLINICAL_SAMPLE {
    tag "$group"
    label 'process_single'

    container params.container_python

    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    input:
        tuple val(group), val(sample_lines)   // newline-separated "subject_id\tsample_id" lines

    output:
        tuple val(group), path("data_clinical_sample.txt"), emit: clinical_sample

    when:
        task.ext.when == null || task.ext.when

    script:
    def cancer_type = params.cancer_type ?: ''
    // heredoc lines sit at column 0: the multi-line value stops Nextflow stripping the indent
    """
    printf 'subject_id\\tsample_id\\n' > linking.tsv
cat << 'EOF' >> linking.tsv
${sample_lines}
EOF

    python3 ${projectDir}/bin/gen_clinical_sample_template.py \\
        --linking     linking.tsv \\
        --cancer-type '${cancer_type}' \\
        --output      data_clinical_sample.txt
    """

    stub:
    """
    touch data_clinical_sample.txt
    """
}
