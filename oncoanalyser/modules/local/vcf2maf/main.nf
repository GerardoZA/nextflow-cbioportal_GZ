// modified to allow vcf.gz and vcf (10/10/2025)
process VCF2MAF {
    tag "$meta.sample"
    label 'process_medium_memory'
    container "${params.container_vcf2maf}"

    input:
        tuple val(meta), path(vcf) // Now accepts both compressed (.vcf.gz) and uncompressed (.vcf) files
        path fasta                 // Required
        path vep_data              // Required for VEP running. A default of /.vep is supplied.

    output:
        tuple val(meta), path("${meta.sample}.maf"), emit: maf
        // path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.sample}"
    def vep_cache_cmd = vep_data       ? "--vep-data $vep_data " : ""
    def VERSION       = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    if [ "$vep_data" ]; then
        # `type -p` is non-fatal here so a missing VEP produces a readable error below
        # rather than `dirname: missing operand` from an empty command substitution.
        VEP_BIN=\$(type -p vep || true)
        if [ -z "\$VEP_BIN" ]; then
            echo "ERROR: vep_data is set but no 'vep' executable was found on PATH inside the container." >&2
            echo "       container_vcf2maf = ${params.container_vcf2maf}" >&2
            echo "       vcf2maf.pl requires ensembl-vep. Point params.container_vcf2maf at an image" >&2
            echo "       that bundles VEP (e.g. vcf2maf_ensembl-vep), or unset params.vep_data." >&2
            exit 1
        fi
        VEP_CMD="--vep-path \$(dirname \$VEP_BIN)"
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_CMD=""
        VEP_VERSION=""
    fi

    # Handle compressed VCF files
    if [[ $vcf == *.gz ]]; then
        tmp=\$(mktemp --suffix=.vcf)
        rm -f "\$tmp"
        gunzip -c "$vcf" > "\$tmp"
        INPUT_VCF="\$tmp"
    else
        INPUT_VCF="$vcf"
    fi

    cat \$INPUT_VCF | grep "#" > tmp.${meta.sample}.somatic.vcf
    cat \$INPUT_VCF | grep PASS >> tmp.${meta.sample}.somatic.vcf

    VEP_PARAMS_RAW="${params.vep_params ?: ''}"

    # detect test mode
    if echo "\$VEP_PARAMS_RAW" | grep -q -- "--test"; then
        echo "VEP test mode detected; removing --test and skipping ID extraction"
        # remove the literal '--test' from the parameter string
        VEP_PARAMS=\$(echo "\$VEP_PARAMS_RAW" | sed 's/--test//g')
        ID_ARGS=""
    else
        # not test mode → extract IDs
        TMP_NORMAL_ID=\$(grep "^#CHROM" \$INPUT_VCF | awk '{print \$10}')
        TMP_TUMOR_ID=\$(grep "^#CHROM" \$INPUT_VCF | awk '{print \$11}')

        ID_ARGS=""
        [ -n "\$TMP_TUMOR_ID" ] && ID_ARGS="\$ID_ARGS --tumor-id \$TMP_TUMOR_ID"
        [ -n "\$TMP_NORMAL_ID" ] && ID_ARGS="\$ID_ARGS --normal-id \$TMP_NORMAL_ID"

        VEP_PARAMS="\$VEP_PARAMS_RAW"
    fi

    vcf2maf.pl \\
        \$ID_ARGS \\
        \$VEP_PARAMS \\
        \$VEP_CMD \\
        $vep_cache_cmd \\
        --ref-fasta $fasta \\
        --input-vcf tmp.${meta.sample}.somatic.vcf \\
        --output-maf tmp.${meta.sample}.maf

    head -2 tmp.${meta.sample}.maf > ${meta.sample}.maf
    tail -n +3 tmp.${meta.sample}.maf | awk -v col16="${meta.sample}" 'BEGIN {FS=OFS="\\t"} {\$16=col16; print}' >> ${meta.sample}.maf
    """
}
