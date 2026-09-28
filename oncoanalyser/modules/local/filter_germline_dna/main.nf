process FILTER_GERMLINE_DNA {
    tag "$meta.sample"
    label 'process_low'

    container "${ workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
        tuple val(meta), path(ger_dna_vcf)

    output:
        tuple val(meta), path("*.vcf.gz")

    script:
    """
    zcat $ger_dna_vcf | grep "#" > tmp.${meta.sample}.vcf
    zgrep -v "^#" $ger_dna_vcf | grep PASS >> tmp.${meta.sample}.vcf

    # SAGE/PAVE germline VCFs leave GT as ./. — CPSR's vcf2tsvpy step drops ./. genotypes, so
    # every variant would be lost before classification. Call GT from the normal's FORMAT/AF
    # (existing genotypes are untouched). INFO/IMPACT (PAVE) clashes with a CPSR tag and aborts CPSR.
    bcftools view \\
        -s ^${meta.sample} \\
        -Ou \\
        tmp.${meta.sample}.vcf \\
    | bcftools annotate -x INFO/IMPACT -Ou \\
    | bcftools +setGT -Ou -- -t q -n c:1/1 -i 'GT="mis" && FMT/AF>=0.9' \\
    | bcftools +setGT -Oz -o ${meta.sample}.vcf.gz -- -t q -n c:0/1 -i 'GT="mis" && FMT/AF>0'
    """
}
