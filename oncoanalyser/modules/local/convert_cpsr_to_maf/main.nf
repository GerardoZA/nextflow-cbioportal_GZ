process CONVERT_CPSR_TO_MAF {
    publishDir { "${params.outdir}/${maf_meta.group}/${maf_meta.subject}" }, mode: 'copy'

    tag { maf_meta.sample }
    label "process_medium_memory"
    container params.container_r

    input:
        tuple val(join_key), val(maf_meta), path(som_dna_rna_maf), val(ger_meta), path(ger_dna_tsv_gz)

    output:
        tuple val(maf_meta), path("${maf_meta.sample}.somatic_rna_germline.maf")


    script:
    """
    zcat $ger_dna_tsv_gz > tmp.tsv
    # keep Pathogenic / Likely_Pathogenic / VUS, locating the classification column by name
    # (position differs between CPSR versions: \$52 is NULL_VARIANT in 2.3.x). Prefer the final
    # CLASSIFICATION (ClinVar-or-CPSR), fall back to CPSR_CLASSIFICATION; header-only if neither exists.
    awk -F"\\t" '
        NR==1 {
            for (i = 1; i <= NF; i++) {
                if (\$i == "CLASSIFICATION") c = i
                if (\$i == "CPSR_CLASSIFICATION") cc = i
            }
            if (!c) c = cc
            print; next
        }
        c { v = tolower(\$c); gsub(/ /, "_", v)
            if (v == "pathogenic" || v == "likely_pathogenic" || v == "vus") print }
    ' tmp.tsv > tmp.germline.cpsr.tsv

    rm tmp.tsv # to reduce size of work dir

    # writes somatic + germline rows to tmp.*.maf (do not overwrite it with the somatic-only input afterwards)
    gen_convert_cpsr_to_maf.R \
       tmp.germline.cpsr.tsv \
       $som_dna_rna_maf \
       tmp.${maf_meta.sample}.somatic_rna_germline.maf

    head -n2 tmp.${maf_meta.sample}.somatic_rna_germline.maf > ${maf_meta.sample}.somatic_rna_germline.maf
    awk -F'\t' 'NR>2{if(\$9!="Intron" && \$9!="IGR"){print \$0}}' tmp.${maf_meta.sample}.somatic_rna_germline.maf >> ${maf_meta.sample}.somatic_rna_germline.maf
    """

}
