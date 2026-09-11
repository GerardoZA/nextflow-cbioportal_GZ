include { VCF2MAF              } from '../../../modules/local/vcf2maf'
include { INTEGRATE_RNA_VARIANTS } from '../../../modules/local/integrate_rna_variants'
include { PCGR                 } from '../../../modules/local/pcgr'
include { CONVERT_CPSR_TO_MAF  } from '../../../modules/local/convert_cpsr_to_maf'
include { DOWNLOAD_VEP_TEST    } from '../../../modules/local/download_vep_test'
include { DOWNLOAD_PCGR        } from '../../../modules/local/download_pcgr'
include { BCFTOOLS_INDEX       } from '../../../modules/nf-core/bcftools/index'
include { FILTER_GERMLINE_DNA  } from '../../../modules/local/filter_germline_dna'

// Case-insensitive: params values set on the command line are immutable, so "--type DNA-only"
// reaches here unnormalised and an exact-match comparison would silently be false.
def isDnaOnly() {
    return (params.type ?: 'both').toString().toLowerCase() == 'dna-only'
}

workflow GENOMIC_MUTATIONS {
    take:
        ger_dna_vcf // tuple (meta, pave.germline.vcf.gz) — pave germline VCF
        som_dna_vcf // tuple (meta, pave.somatic.vcf.gz)  — pave somatic VCF
        som_rna_vcf // tuple (meta, sage.append.vcf.gz)   — SAGE RNA-append VCF
        fasta
        vep_data
        pcgr_data
        needs_vep
        needs_pcgr

    main:
        ch_vep_data  = needs_vep  ? DOWNLOAD_VEP_TEST().cache_dir.first() : vep_data.first()
        ch_pcgr_data = needs_pcgr ? DOWNLOAD_PCGR().data_dir.first()      : pcgr_data.first()

        ger_dna_filtered = FILTER_GERMLINE_DNA(ger_dna_vcf)

        ger_dna_index = BCFTOOLS_INDEX(ger_dna_filtered).tbi

        ger_dna_vcf_with_index = ger_dna_filtered
            .join(ger_dna_index)
            .map { meta, filepath, index -> tuple(meta, filepath, index) }

        ger_dna_tsv = PCGR(
            ger_dna_vcf_with_index,
            ch_vep_data,
            ch_pcgr_data
        )

        // add germinal_sample to somatic meta so vcf2maf can set the normal sample column
        som_dna_vcf_input = som_dna_vcf
            .map { meta, vcf -> tuple(meta.subject, meta, vcf) }
            .join(
                ger_dna_vcf.map { meta, vcf -> tuple(meta.subject, meta) }
            )
            .map { subject, som_meta, som_vcf, ger_meta ->
                def meta = som_meta + [germinal_sample: ger_meta.sample]
                return tuple(meta, som_vcf)
            }

        VCF2MAF(
            som_dna_vcf_input,
            fasta,
            ch_vep_data
        )

        som_dna_maf = VCF2MAF.out.maf.map { meta, maf ->
            return tuple(meta, maf)
        }

        // dna-only mode: no RNA VCF exists for any subject, so skip the RNA integration
        // join entirely and pass the DNA MAF straight through (an inner join against an
        // always-empty RNA channel would otherwise silently drop every subject's mutations)
        if (isDnaOnly()) {
            som_dna_rna_maf = som_dna_maf
        } else {
            // join RNA-append VCF with somatic DNA MAF on subject
            som_rna_dna_tuple = som_rna_vcf
                .map { meta, file -> tuple(meta.subject, meta, file) }
                .join(
                    som_dna_maf.map { meta, file -> tuple(meta.subject, meta, file) }
                )

            som_dna_rna_maf = INTEGRATE_RNA_VARIANTS(som_rna_dna_tuple)
        }

        som_dna_maf_tsv = som_dna_rna_maf
            .map { meta, file -> return tuple(meta.subject, meta, file) }
            .join(
                ger_dna_tsv.map { meta, file -> return tuple(meta.subject, meta, file) }
            )

        cbioportal_genomic_mutation_files = CONVERT_CPSR_TO_MAF(som_dna_maf_tsv)

    emit:
        out = cbioportal_genomic_mutation_files
}
