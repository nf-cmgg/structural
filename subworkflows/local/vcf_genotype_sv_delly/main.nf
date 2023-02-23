//
// Genotype using DELLY
//

include { DELLY_CALL  } from '../../../modules/nf-core/delly/call/main'
include { TABIX_TABIX } from '../../../modules/nf-core/tabix/tabix/main'

workflow VCF_GENOTYPE_SV_DELLY {
    take:
        vcfs                    // channel: [mandatory] [ meta, vcf, tbi ] VCFs containing the called structural variants
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The CRAM files used to create the VCF files
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    crams.
        join(vcfs)
        .map { it + [[]] }
        .set { delly_input }

    DELLY_CALL(
        delly_input,
        fasta,
        fasta_fai
    )
    ch_versions = ch_versions.mix(DELLY_CALL.out.versions)

    TABIX_TABIX(
        DELLY_CALL.out.bcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    DELLY_CALL.out.bcf
        .join(TABIX_TABIX.out.tbi)
        .view()
        .set { genotyped_vcfs }

    emit:
    genotyped_vcfs
    versions = ch_versions
}
