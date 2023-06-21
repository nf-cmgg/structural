//
// Run QDNAseq
//

include { QDNASEQ           } from '../../../modules/local/qdnaseq/main'
include { SAMTOOLS_CONVERT  } from '../../../modules/nf-core/samtools/convert/main'

workflow BAM_VARIANT_CALLING_QDNASEQ {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_reference    // channel: [mandatory] [ meta, qdnaseq_reference ] => The reference to be used for qDNAseq

    main:

    ch_versions     = Channel.empty()

    SAMTOOLS_CONVERT(
        ch_crams,
        ch_fasta.map { it[1] },
        ch_fai.map { it[1] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions.first())

    QDNASEQ(
        SAMTOOLS_CONVERT.out.alignment_index,
        ch_qdnaseq_reference
    )
    ch_versions = ch_versions.mix(QDNASEQ.out.versions.first())

    emit:
    qdnaseq_beds  = QDNASEQ.out.bed  // channel: [ val(meta), path(bed) ]

    versions    = ch_versions
}
