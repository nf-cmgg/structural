//
// Call copy number variants
//

// Import subworkflows
include { BAM_VARIANT_CALLING_QDNASEQ                   } from '../bam_variant_calling_qdnaseq/main'

workflow BAM_CNV_CALLING {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_reference    // channel: [mandatory] [ meta, qdnaseq_reference ] => The reference for qDNAseq

    main:

    val_callers     = params.callers.tokenize(",").intersect(params.cnv_callers)

    ch_versions     = Channel.empty()
    ch_reports      = Channel.empty()
    ch_called_vcfs  = Channel.empty()

    BAM_VARIANT_CALLING_QDNASEQ(
        ch_crams,
        ch_fasta,
        ch_fai,
        ch_qdnaseq_reference
    )
    ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_QDNASEQ.out.versions)

    emit:
    vcfs                = Channel.empty()    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions            = ch_versions
    reports             = ch_reports
}
