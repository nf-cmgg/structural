//
// Call copy number variants
//

// Import subworkflows
include { BAM_VARIANT_CALLING_QDNASEQ       } from '../bam_variant_calling_qdnaseq/main'
include { BAM_VARIANT_CALLING_WISECONDORX   } from '../bam_variant_calling_wisecondorx/main'

workflow BAM_CNV_CALLING {
    take:
        ch_crams                    // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_male             // channel: [mandatory] [ meta, qdnaseq_reference ] => The male reference for qDNAseq
        ch_qdnaseq_female           // channel: [mandatory] [ meta, qdnaseq_reference ] => The female reference for qDNAseq
        ch_wisecondorx_reference    // channel: [mandatory] [ meta, wisecondorx_reference ] => The reference for WisecondorX
        ch_blacklist                // channel: [optional]  [ meta, bed ] => The blacklist regions to be excluded from the Wisecondorx analysis

    main:

    val_callers     = params.callers.tokenize(",").intersect(GlobalVariables.cnvCallers)

    ch_versions     = Channel.empty()
    ch_reports      = Channel.empty()
    ch_called_vcfs  = Channel.empty()

    if("qdnaseq" in val_callers) {
        BAM_VARIANT_CALLING_QDNASEQ(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_qdnaseq_male,
            ch_qdnaseq_female
        )
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_QDNASEQ.out.versions)
    }

    if("wisecondorx" in val_callers) {
        BAM_VARIANT_CALLING_WISECONDORX(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_wisecondorx_reference,
            ch_blacklist
        )
    }

    emit:
    versions            = ch_versions
    reports             = ch_reports
}
