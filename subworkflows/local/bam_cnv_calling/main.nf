//
// Call copy number variants
//

// Import modules
include { JASMINESV     } from '../../../modules/nf-core/jasminesv/main'
include { TABIX_TABIX   } from '../../../modules/nf-core/tabix/tabix/main'

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
        ch_called_vcfs = ch_called_vcfs.mix(BAM_VARIANT_CALLING_QDNASEQ.out.vcf)
    }

    if("wisecondorx" in val_callers) {
        BAM_VARIANT_CALLING_WISECONDORX(
            ch_crams,
            ch_fasta,
            ch_fai,
            ch_wisecondorx_reference,
            ch_blacklist
        )
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_WISECONDORX.out.versions)
        ch_called_vcfs = ch_called_vcfs.mix(BAM_VARIANT_CALLING_WISECONDORX.out.vcf)
    }

    if(val_callers.size() > 1) {
        ch_called_vcfs
            .groupTuple(size:val_callers.size())
            .map { meta, vcfs ->
                [ meta, vcfs, [], []]
            }
            .set { ch_jasmine_input }

        JASMINESV(
            ch_jasmine_input,
            ch_fasta.map { it[1] },
            ch_fai.map { it[1] },
            []
        )
        ch_versions = ch_versions.mix(JASMINESV.out.versions.first())

        JASMINESV.out.vcf
            .set { ch_vcfs }
    } else {
        ch_called_vcfs
            .set { ch_vcfs }
    }

    TABIX_TABIX(
        ch_vcfs
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_vcfs
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { ch_vcfs_out }

    emit:
    versions            = ch_versions
    reports             = ch_reports
    vcfs                = ch_vcfs_out
}
