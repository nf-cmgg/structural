//
// Call copy number variants
//

// Import modules
include { TABIX_TABIX                       } from '../../../modules/nf-core/tabix/tabix/main'

// Import subworkflows
include { BAM_VARIANT_CALLING_QDNASEQ       } from '../bam_variant_calling_qdnaseq/main'
include { BAM_VARIANT_CALLING_WISECONDORX   } from '../bam_variant_calling_wisecondorx/main'
include { VCF_MERGE_CALLERS_JASMINE         } from '../vcf_merge_callers_jasmine/main'

workflow BAM_CNV_CALLING {
    take:
        ch_crams                    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                    // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                      // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_qdnaseq_male             // channel: [mandatory] [ meta, qdnaseq_reference ] => The male reference for qDNAseq
        ch_qdnaseq_female           // channel: [mandatory] [ meta, qdnaseq_reference ] => The female reference for qDNAseq
        ch_wisecondorx_reference    // channel: [mandatory] [ meta, wisecondorx_reference ] => The reference for WisecondorX
        ch_blacklist                // channel: [optional]  [ meta, bed ] => The blacklist regions to be excluded from the Wisecondorx analysis
        val_callers                 // value:   [mandatory] => List of all CNV callers to use

    main:

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
        VCF_MERGE_CALLERS_JASMINE(
            ch_called_vcfs.map { it + [[]] },
            ch_fasta,
            ch_fai,
            val_callers,
            "cnv"
        )
        ch_versions = ch_versions.mix(VCF_MERGE_CALLERS_JASMINE.out.versions)
        VCF_MERGE_CALLERS_JASMINE.out.vcfs
            .set { ch_merged_vcfs }

    } else {
        TABIX_TABIX(
            ch_called_vcfs
        )
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

        ch_called_vcfs
            .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
            .map { meta, vcf, tbi ->
                def new_meta = meta - meta.subMap("caller") + [variant_type:"cnv"]
                [ new_meta, vcf, tbi ]
            }
            .set { ch_merged_vcfs }
    }


    emit:
    versions            = ch_versions
    reports             = ch_reports
    vcfs                = ch_merged_vcfs
}
