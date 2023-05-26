//
// Run Delly
//

include { DELLY_CALL         } from '../../../modules/nf-core/delly/call/main'
include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
        ch_crams    // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta    // channel: [mandatory] [ fasta ] => The fasta reference file
        ch_fai      // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Calling variants using Delly
    //

    ch_crams
        .map { meta, cram, crai, ->
            [ meta, cram, crai, [], [], [] ]
        }
        .dump(tag: 'delly_input', pretty: true)
        .set { ch_delly_input }

    DELLY_CALL(
        ch_delly_input,
        ch_fasta,
        ch_fai
    )

    ch_versions = ch_versions.mix(DELLY_CALL.out.versions.first())

    DELLY_CALL.out.bcf
        .join(DELLY_CALL.out.csi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi ->
            new_meta = meta + [caller:"delly"]
            [ new_meta, vcf, tbi ]
        }
        .dump(tag: 'delly_vcfs', pretty: true)
        .set { ch_delly_vcfs }

    emit:
    delly_vcfs  = ch_delly_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}
