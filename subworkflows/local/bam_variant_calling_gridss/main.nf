//
// Run Gridss
//

include { SIMPLE_EVENT_ANNOTATION   } from '../../../modules/local/gridss/simple_event_annotation/main'

include { GRIDSS_GRIDSS             } from '../../../modules/nf-core/gridss/gridss/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_GRIDSS {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fai               // channel: [mandatory] [ fai ] => The index of the fasta reference file
        bwa_index               // channel: [mandatory] [ index ] => The BWA MEM index

    main:

    ch_versions     = Channel.empty()

    GRIDSS_GRIDSS(
        crams.map {meta, cram, crai -> [meta, cram, []]},
        fasta.map {[[],it]},
        fai.map {[[],it]},
        bwa_index
    )
    ch_versions = ch_versions.mix(GRIDSS_GRIDSS.out.versions)

    SIMPLE_EVENT_ANNOTATION(
        GRIDSS_GRIDSS.out.vcf
    )
    ch_versions = ch_versions.mix(SIMPLE_EVENT_ANNOTATION.out.versions)

    TABIX_TABIX(
        SIMPLE_EVENT_ANNOTATION.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    SIMPLE_EVENT_ANNOTATION.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .map(
            { meta, vcf, tbi ->
                new_meta = meta + [caller:"gridss"]
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'gridss_vcfs', pretty: true)
        .set { gridss_vcfs }

    emit:
    gridss_vcfs
    versions = ch_versions
}
