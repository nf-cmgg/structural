//
// Run Gridss
//

include { GRIDSS_GRIDSS      } from '../../../modules/nf-core/gridss/gridss/main'

include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_GRIDSS {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file
        bwa_index               // channel: [mandatory] [ index ] => The BWA MEM index

    main:

    ch_versions     = Channel.empty()

    GRIDSS_GRIDSS(
        crams.map {meta, cram, crai -> [meta, cram, []]},
        fasta.map {[[],it]},
        fasta_fai.map {[[],it]},
        bwa_index
    )

    ch_versions = ch_versions.mix(GRIDSS_GRIDSS.out.versions)

    TABIX_TABIX(
        GRIDSS_GRIDSS.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    GRIDSS_GRIDSS.out.vcf
        .join(TABIX_TABIX.out.tbi)
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
