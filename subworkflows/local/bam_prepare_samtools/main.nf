//
// Prepare the CRAM files
//

// Import modules
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_PREPARE_SAMTOOLS {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file

    main:

    def ch_versions     = channel.empty()

    def ch_merge_input = ch_crams
        .groupTuple() // no size needed here as no process has been run before this
        .branch { meta, cram, crai ->
            multiple: cram.size() > 1
                return [ meta, cram ]
            single: cram.size() == 1
                return [ meta, cram[0], crai[0] ]
        }

    SAMTOOLS_MERGE(
        ch_merge_input.multiple,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    def ch_index_input = ch_merge_input.single
        .mix(SAMTOOLS_MERGE.out.cram)
        .branch { meta, cram, crai=[] ->
            index: crai
            no_index: !crai
                return [ meta, cram ]
        }

    SAMTOOLS_INDEX(
        ch_index_input.no_index
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    def ch_crams_ready = ch_index_input.no_index
        .join(SAMTOOLS_INDEX.out.crai, failOnMismatch:true, failOnDuplicate:true)
        .mix(ch_index_input.index)

    emit:
    crams    = ch_crams_ready // channel: [ val(meta), path(cram), path(crai) ]

    versions = ch_versions
}
