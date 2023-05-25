//
// Run Expansionhunter
//

include { EXPANSIONHUNTER      } from '../../../modules/nf-core/expansionhunter/main'
include { TABIX_TABIX          } from '../../../modules/nf-core/tabix/tabix/main'
include { NGSBITS_SAMPLEGENDER } from '../../../modules/nf-core/ngsbits/samplegender/main'

workflow BAM_REPEATS_ESTIMATION_EXPANSIONHUNTER {
    take:
        ch_crams                // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta                // channel: [mandatory] [ meta2, fasta ] => The fasta reference file
        ch_fai                  // channel: [mandatory] [ meta3, fai ] => The index of the fasta reference file
        ch_variant_catalogue    // channel: [mandatory] [ meta4, variant_catalogue ] => The variant catalogue for expansionhunter

    main:

    ch_versions     = Channel.empty()

    NGSBITS_SAMPLEGENDER(
        ch_crams,
        ch_fasta,
        ch_fai,
        "xy"
    )
    ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions.first())

    NGSBITS_SAMPLEGENDER.out.tsv
        .join(ch_crams, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, tsv, cram, crai ->
            new_meta = meta + [gender:define_gender(tsv)]
            [ new_meta, cram, crai ]
        }
        .set { ch_expansionhunter_input }

    EXPANSIONHUNTER(
        ch_expansionhunter_input,
        ch_fasta,
        ch_fai,
        ch_variant_catalogue
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())

    TABIX_TABIX(
        EXPANSIONHUNTER.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    EXPANSIONHUNTER.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_expansionhunter_vcfs }

    emit:
    expansionhunter_vcfs = ch_expansionhunter_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]

    versions    = ch_versions
}

def define_gender(tsv) {
    gender = tsv.splitCsv(sep:"\t", header:true)[0].gender
    result = gender in ["male", "female"] ? gender : null
    return result
}