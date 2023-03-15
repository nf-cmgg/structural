//
// Run Whamg
//
include { WHAMG                       } from '../../../modules/nf-core/whamg/main'
include { BCFTOOLS_CONCAT             } from '../../../modules/nf-core/bcftools/concat/main'
include { BEDTOOLS_MERGE              } from '../../../modules/nf-core/bedtools/merge/main'
include { SAMTOOLS_CONVERT            } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX as TABIX_CONCAT } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_WHAMG  } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_WHAMG {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions      = Channel.empty()

    //
    // Convert the CRAMs to BAMs
    //

    SAMTOOLS_CONVERT(
        crams,
        fasta,
        fasta_fai
    )

    SAMTOOLS_CONVERT.out.alignment_index.set { bams }
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //
    // Calling variants using Whamg
    //

    WHAMG(
        bams,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(WHAMG.out.versions)

    WHAMG.out.vcf
        .join(WHAMG.out.tbi, failOnMismatch:true, failOnDuplicate:true)
        .set { whamg_vcfs }

    whamg_vcfs
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.findAll {true} + [caller:"whamg"]
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'whamg_vcfs', pretty: true)
        .set { whamg_vcfs }

    emit:
    whamg_vcfs
    versions = ch_versions
}
