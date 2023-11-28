//
// Estimate repeat size with Expanionhunter
//

include { EXPANSIONHUNTER   } from '../../../modules/nf-core/expansionhunter/main'
include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'
include { BCFTOOLS_REHEADER } from '../../../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_REPEAT_ESTIMATION_EXPANSIONHUNTER {
    take:
        ch_crams        // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        ch_fasta        // channel: [mandatory] [ meta, fasta ] => The fasta reference file
        ch_fai          // channel: [mandatory] [ meta, fai ] => The index of the fasta reference file
        ch_catalogue    // channel: [mandatory] [ meta, catalogue ] => The expansionhunter catalogue

    main:

    ch_versions = Channel.empty()

    ch_crams
        .map { meta, cram, crai ->
            def new_meta = meta + ["variant_type":"repeats"]
            [ new_meta, cram, crai ]
        }
        .set { ch_expansionhunter_input }

    EXPANSIONHUNTER(
        ch_expansionhunter_input,
        ch_fasta,
        ch_fai,
        ch_catalogue
    )
    ch_versions = ch_versions.mix(EXPANSIONHUNTER.out.versions.first())

    Channel.of(
        '##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">',
        '##INFO=<ID=RU,Number=1,Type=String,Description="The repeat unit of the STR">',
        '##INFO=<ID=REPREF,Number=1,Type=Integer,Description="How many repeats are in the reference">',
        '##INFO=<ID=RL,Number=1,Type=Integer,Description="The repeat length in basepairs">',
        '##INFO=<ID=REPID,Number=1,Type=String,Description="The ID of the repeat">',
        '##INFO=<ID=VARID,Number=1,Type=String,Description="Variant identifier as specified in the variant catalog">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
        '##FORMAT=<ID=LC,Number=1,Type=Float,Description="Locus coverage">',
        '##FORMAT=<ID=SO,Number=1,Type=String,Description="The type of reads used to estimate the repeat expansions">',
        '##FORMAT=<ID=REPCN,Number=1,Type=String,Description="The mean copy number in this repeat expansion">',
        '##FORMAT=<ID=REPCI,Number=1,Type=String,Description="The confidence interval for the size of the expanded allele">',
        '##FORMAT=<ID=ADSP,Number=1,Type=String,Description="The amount of spanning reads consistent with the repeat allele">',
        '##FORMAT=<ID=ADFL,Number=1,Type=String,Description="The amount of flanking reads that overlap the repeat allele">',
        '##FORMAT=<ID=ADIR,Number=1,Type=String,Description="The amount of in-repeat reads that are consistent with the repeat allele">'
        )
        .collectFile(name:'header.txt', newLine:true)
        .set { ch_ref_header }

    EXPANSIONHUNTER.out.vcf
        .combine(ch_ref_header)
        .map { meta, vcf, ref_header ->
            [ meta, vcf, [], [], [], ref_header]
        }
        .set { ch_annotate_input }

    BCFTOOLS_ANNOTATE(
        ch_annotate_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())

    Channel.fromPath("${projectDir}/assets/header.txt")
        .collect()
        .set { ch_new_header }

    BCFTOOLS_ANNOTATE.out.vcf
        .combine(ch_new_header)
        .map { meta, vcf, header ->
            [ meta, vcf, header, [] ]
        }
        .set { ch_reheader_input }

    BCFTOOLS_REHEADER(
        ch_reheader_input,
        ch_fai
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())

    TABIX_TABIX(
        BCFTOOLS_REHEADER.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    BCFTOOLS_REHEADER.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_vcfs }

    emit:
    vcfs                = ch_vcfs    // channel: [ val(meta), path(vcf), path(tbi) ]

    versions            = ch_versions
}
