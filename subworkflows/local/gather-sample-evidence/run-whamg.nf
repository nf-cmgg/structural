//
// Run Whamg
//
include { WHAMG                     } from '../../../modules/nf-core/whamg/main'
include { BCFTOOLS_CONCAT           } from '../../../modules/nf-core/bcftools/concat/main'
include { SAMTOOLS_CONVERT          } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_BGZIPTABIX          } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_TABIX               } from '../../../modules/nf-core/tabix/tabix/main'

workflow RUN_WHAMG {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions      = Channel.empty()
    include_bed_file = params.whamg_include_bed_file

    //
    // Convert the CRAMs to BAMs
    //

    SAMTOOLS_CONVERT(
        crams,
        fasta,
        fasta_fai
    )

    bams        = SAMTOOLS_CONVERT.out.alignment_index
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //
    // Calling variants using Whamg (BED file support isn't present in Whamg, the `include_bed_file` option splits the BED file and runs the caller for every region)
    //

    if(include_bed_file){
        stringified_beds = beds.map({ meta, bed, bed_gz, bed_gz_tbi -> [ meta, bed ]})
                               .splitText()

        whamg_input = bams.combine(stringified_beds, by: 0)
                                .map({ meta, bam, bai, region ->
                                    new_meta = meta.clone()
                                    new_meta.region = region.replace("\n","").replaceFirst("\t",":").replace("\t","-")
                                    [ new_meta, bam, bai ]
                                })
    } else {
        whamg_input = bams
    }

    WHAMG(
        whamg_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(WHAMG.out.versions)

    //
    // Gzip and index the resulting VCF
    //

    TABIX_BGZIPTABIX(
        WHAMG.out.vcf
    )

    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    //
    // Concatenate and re-index the VCFs if the `include_bed_file` option was used
    //

    if(include_bed_file){
        concat_input = TABIX_BGZIPTABIX.out.gz_tbi.map({ meta, vcf, tbi ->
                            new_meta = meta.clone()
                            new_meta.remove('region')
                            [ new_meta, vcf, tbi ]
                        }).groupTuple()

        BCFTOOLS_CONCAT(
            concat_input
        )

        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

        TABIX_TABIX(
            BCFTOOLS_CONCAT.out.vcf
        )

        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

        whamg_vcfs = BCFTOOLS_CONCAT.out.vcf.combine(TABIX_TABIX.out.tbi, by: 0)
    } else {
        whamg_vcfs = TABIX_BGZIPTABIX.out.gz_tbi
    }

    whamg_vcfs = whamg_vcfs.map({ meta, vcf, tbi ->
        new_meta = meta.clone()
        new_meta.caller = "whamg"
        [ new_meta, vcf, tbi ]
    })

    emit:
    whamg_vcfs
    versions = ch_versions
}
