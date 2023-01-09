//
// Run Whamg
//
include { WHAMG                       } from '../../../modules/nf-core/whamg/main'
include { BCFTOOLS_CONCAT             } from '../../../modules/nf-core/bcftools/concat/main'
include { BEDTOOLS_MERGE              } from '../../../modules/nf-core/bedtools/merge/main'
include { SAMTOOLS_CONVERT            } from '../../../modules/nf-core/samtools/convert/main'
include { TABIX_TABIX as TABIX_CONCAT } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_WHAMG  } from '../../../modules/nf-core/tabix/tabix/main'

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

    SAMTOOLS_CONVERT.out.alignment_index.set { bams }
    ch_versions = ch_versions.mix(SAMTOOLS_CONVERT.out.versions)

    //
    // Calling variants using Whamg (BED file support isn't present in Whamg, the `include_bed_file` option splits the BED file and runs the caller for every region)
    //

    if(include_bed_file){

        BEDTOOLS_MERGE{
            beds.map({ meta, bed, bed_gz, bed_gz_tbi -> [ meta, bed ]})
        }

        ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions)

        BEDTOOLS_MERGE.out.bed
            .splitText()
            .dump(tag:'stringified_beds', pretty:true)
            .tap { stringified_beds }
            .groupTuple()
            .map(
                { meta, regions ->
                    [ meta, regions.size() ]
                }
            )
            .set { region_count }

        bams
            .combine(stringified_beds, by: 0)
            .map(
                { meta, bam, bai, region ->
                    new_meta = meta.clone() + [region:region.replace("\n","").replaceFirst("\t",":").replace("\t","-")]
                    [ new_meta, bam, bai ]
                }
            )
            .set { whamg_input }

    } else {
        bams.set { whamg_input }
    }

    whamg_input.dump(tag: 'whamg_input', pretty: true)

    WHAMG(
        whamg_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(WHAMG.out.versions)

    //
    // Concatenate and re-index the VCFs if the `include_bed_file` option was used
    //

    if(include_bed_file){
        WHAMG.out.vcf
            .join(WHAMG.out.tbi)
            .map(
                { meta, vcf, tbi ->
                    new_meta = meta.clone()
                    new_meta.remove('region')
                    [ new_meta, vcf, tbi ]
                }
            )
            .combine(region_count, by:0)
            .map(
                { meta, vcf, tbi, region_count ->
                    [ groupKey(meta, region_count), vcf, tbi ]
                }
            )
            .groupTuple()
            .set { concat_input }

        BCFTOOLS_CONCAT(
            concat_input
        )

        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

        TABIX_CONCAT(
            BCFTOOLS_CONCAT.out.vcf
        )

        ch_versions = ch_versions.mix(TABIX_CONCAT.out.versions)

        BCFTOOLS_CONCAT.out.vcf
            .join(TABIX_CONCAT.out.tbi)
            .set { whamg_vcfs }
    } else {
        WHAMG.out.vcf
            .join(WHAMG.out.tbi)
            .set { whamg_vcfs }
    }

    whamg_vcfs
        .map(
            { meta, vcf, tbi ->
                new_meta = meta.findAll {true}[0] + [caller:"whamg"]
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'whamg_vcfs', pretty: true)
        .set { whamg_vcfs }

    emit:
    whamg_vcfs
    versions = ch_versions
}
