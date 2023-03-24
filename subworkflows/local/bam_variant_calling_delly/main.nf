//
// Run Delly
//

include { REVERSE_BED        } from '../../../modules/local/reversebed/main'
include { SCATTER_BEDS       } from '../../../modules/local/scatter_beds/main'

include { DELLY_CALL         } from '../../../modules/nf-core/delly/call/main'
include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fai               // channel: [mandatory] [ fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    //
    // Split the BED files if the scatter count is bigger than 1
    //

    beds
        .map(
            { meta, bed, bed_gz, bed_gz_tbi ->
                [ meta, bed ]
            }
        )
        .set { single_beds }

    SCATTER_BEDS(
        single_beds,
        params.delly_scatter_size
    )

    ch_versions = ch_versions.mix(SCATTER_BEDS.out.versions)

    SCATTER_BEDS.out.scatter
        .transpose()
        .map{ meta, bed ->
                new_meta = meta + [id:bed.baseName]
                [ new_meta, bed ]
        }
        .set { split_beds }

    split_beds.dump(tag: 'split_beds', pretty: true)

    //
    // Reverse the BED file (It will only contain the regions that aren't of interest now)
    //

    REVERSE_BED(
        split_beds,
        fai
    )

    ch_versions = ch_versions.mix(REVERSE_BED.out.versions)

    //
    // Calling variants using Delly
    //

    crams
        .combine(
            REVERSE_BED.out.bed
                .map(
                    { meta, bed ->
                        new_meta = meta + [id:meta.sample]
                        [ new_meta, bed ]
                    }
                )
        , by:0)
        .map(
            { meta, cram, crai, bed ->
                new_meta = meta + [id:bed.baseName.replace("_reversed","")]
                [ new_meta, cram, crai, [], [], bed ]
            }
        )
        .dump(tag: 'delly_input', pretty: true)
        .set { delly_input }

    DELLY_CALL(
        delly_input,
        fasta,
        fai
    )

    ch_versions = ch_versions.mix(DELLY_CALL.out.versions)

    //
    // Concatenate the BCF files if the scatter count is bigger than 1 and convert the BCF to VCF
    //

    DELLY_CALL.out.bcf
        .join(DELLY_CALL.out.csi, failOnMismatch:true, failOnDuplicate:true)
        .map { meta, vcf, csi ->
            new_meta = meta + [id:meta.sample]
            [ new_meta, vcf, csi ]
        }
        .combine(SCATTER_BEDS.out.scatter, by:0)
        .map(
            { meta, vcf, csi, beds ->
                count = beds instanceof ArrayList ? beds.size() : 1
                [ groupKey(meta, count), vcf, csi ]
            }
        )
        .groupTuple()
        .dump(tag: 'bcftools_input', pretty: true)
        .set { bcftools_input }


    BCFTOOLS_CONCAT(
        bcftools_input
    )

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    BCFTOOLS_CONCAT.out.vcf.set { sort_input }

    //
    // Index the VCF file
    //

    BCFTOOLS_SORT(
        sort_input
    )

    TABIX_TABIX(
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    BCFTOOLS_SORT.out.vcf
        .combine(TABIX_TABIX.out.tbi, by:0)
        .map(
            { meta, vcf, tbi ->
                new_meta = meta + [caller:"delly"]
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'delly_vcfs', pretty: true)
        .set { delly_vcfs }

    emit:
    delly_vcfs
    versions = ch_versions
}
