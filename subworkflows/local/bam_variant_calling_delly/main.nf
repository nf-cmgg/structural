//
// Run Delly
//

include { REVERSE_BED        } from '../../../modules/local/reversebed/main'

include { DELLY_CALL         } from '../../../modules/nf-core/delly/call/main'
include { BCFTOOLS_CONCAT    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT      } from '../../../modules/nf-core/bcftools/sort/main'
include { BCFTOOLS_CONVERT   } from '../../../modules/nf-core/bcftools/convert/main'
include { BEDTOOLS_SPLIT     } from '../../../modules/nf-core/bedtools/split/main'
include { TABIX_TABIX        } from '../../../modules/nf-core/tabix/tabix/main'

workflow BAM_VARIANT_CALLING_DELLY {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    scatter_count   = params.delly_scatter_count

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

    if(scatter_count > 1){
        BEDTOOLS_SPLIT(
            single_beds
        )

        BEDTOOLS_SPLIT.out.beds
            .transpose()
            .map(
                { meta, bed ->
                    new_meta = meta.clone()
                    new_meta.id = bed.baseName
                    [ new_meta, bed ]
                }
            )
            .set { split_beds }

        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

    } else {
        split_beds = single_beds
    }

    split_beds.dump(tag: 'split_beds', pretty: true)

    //
    // Reverse the BED file (It will only contain the regions that aren't of interest now)
    //

    REVERSE_BED(
        split_beds,
        fasta_fai
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
                        new_meta = meta.clone()
                        new_meta.id = meta.sample
                        [ new_meta, bed ]
                    }
                )
        , by:0)
        .map(
            { meta, cram, crai, bed ->
                new_meta = meta.clone()
                new_meta.id = bed.baseName.replace("_reversed","")
                [ new_meta, cram, crai, bed ]
            }
        )
        .dump(tag: 'delly_input', pretty: true)
        .set { delly_input }

    DELLY_CALL(
        delly_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(DELLY_CALL.out.versions)

    //
    // Concatenate the BCF files if the scatter count is bigger than 1 and convert the BCF to VCF
    //

    // TODO: fix the group size for when the BED file couldn't be split enough to reach the scatter count (see germline pipeline)

    DELLY_CALL.out.bcf
        .combine(DELLY_CALL.out.csi, by:0)
        .map(
            { meta, bcf, csi ->
                new_meta = meta.clone()
                new_meta.id = meta.sample
                [ new_meta, bcf, csi ]
            }
        )
        .groupTuple(size: scatter_count, remainder: true)
        .dump(tag: 'bcftools_input', pretty: true)
        .set { bcftools_input }

    if(scatter_count > 1){

        BCFTOOLS_CONCAT(
            bcftools_input
        )

        ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

        BCFTOOLS_CONCAT.out.vcf.set { sort_input }

    } else {

        BCFTOOLS_CONVERT(
            bcftools_input,
            [],
            fasta
        )

        ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions)

        BCFTOOLS_CONVERT.out.vcf_gz.set { sort_input }

    }

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
                new_meta = meta.clone()
                new_meta.caller = "delly"
                [ new_meta, vcf, tbi ]
            }
        )
        .dump(tag: 'delly_vcfs', pretty: true)
        .set { delly_vcfs }

    emit:
    delly_vcfs
    versions = ch_versions
}
