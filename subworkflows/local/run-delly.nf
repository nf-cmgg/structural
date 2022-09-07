//
// Gather Sample Evidence
//
include { DELLY_CALL         } from '../../modules/nf-core/modules/delly/call/main'
include { REVERSE_BED        } from '../../modules/local/reversebed/main'
include { BCFTOOLS_CONCAT    } from '../../modules/nf-core/modules/bcftools/concat/main'
include { BEDTOOLS_SPLIT     } from '../../modules/nf-core/modules/bedtools/split/main'
include { TABIX_TABIX        } from '../../modules/nf-core/modules/tabix/tabix/main'

workflow RUN_DELLY {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    sv_types        = params.sv_types.split(",")
    scatter_count   = params.delly_scatter_count

    single_beds     = beds.map({ meta, bed, bed_gz, bed_gz_tbi -> [ meta, bed ]})


    if(scatter_count > 1){
        BEDTOOLS_SPLIT(
            single_beds,
            scatter_count
        )

        split_beds = BEDTOOLS_SPLIT.out.beds.transpose().map({ meta, bed ->
                            new_meta = meta.clone()
                            new_meta.id = bed.baseName
                            [ new_meta, bed ]
                        })
        ch_versions = ch_versions.mix(BEDTOOLS_SPLIT.out.versions)

    } else {
        split_beds = single_beds
    }

    REVERSE_BED(
        split_beds,
        fasta_fai
    )
    ch_versions = ch_versions.mix(REVERSE_BED.out.versions)

    delly_input = crams.combine(REVERSE_BED.out.bed.map({ meta, bed -> 
                                new_meta = meta.clone()
                                new_meta.id = meta.sample
                                [ new_meta, bed ]
                            }), by:0)
                       .combine(sv_types)
                       .map({ meta, cram, crai, bed, type ->
                            new_meta = meta.clone()
                            new_meta.sv_type = type
                            new_meta.id = bed.baseName.replace("_reversed","")
                            [ new_meta, cram, crai, bed ]
                        })

    DELLY_CALL(
        delly_input,
        fasta,
        fasta_fai
    )
    ch_versions = ch_versions.mix(DELLY_CALL.out.versions)

    concat_input = DELLY_CALL.out.bcf.combine(DELLY_CALL.out.csi, by:0)
                            .map({ meta, bcf, csi ->
                                new_meta = meta.clone()
                                new_meta.remove('sv_type')
                                [ new_meta, bcf, csi ]
                            }).groupTuple()

    BCFTOOLS_CONCAT(
        concat_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    TABIX_TABIX(
        BCFTOOLS_CONCAT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    delly_vcfs = BCFTOOLS_CONCAT.out.vcf.combine(TABIX_TABIX.out.tbi, by:0)
                                .map({ meta, vcf, tbi -> 
                                    new_meta = meta.clone()
                                    new_meta.caller = "delly"
                                    [ new_meta, vcf, tbi ]
                                })

    emit:
    delly_vcfs
    versions = ch_versions
}