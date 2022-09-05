//
// Gather Sample Evidence
//
include { DELLY_CALL         } from '../../modules/nf-core/modules/delly/call/main'
include { REVERSE_BED        } from '../../modules/local/reversebed/main'
include { BCFTOOLS_CONCAT    } from '../../modules/nf-core/modules/bcftools/concat/main'
include { TABIX_TABIX        } from '../../modules/nf-core/modules/tabix/tabix/main'

workflow RUN_DELLY {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    sv_types = params.sv_types ? params.sv_types.split(",") : ["DEL","DUP","INV"]

    REVERSE_BED(
        beds.map({ meta, bed, bed_gz, bed_gz_tbi -> [ meta, bed ]}),
        fasta_fai
    )

    delly_input = crams.combine(REVERSE_BED.out.bed, by:0)
                       .combine(sv_types)
                       .map({ meta, cram, crai, bed, type ->
                            new_meta = meta.clone()
                            new_meta.sv_type = type
                            [ new_meta, cram, crai, bed ]
                        })

    DELLY_CALL(
        delly_input,
        fasta,
        fasta_fai
    )

    concat_input = DELLY_CALL.out.bcf.combine(DELLY_CALL.out.csi, by:0)
                            .map({ meta, bcf, csi ->
                                new_meta = [:]
                                new_meta.id = meta.id
                                [ new_meta, bcf, csi ]
                            }).groupTuple()

    BCFTOOLS_CONCAT(
        concat_input
    )

    TABIX_TABIX(
        BCFTOOLS_CONCAT.out.vcf
    )

    delly_vcfs = BCFTOOLS_CONCAT.out.vcf.combine(TABIX_TABIX.out.tbi, by:0)
                                .map({ meta, vcf, tbi -> 
                                    new_meta = meta.clone()
                                    new_meta.caller = "delly"
                                    [ new_meta, vcf, tbi ]
                                })

    emit:
    delly_vcfs
}