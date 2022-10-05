//
// Run Manta
//
include { MANTA_GERMLINE         } from '../../../modules/nf-core/manta/germline/main'
include { MANTA_CONVERTINVERSION } from '../../../modules/nf-core/manta/convertinversion/main'
include { BCFTOOLS_REHEADER      } from '../../../modules/nf-core/bcftools/reheader/main'

workflow RUN_MANTA {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file

    main:

    ch_versions     = Channel.empty()

    gzipped_beds = beds.map({ meta, bed, bed_gz, bed_gz_tbi -> [ meta, bed_gz, bed_gz_tbi ]})

    manta_input = crams.combine(gzipped_beds, by: 0)

    MANTA_GERMLINE(
        manta_input,
        fasta,
        fasta_fai
    )

    ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions)

    MANTA_CONVERTINVERSION(
        MANTA_GERMLINE.out.diploid_sv_vcf,
        fasta
    )

    ch_versions = ch_versions.mix(MANTA_CONVERTINVERSION.out.versions)

    BCFTOOLS_REHEADER(
        MANTA_CONVERTINVERSION.out.vcf,
        [],
        []
    )

    manta_vcfs = BCFTOOLS_REHEADER.out.vcf.combine(MANTA_CONVERTINVERSION.out.tbi, by:0)
                                               .map({ meta, vcf, tbi ->
                                                   new_meta = meta.clone()
                                                   new_meta.caller = "manta"
                                                   [ new_meta, vcf, tbi ]
                                               })

    emit:
    manta_vcfs
    versions = ch_versions
}
