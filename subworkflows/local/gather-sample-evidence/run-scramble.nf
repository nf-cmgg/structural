//
// Run Scramble
//
include { SCRAMBLE_CLUSTERIDENTIFIER  } from '../../../modules/nf-core/scramble/clusteridentifier/main'
include { SCRAMBLE_CLUSTERANALYSIS    } from '../../../modules/nf-core/scramble/clusteranalysis/main'

workflow RUN_SCRAMBLE {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file

    main:

    ch_versions     = Channel.empty()

    SCRAMBLE_CLUSTERIDENTIFIER(
        crams.view(),
        fasta
    )

    ch_versions = ch_versions.mix(SCRAMBLE_CLUSTERIDENTIFIER.out.versions)

    SCRAMBLE_CLUSTERIDENTIFIER.out.clusters
        .filter(
            { meta, clusters ->
                !clusters.isEmpty()
            }
        )
        .set { clusteranalysis_input }

    SCRAMBLE_CLUSTERANALYSIS(
        clusteranalysis_input,
        fasta,
        []
    )

    ch_versions = ch_versions.mix(SCRAMBLE_CLUSTERANALYSIS.out.versions)

    SCRAMBLE_CLUSTERANALYSIS.out.meis_tab.combine(SCRAMBLE_CLUSTERANALYSIS.out.dels_tab, by:0)

    // Not finished yet, find a better dataset for this part

    emit:
    scramble_vcfs = Channel.empty()
    versions = ch_versions
}
