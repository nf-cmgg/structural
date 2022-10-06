//
// CNMOPS (Copy Number estimation by a Mixture Of PoissonS)
//

// Import modules
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE  } from '../../../modules/nf-core/gatk4/printsvevidence/main'

workflow CNMOPS {
    take:
        BAF_files               // channel: [optional]  [ baf_files ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file, index ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file, index ] => The split read evidence files created in GatherSampleEvidence

        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome

    main:

    ch_versions             = Channel.empty()

    // TODO add later (is mostly CNV focused)

    emit:
    versions                        = ch_versions
}
