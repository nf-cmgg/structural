//
// Batch Evidence Merging
//

// Import subworkflows
include { MAKE_BINCOV_MATRIX                            } from '../common-workflows/make-bincov-matrix'

// Import modules
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE_SR  } from '../../../modules/nf-core/modules/gatk4/printsvevidence/main'
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE_PE  } from '../../../modules/nf-core/modules/gatk4/printsvevidence/main'
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE_BAF } from '../../../modules/nf-core/modules/gatk4/printsvevidence/main'

workflow BATCH_EVIDENCE_MERGING {
    take:
        BAF_files               // channel: [optional]  [ baf_files ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file ] => The split read evidence files created in GatherSampleEvidence

        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome

    main:

    ch_versions         = Channel.empty()

    // Remove the last line from the merged file
    merged_pe_files     = PE_files.map({ meta, pe_file -> return pe_file }).collectFile(name: 'merged.pe.txt')
                                  .map({ merged_file -> [ [id:"pe"], merged_file, [] ]})
    // merged_sr_files     = SR_files.collectFile(name: 'merged.sr.txt.gz')
    // merged_baf_files    = BAF_files.collectFile(name: 'merged.baf.txt.gz')

    printsvevidence_input = merged_pe_files

    PRINTSVEVIDENCE_PE(
        printsvevidence_input,
        [],
        [],
        [],
        dict
    )

    emit:
    versions            = ch_versions
}