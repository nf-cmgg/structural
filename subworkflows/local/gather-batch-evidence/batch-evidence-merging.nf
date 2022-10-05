//
// Batch Evidence Merging
//

// Import modules
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE  } from '../../../modules/nf-core/modules/gatk4/printsvevidence/main'
include { TABIX_TABIX as TABIX                      } from '../../../modules/nf-core/modules/tabix/tabix/main'
include { BEDTOOLS_SORT                             } from '../../../modules/nf-core/modules/bedtools/sort/main'

workflow BATCH_EVIDENCE_MERGING {
    take:
        BAF_files               // channel: [optional]  [ baf_files ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file, index ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file, index ] => The split read evidence files created in GatherSampleEvidence

        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome

    main:

    ch_versions         = Channel.empty()

    printsvevidence_input = Channel.empty()

    printsvevidence_input  = printsvevidence_input.mix(PE_files.map({ meta, file, index ->
                                        [ [type:'pe'], file, index ]
                                    })
                                    .groupTuple()
                                    .map({ meta, files, indices -> [ files, indices ]})
                                )
    printsvevidence_input  = printsvevidence_input.mix(SR_files.map({ meta, file, index ->
                                        [ [type:'sr'], file, index ]
                                    })
                                    .groupTuple()
                                    .map({ meta, files, indices -> [ files, indices ]})
                                )

    // merged_baf_files    = BAF_files.collectFile(name: 'merged.baf.txt.gz')

    PRINTSVEVIDENCE(
        printsvevidence_input,
        [],
        [],
        [],
        dict
    ).printed_evidence.view()

    emit:
    versions            = ch_versions
}
