//
// Batch Evidence Merging
//

// Import modules
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE  } from '../../../modules/nf-core/gatk4/printsvevidence/main'
include { TABIX_TABIX as TABIX_EVIDENCE             } from '../../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_SORT                             } from '../../../modules/nf-core/bedtools/sort/main'

workflow BATCH_EVIDENCE_MERGING {
    take:
        BAF_files               // channel: [optional]  [ baf_files ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file, index ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file, index ] => The split read evidence files created in GatherSampleEvidence

        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome

    main:

    ch_versions             = Channel.empty()

    printsvevidence_input   = Channel.empty()

    pe_files = PE_files.map({ meta, file, index ->
                            new_meta = meta.clone()
                            new_meta.type = 'paired_end_evidence'
                            [ new_meta, file, index ]
                        })

    sr_files = SR_files.map({ meta, file, index ->
                            new_meta = meta.clone()
                            new_meta.type = 'split_read_evidence'
                            [ new_meta, file, index ]
                        })

    printsvevidence_input = pe_files.mix(sr_files)
                                    .map({ meta, file, index ->
                                        new_meta = [:]
                                        new_meta.id = meta.type
                                        [ new_meta, file, index ]
                                    })
                                    .groupTuple()

    // TODO add the BAF file support

    // merged_baf_files    = BAF_files.collectFile(name: 'merged.baf.txt.gz')

    // TODO find a way to fix the 'file not sorted' error in PRINTSVEVIDENCE

    // PRINTSVEVIDENCE(
    //     printsvevidence_input,
    //     [],
    //     [],
    //     [],
    //     dict
    // )

    // TODO SiteDepthtoBAF isn't in a released version of GATK => implement this later

    emit:
    // merged_evidence_files           = PRINTSVEVIDENCE.out.printed_evidence
    // merged_evidence_files_indices   = PRINTSVEVIDENCE.out.printed_evidence_index
    versions                        = ch_versions
}
