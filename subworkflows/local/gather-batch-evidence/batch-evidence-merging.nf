//
// Batch Evidence Merging
//

// Import modules
include { GATK4_PRINTSVEVIDENCE as PRINTSVEVIDENCE  } from '../../../modules/nf-core/gatk4/printsvevidence/main'
include { GATK4_SITEDEPTHTOBAF as SITEDEPTHTOBAF    } from '../../../modules/nf-core/gatk4/sitedepthtobaf/main'
include { TABIX_TABIX as TABIX_EVIDENCE             } from '../../../modules/nf-core/tabix/tabix/main'
include { BEDTOOLS_SORT                             } from '../../../modules/nf-core/bedtools/sort/main'

workflow BATCH_EVIDENCE_MERGING {
    take:
        BAF_files               // channel: [optional]  [ baf_files ] => The BAF files
        PE_files                // channel: [mandatory] [ meta, pe_file, index ] => The paired end evidence files created in GatherSampleEvidence
        SR_files                // channel: [mandatory] [ meta, sr_file, index ] => The split read evidence files created in GatherSampleEvidence

        SD_files                // channel: [optional]  [ meta, sd_file, index ] => The site depth evidence files created in GatherSampleEvidence
        allele_loci_vcf         // channel: [optional]  [ vcf, tbi ] => VCF of SNPs marking loci for allele count

        fasta                   // channel: [mandatory] [ fasta ] => The reference FASTA file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the FASTA reference file
        dict                    // channel: [mandatory] [ dict ] => The sequence dictionary of the reference genome

    main:

    ch_versions             = Channel.empty()

    printsvevidence_input   = Channel.empty()

    printsvevidence_input = printsvevidence_input.mix(
                            PE_files.map({ meta, file, index ->
                                new_meta = [:]
                                new_meta.id = 'paired_end_evidence'
                                [ new_meta, file, index ]
                            })
                        )

    printsvevidence_input = printsvevidence_input.mix(
                            SR_files.map({ meta, file, index ->
                                new_meta = [:]
                                new_meta.id = 'split_read_evidence'
                                [ new_meta, file, index ]
                            })
                        )

    if(BAF_files){

        printsvevidence_input = printsvevidence_input.mix(
                                BAF_files.map({ meta, file, index ->
                                    new_meta = [:]
                                    new_meta.id = 'baf_evidence'
                                    [ new_meta, file, index ]
                                })
                            )

    }

    // TODO find a way to fix the 'file not sorted' error in PRINTSVEVIDENCE

    // PRINTSVEVIDENCE(
    //     printsvevidence_input.groupTuple(),
    //     [],
    //     [],
    //     [],
    //     dict
    // )

    if(!SD_files){

        sitedepthtobaf_input = SD_files.map({ meta, file, index ->
                                    new_meta = [:]
                                    new_meta.id = 'site_depth_evidence'
                                    [ new_meta, file, index ]
                                }).groupTuple()

        SITEDEPTHTOBAF(
            sitedepthtobaf_input,
            allele_loci_vcf,
            fasta,
            fasta_fai,
            dict
        )

    }


    emit:
    // merged_evidence_files           = PRINTSVEVIDENCE.out.printed_evidence
    // merged_evidence_files_indices   = PRINTSVEVIDENCE.out.printed_evidence_index
    versions                        = ch_versions
}
