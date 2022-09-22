//
// Gather Sample Evidence
//

// Import subworkflows
include { RUN_MANTA                                     } from './run-manta'
include { RUN_DELLY                                     } from './run-delly'
include { RUN_WHAMG                                     } from './run-whamg'
include { RUN_SCRAMBLE                                  } from './run-scramble'
include { GATHER_SAMPLE_EVIDENCE_METRICS                } from './gather-metrics'

// Import modules
include { GATK4_COLLECTREADCOUNTS as COLLECTREADCOUNTS  } from '../../../modules/nf-core/modules/gatk4/collectreadcounts/main'
include { GATK4_COLLECTSVEVIDENCE as COLLECTSVEVIDENCE  } from '../../../modules/nf-core/modules/gatk4/collectsvevidence/main'
include { TABIX_TABIX                                   } from '../../../modules/nf-core/modules/tabix/tabix/main'

workflow GATHER_SAMPLE_EVIDENCE {
    take:
        crams                   // channel: [mandatory] [ meta, cram, crai, bed ] => The aligned CRAMs per sample with the regions they should be called on
        beds                    // channel: [optional]  [ meta, bed, bed_gz, bed_gz_tbi ] => A channel containing the normal BED, the bgzipped BED and its index file
        allele_loci_vcf         // channel: [optional]  [ vcf ] => A channel containing the VCF and its index for counting the alleles
        fasta                   // channel: [mandatory] [ fasta ] => The fasta reference file
        fasta_fai               // channel: [mandatory] [ fasta_fai ] => The index of the fasta reference file
        dict                    // channel: [mandatory] [ dict ] => The dictionary of the fasta reference file

    main:

    callers = params.callers.split(",")

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    ch_metrics  = Channel.empty()
    called_vcfs = Channel.empty()

    collectreadcounts_input = crams.combine(
        beds.map({meta, bed, bed_gz, bed_gz_tbi -> [meta, bed]}), by:0
    )

    COLLECTREADCOUNTS(
        collectreadcounts_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COLLECTREADCOUNTS.out.versions)

    if("manta" in callers){
        RUN_MANTA(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(RUN_MANTA.out.manta_vcfs)
        ch_versions = ch_versions.mix(RUN_MANTA.out.versions)
    }

    if("delly" in callers){
        RUN_DELLY(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(RUN_DELLY.out.delly_vcfs)
        ch_versions = ch_versions.mix(RUN_DELLY.out.versions)
    }

    // Whamg needs some reheadering (like done in https://github.com/broadinstitute/gatk-sv/blob/90e3e9a221bdfe7ab2cfedeffb704bc6f0e99aa9/wdl/Whamg.wdl#L209)
    if("whamg" in callers){
        RUN_WHAMG(
            crams,
            beds,
            fasta,
            fasta_fai
        )

        called_vcfs = called_vcfs.mix(RUN_WHAMG.out.whamg_vcfs)
        ch_versions = ch_versions.mix(RUN_WHAMG.out.versions)
    }

    // Scramble is unfinished. It needs a lot of improvements if we were to add it

    // if("scramble" in callers){
    //     RUN_SCRAMBLE(
    //         crams,
    //         beds,
    //         fasta
    //     )

    //     called_vcfs = called_vcfs.mix(RUN_SCRAMBLE.out.scramble_vcfs)
    //     ch_versions = ch_versions.mix(RUN_SCRAMBLE.out.versions)
    // }

    if(allele_loci_vcf){
        TABIX_TABIX(
            [ [], allele_loci_vcf ]
        )

        collectsvevidence_input = crams.combine(TABIX_TABIX.out.tbi)
                                    .map({ meta, cram, crai, tbi ->
                                        [ meta, cram, crai, allele_loci_vcf, tbi ]
                                    })
    } else {
        collectsvevidence_input = crams.map({ meta, cram, crai -> [ meta, cram, crai, [], [] ] })
    }


    COLLECTSVEVIDENCE(
        collectsvevidence_input,
        fasta,
        fasta_fai,
        dict
    )

    ch_versions = ch_versions.mix(COLLECTSVEVIDENCE.out.versions)

    if(params.run_module_metrics) {
        GATHER_SAMPLE_EVIDENCE_METRICS(
            called_vcfs,
            COLLECTSVEVIDENCE.out.split_read_evidence,
            COLLECTSVEVIDENCE.out.paired_end_evidence,
            COLLECTSVEVIDENCE.out.allele_counts,
            fasta_fai
        )

        ch_metrics  = ch_reports.mix(GATHER_SAMPLE_EVIDENCE_METRICS.out.metrics)
        ch_versions = ch_versions.mix(GATHER_SAMPLE_EVIDENCE_METRICS.out.versions)
    }

    emit:
    vcfs                = called_vcfs

    coverage_counts     = COLLECTREADCOUNTS.out.tsv

    split_reads         = COLLECTSVEVIDENCE.out.split_read_evidence
    split_reads_index   = COLLECTSVEVIDENCE.out.split_read_evidence_index
    read_pairs          = COLLECTSVEVIDENCE.out.paired_end_evidence
    read_pairs_index    = COLLECTSVEVIDENCE.out.paired_end_evidence_index
    allele_counts       = COLLECTSVEVIDENCE.out.allele_counts
    allele_counts_index = COLLECTSVEVIDENCE.out.allele_counts_index

    versions            = ch_versions
    metrics             = ch_metrics
    reports             = ch_reports
}